.SECONDARY:

BOOSTPREFIX = /usr
ifeq (,$(wildcard $(BOOSTPREFIX)/include/boost/regex.h))
BOOSTPREFIX = /usr/local
ifeq (,$(wildcard $(BOOSTPREFIX)/include/boost/regex.h))
BOOSTPREFIX =
endif
endif

BOOSTFLAGS =
BOOSTLIBS =
ifneq (,$(BOOSTPREFIX))
BOOSTFLAGS := -I$(BOOSTPREFIX)/include
BOOSTLIBS := -L$(BOOSTPREFIX)/lib -lboost_regex -lboost_program_options
endif

# install dir
PREFIX = /usr/local

# other flags
ifneq (,$(findstring debug,$(MAKECMDGOALS)))
CPPFLAGS = -std=c++11 -g -DUSE_VECTOR_GUARDS -DDEBUG $(BOOSTFLAGS)
else
CPPFLAGS = -std=c++11 -g -O3 $(BOOSTFLAGS)
endif
LIBFLAGS = -lstdc++ -lz $(BOOSTLIBS)

CPPFILES = $(wildcard src/*.cpp)
OBJFILES = $(subst src/,obj/,$(subst .cpp,.o,$(CPPFILES)))

# try clang++, fall back to g++
CPP = clang++
ifeq (, $(shell which $(CPP)))
CPP = g++
endif

# pwd
PWD = $(shell pwd)

# /bin/sh
SH = /bin/sh

# Targets

MAIN = dnastore

all: $(MAIN)

# Main build rules
bin/%: $(OBJFILES) obj/%.o
	@test -e bin || mkdir bin
	$(CPP) $(LIBFLAGS) -o $@ obj/$*.o $(OBJFILES)

obj/%.o: src/%.cpp
	@test -e obj || mkdir obj
	$(CPP) $(CPPFLAGS) -c -o $@ $<

obj/%.o: t/%.cpp
	@test -e obj || mkdir obj
	$(CPP) $(CPPFLAGS) -c -o $@ $<

$(MAIN): bin/$(MAIN)

clean:
	rm -rf bin/$(MAIN) obj/*

debug: all

# Codes
codes: hamming_codes sync_codes mixradar_codes watermark_codes

# HAMMING
hamming_codes: data/hamming74.json

data/hamming74.json: bin/hamming74.pl
	perl $< --json >$@

# MIXRADAR
mixradar_codes: data/mixradar2.json data/mixradar6.json

data/mixradar%.json: bin/mixradar.pl
	perl $< --flush --json --verbose $* .001 >$@

# SYNC
sync_codes: data/sync16.json data/sync64.json data/sync128.json

data/sync%.json: bin/syncer.pl
	perl $< $* >$@

# WMARK
WATER128 = $(addprefix data/water,$(addsuffix .json,128 128.16 128.16b 128n 128a2 128.1))
WATER64 = $(addprefix data/water,$(addsuffix .json,64 64.16 64.16b 64n 64a2 64.1))
watermark_codes: $(WATER128) $(WATER64)

data/water%.1.json:
	bin/wmark.pl $* -sub 1 >$@

data/water%.16.json:
	bin/wmark.pl $* -sub 16 >$@

data/water%.16b.json:
	bin/wmark.pl $* -sub 16 -word >$@

data/water%n.json:
	bin/wmark.pl $* -nomixwater >$@

data/water%a2.json:
	bin/wmark.pl $* -copies 2 >$@

data/water%.json:
	bin/wmark.pl $* >$@

# Tests
TEST = t/testexpect.pl
NOSUBS = --error-sub-prob 0
NODUPS = --error-dup-prob 0
NODELS = --error-del-open 0
GLOBAL = --error-global
NOERRS = $(NOSUBS) $(NODUPS) $(NODELS) $(GLOBAL)
ONLYDUPS = $(NOSUBS) $(NODELS) $(GLOBAL)

test: testpattern testdist testmachine testencode testdecode testviterbi testcompose testham testsync testsyncham testcount testfit

testpattern: bin/testpattern
	$<

testdist: bin/editdist
	@$(TEST) $< ABCDEF ADEF 2
	@$(TEST) $< '""' '""' 0
	@$(TEST) $< '""' 'ABC' 3

testmachine: $(MAIN)
	@$(TEST) bin/$(MAIN) -v0 --length 4 --controls 4 --save-machine - data/l4c4.json
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/l4c4.json --save-machine - data/l4c4.json

testencode: $(MAIN)
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/l4c4.json --encode-file data/hello.txt data/hello.fa
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/l4c4.json --raw --encode-string HELLO data/hello.dna

testdecode: $(MAIN)
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/l4c4.json --decode-file data/hello.fa data/hello.txt
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/l4c4.json --decode-string `cat data/hello.dna` data/hello.txt
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/l4c4.json --decode-bits `cat data/hello.dna` data/hello.padded.bits

testviterbi: $(MAIN)
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/l4c4.json --decode-viterbi data/hello.fa $(NOERRS) --raw data/hello.padded.bits
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/l4c4.json --decode-viterbi data/hello.dup.fa $(ONLYDUPS) --raw data/hello.padded.bits

testcompose: $(MAIN) data/mixradar2.json
	@$(TEST) bin/$(MAIN) -v0 --compose-machine data/mixradar2.json --load-machine data/l4c4.json --save-machine - data/mr2l4c4.json
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/mr2l4c4.json --encode-file data/hello.txt data/hello.mr2.fa
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/mr2l4c4.json --decode-file data/hello.mr2.fa data/hello.txt
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/mr2l4c4.json --decode-viterbi data/hello.mr2.fa $(NOERRS) --raw data/hello.exact.bits

testcount: $(MAIN)
	@$(TEST) bin/$(MAIN) -v0 -l6 --error-sub-prob 1e-9 --error-dup-prob 1e-9 --error-del-open 1e-9 --error-counts data/dup.stk data/dup.counts.json
	@$(TEST) bin/$(MAIN) -v0 -l6 --error-sub-prob 1e-9 --error-dup-prob 1e-9 --error-del-open 1e-9 --error-counts data/dup.sub.stk data/dup.sub.counts.json
	@$(TEST) bin/$(MAIN) -v0 -l6 --error-sub-prob 1e-9 --error-dup-prob 1e-9 --error-del-open 1e-9 --error-counts data/dup.sub.misaligned.stk data/dup.sub.counts.misaligned.json

testfit: $(MAIN)
	@$(TEST) bin/$(MAIN) -v0 --fit-error data/tiny.stk --strict-guides data/tiny.params.json
	@$(TEST) bin/$(MAIN) -v0 --fit-error data/test.stk --strict-guides data/test.params.json

testham: $(MAIN) data/hamming74.json
	@$(TEST) bin/$(MAIN) -v0 --compose-machine data/hamming74.json --load-machine data/l4c4.json --save-machine - data/h74l4c4.json
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/h74l4c4.json --encode-file data/hello.txt data/hello.h74.fa
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/h74l4c4.json --decode-file data/hello.h74.fa data/hello.txt
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/h74l4c4.json --decode-viterbi data/hello.h74.fa $(NOERRS) --raw data/hello.exact.bits
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/h74l4c4.json --decode-viterbi data/hello.h74.fa --raw data/hello.exact.bits
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/h74l4c4.json --decode-viterbi data/hello.h74.sub.fa --raw data/hello.exact.bits

testsync: $(MAIN) data/sync16.json
	@$(TEST) bin/$(MAIN) -v0 --compose-machine data/sync16.json --compose-machine data/flusher.json --compose-machine data/mixradar2.json --load-machine data/l4c4.json --save-machine - data/s16mr2l4c4.json
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/s16mr2l4c4.json --encode-file data/hello.txt data/hello.s16mr2.fa
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/s16mr2l4c4.json --decode-file data/hello.s16mr2.fa data/hello.txt
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/s16mr2l4c4.json --decode-viterbi data/hello.s16mr2.fa $(NOERRS) --raw data/hello.exact.bits
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/s16mr2l4c4.json --decode-viterbi data/hello.s16mr2.fa --raw data/hello.exact.bits

testsyncham: $(MAIN) data/sync16.json data/hamming74.json
	@$(TEST) bin/$(MAIN) -v0 --compose-machine data/sync16.json --compose-machine data/flusher.json --compose-machine data/hamming74.json --load-machine data/l4c4.json --save-machine - data/s16h74l4c4.json
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/s16h74l4c4.json --encode-file data/hello.txt data/hello.s16h74.fa
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/s16h74l4c4.json --decode-file data/hello.s16h74.fa data/hello.txt
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/s16h74l4c4.json --decode-viterbi data/hello.s16h74.fa $(NOERRS) --raw data/hello.exact.bits
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/s16h74l4c4.json --decode-viterbi data/hello.s16h74.fa --raw data/hello.exact.bits
	@$(TEST) bin/$(MAIN) -v0 --load-machine data/s16h74l4c4.json --decode-viterbi data/hello.s16h74.del.fa --raw data/hello.exact.bits
