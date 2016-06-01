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
codes: data/hamming74.json data/mixradar2.json data/mixradar6.json data/sync16.json data/sync64.json data/sync128.json

data/hamming74.json: bin/hamming74.pl
	perl $< --json >$@

data/mixradar%.json: bin/mixradar.pl
	perl $< --flush --json --verbose $* .001 >$@

data/sync%.json: bin/syncer.pl
	perl $< $* >$@

# Tests
TEST = t/testexpect.pl
NOSUBS = --error-sub-prob 0
NODUPS = --error-dup-prob 0
NODELS = --error-del-open 0
GLOBAL = --error-global
NOERRS = $(NOSUBS) $(NODUPS) $(NODELS) $(GLOBAL)
ONLYDUPS = $(NOSUBS) $(NODELS) $(GLOBAL)

test: testpattern testmachine testencode testdecode testviterbi testcompose testham testsync testsyncham testfit

testpattern: bin/testpattern
	$<

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

testfit: $(MAIN)
	@$(TEST) bin/$(MAIN) -v0 --fit-error data/tiny.stk data/tiny.params.json
	@$(TEST) bin/$(MAIN) -v0 --fit-error data/test.stk data/test.params.json

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
