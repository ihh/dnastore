# DNAstore
Software for error-tolerant coding of information into DNA sequences using finite-state transducers.

For license, please see LICENSE.txt

To build:

    make dnastore
    make test

I have built using Apple LLVM version 7.3.0 (clang-703.0.31).
You will also need Boost: http://www.boost.org/

Description of method is in doc/ directory:

    cd doc
    make
    open trans.ps

The tests include a few examples.

To generate a code transducer with 4 nucleotides of context (this will avoid all homopolymer and dinucleotide repeats):

    bin/dnastore -l 4

To encode the string "Hello World!" in DNA using this transducer:

    bin/dnastore -l 4 -E "Hello World!" > HelloWorld.fasta

To decode this DNA:

    bin/dnastore -l 4 -d HelloWorld.fasta

To construct a composite machine using the MIXRADAR6 block code for mixed-radix conversion:

    bin/dnastore -l 4 --compose-machine data/flusher.json  --compose-machine data/mixradar6.json --save-machine mixradar6-dnastore4.json

(The <code>flusher.json</code> file describes an outer transducer that automatically flushes MIXRADAR6 whenever the end of the file is reached.)

To encode and decode using this transducer:

    bin/dnastore --load-machine mixradar6-dnastore4.json -E "Hello World!" >HelloWorld46.fasta
    bin/dnastore --load-machine mixradar6-dnastore4.json -d HelloWorld46.fasta

To use Viterbi decoding instead of exact decoding:

    bin/dnastore --load-machine mixradar6-dnastore4.json -V HelloWorld46.fasta --error-global

To create a 64-bit watermark synchronization code with 1 watermark bit per signal bit:

    bin/dnastore -l 4 --compose-machine data/water64.1.json --save-machine watmark64-dnastore4.json

Strings encoded using this code must be a multiple of 64 bits in length:

    bin/dnastore --load-machine watmark64-dnastore4.json -E "Hello World! (192 bits.)" >hw64.fa
    bin/dnastore --load-machine watmark64-dnastore4.json -d hw64.fa

For a list of more options:

    bin/dnastore -h

You can also look at the script <code>doc/errdecode.pl</code> for examples of usage (this script does a lot of benchmarking).
