#ifndef ENCODER_INCLUDED
#define ENCODER_INCLUDED

#include "trans.h"

struct Encoder {
  const Machine& machine;
  State current;
  bool bigEndian, msb0;
  Encoder (const Machine&)
    : machine(machine),
      current(0),
      bigEndian(false),
      msb0(false)
  { }

  template<class Writer>
  void encodeSymbol (char sym, Writer& outs) {
    const MachineTransition* t = machine.state[current].transFor (sym);
    Assert (t != NULL, "Couldn't encode bit");
    while (true) {
      if (t->out)
	(void) outs.write (&t->out, 1);
      current = t->dest;
      const MachineState& ms = machine.state[current];
      if (ms.isInput())
	break;
      t = &ms.next();
    }
  }

  template<class Writer>
  void encodeBit (bool bit, Writer& outs) {
    encodeSymbol (bit ? '1' : '0', outs);
  }

  template<class Writer>
  void encodeMeta (ControlIndex control, Writer& outs) {
    encodeSymbol (Machine::controlChar (control), outs);
  }

  template<class Writer, CharType>
  void encodeByte (CharType byte, Writer& outs) {
    if (msb0)
      for (int n = 7; n >= 0; --n)
	encodeBit (byte & (1 << n), outs);
    else
      for (int n = 0; n <= 7; ++n)
	encodeBit (byte & (1 << n), outs);
  }

  template<class Writer, WordType>
  void encodeWord (WordType word, Writer& outs) {
    if (bigEndian)
      for (int n = 24; n >= 0; n -= 8)
	encode ((word >> n) & 0xff, outs);
    else
      for (int n = 0; n <= 24; n += 8)
	encode ((word >> n) & 0xff, outs);
  }
};

struct FastaWriter {
  ostream& outs;
  size_t col, maxCols;
  
  FastaWriter (ostream& outs, char* seqname = "SEQ")
    : outs(outs),
      col(0),
      colsPerLine(50)
  {
    outs << ">" << seqname << endl;
  }

  ~FastaWriter() {
    if (col > 0)
      outs << endl;
  }
  
  void write (char* buf, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      outs << buf[i];
      if (++col >= maxCols) {
	outs << endl;
	col = 0;
      }
    }
  }
};

#endif /* ENCODER_INCLUDED */
