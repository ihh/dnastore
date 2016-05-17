#ifndef ENCODER_INCLUDED
#define ENCODER_INCLUDED

#include "trans.h"

struct Encoder {
  const Machine& machine;
  State current;
  Encoder (const Machine&)
    : machine(machine),
      current(0)
  { }

  template<class Writer>
  void write (bool bit, Writer& outs) {
    const MachineTransition* t = machine.state[current].transFor (bit ? '1' : '0');
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
