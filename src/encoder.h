#ifndef ENCODER_INCLUDED
#define ENCODER_INCLUDED

#include "trans.h"

struct Encoder {
  const Machine& machine;
  State current;
  bool msb0;  // set this to encode MSB first, instead of LSB first
  Encoder (const Machine& machine)
    : machine(machine),
      current(0),
      msb0(false)
  { }

  template<class Writer>
  void write (char outc, Writer& outs) {
    if (outc)
      (void) outs.write (&outc, 1);
  }

  template<class Writer>
  void encodeSymbol (char sym, Writer& outs) {
    LogThisAt(8,"Encoding " << sym << endl);
    while (!machine.state[current].isInput()) {
      const MachineTransition& tn = machine.state[current].next();
      write (tn.out, outs);
      LogThisAt(9,"Transition " << machine.stateName(current)
		<< " -> " << machine.stateName(tn.dest)
		<< ": output " << tn.out
		<< endl);
      current = tn.dest;
    }
    const MachineTransition* t = machine.state[current].transFor (sym);
    Assert (t != NULL, "Couldn't encode symbol %c in state %s", sym, machine.stateName(current).c_str());
    while (true) {
      write (t->out, outs);
      LogThisAt(9,"Transition " << machine.stateName(current)
		<< " -> " << machine.stateName(t->dest) << ": "
		<< (t->in ? (string("input ") + t->in + ", ") : string())
		<< "output " << t->out
		<< endl);
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

  template<class Writer, typename CharType>
  void encodeByte (CharType byte, Writer& outs) {
    LogThisAt(7,"Encoding '" << (char)byte << "' (\\x" << hex << (int)byte << ")" << endl);
    if (msb0)
      for (int n = 7; n >= 0; --n)
	encodeBit (byte & (1 << n), outs);
    else
      for (int n = 0; n <= 7; ++n)
	encodeBit (byte & (1 << n), outs);
  }

  template<class Writer>
  void encodeStream (istream& in, Writer& outs) {
    istreambuf_iterator<char> iter(in), iterEnd;
    while (iter != iterEnd) {
      encodeByte (*iter, outs);
      ++iter;
    }
  }

  template<class Writer>
  void encodeString (const string& s, Writer& outs) {
    for (auto c: s)
      encodeByte (c, outs);
  }

  template<class Writer>
  void encodeSymbolString (const string& s, Writer& outs) {
    for (auto c: s)
      encodeSymbol (c, outs);
  }
};

struct FastaWriter {
  ostream& outs;
  size_t col, colsPerLine;
  
  FastaWriter (ostream& outs, const char* seqname = "SEQ")
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
      if (++col >= colsPerLine) {
	outs << endl;
	col = 0;
      }
    }
  }
};

#endif /* ENCODER_INCLUDED */
