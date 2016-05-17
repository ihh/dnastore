#ifndef DECODER_INCLUDED
#define DECODER_INCLUDED

#include "trans.h"

struct Decoder {
  const Machine& machine;
  map<State,string> current;
  Decoder (const Machine& machine)
    : machine(machine)
  {
    current[machine.startState()] = string();
    expand();
  }

  void expand() {
    map<State,string> next;
    bool foundNew;
    do {
      foundNew = false;
      for (const auto& ss: current) {
	const State state = ss.first;
	const string& str = ss.second;
	bool hasOutput = false;
	for (const auto& t: machine.state[state].trans)
	  if (t.out)
	    hasOutput = true;
	  else {
	    next[t.dest] = t.in ? (str + t.in) : str;
	    foundNew = true;
	  }
	if (hasOutput)
	  next[state] = str;
      }
      current.swap (next);
    } while (foundNew);
  }

  template<class Writer>
  void decodeBase (char base, Writer& outs) {
    base = toupper(base);
    LogThisAt(8,"Decoding " << base << endl);
    map<State,string> next;
    for (const auto& ss: current) {
      const State state = ss.first;
      const string& str = ss.second;
      for (const auto& t: machine.state[state].trans)
	if (t.out == base) {
	  const State nextState = t.dest;
	  const string nextStr = t.in ? (str + t.in) : str;
	  Assert (!next.count(nextState) || next.at(nextState) == nextStr,
		  "Multiple outputs decode to single state");
	  next[nextState] = nextStr;
	  LogThisAt(9,"Transition " << machine.stateName(state) << " -> " << machine.stateName(nextState)
		    << ": "
		    << (nextStr.empty() ? string() : (string("input queue ") + nextStr + ", "))
		    << "output " << t.out
		    << endl);
	}
    }
    Assert (!next.empty(), "No inputs consistent with given output sequence");
    current.swap (next);
    expand();
    if (current.size() == 1) {
      auto iter = current.begin();
      string& str = iter->second;
      if (str.size()) {
	LogThisAt(9,"Flushing input queue: " << str << endl);
	char* buf = (char*) malloc (sizeof(char) * (str.size() + 1));
	strcpy (buf, str.c_str());
	(void) outs.write (buf, str.size());
	free (buf);
	str.clear();
      }
    }
  }

  template<class Writer>
  void decodeString (const string& seq, Writer& outs) {
    for (char c: seq)
      decodeBase (c, outs);
  }
};

struct BinaryWriter {
  ostream& outs;
  bool msb0;
  vguard<bool> outbuf;

  BinaryWriter (ostream& outs)
    : outs(outs),
      msb0(false)
  { }

  ~BinaryWriter() {
    if (!outbuf.empty())
      flush();
  }

  void flush() {
    unsigned char c = 0;
    for (size_t n = 0; n < outbuf.size(); ++n)
      if (outbuf[n])
	c = c | (1 << (msb0 ? (7-n) : n));
    outs << c;
    outbuf.clear();
  }
  
  void write (char* buf, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      const char c = buf[n];
      if (c == '0' || c == '1') {
	outbuf.push_back (c == '1');
	if (outbuf.size() == 8)
	  flush();
      }
    }
  }
};

#endif /* DECODER_INCLUDED */

