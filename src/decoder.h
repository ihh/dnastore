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

  ~Decoder() {
    bool unresolved = false;
    for (const auto& ss: current)
      if (!ss.second.empty()) {
	unresolved = true;
	break;
      }
    if (unresolved) {
      Warn ("Decoder unresolved: %d state(s) remaining with symbols on input queue", current.size());
      for (const auto& ss: current)
	Warn ("State %s: input queue %s", machine.stateName(ss.first).c_str(), ss.second.c_str());
    }
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
      const StateType type = machine.state[iter->first].type;
      if (type == ControlState || type == CodeState) {
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
      Warn ("%d bits (%s) remain on output", outbuf.size(), to_string_join(outbuf,"").c_str());
  }

  void flush() {
    unsigned char c = 0;
    for (size_t n = 0; n < outbuf.size(); ++n)
      if (outbuf[n])
	c = c | (1 << (msb0 ? (7-n) : n));
    LogThisAt(7,"Decoding '" << (char)c << "' (\\x" << hex << (int)c << ")" << endl);
    outs << c;
    outbuf.clear();
  }
  
  void write (char* buf, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      const char c = buf[i];
      if (c == '0' || c == '1') {
	outbuf.push_back (c == '1');
	if (outbuf.size() == 8)
	  flush();
      } else {
	const ControlIndex idx = Machine::controlIndex(c);
	if (idx >= 0)
	  Warn("Ignoring control character #%d (%c)",idx,c);
	else
	  Warn("Ignoring unknown character '%c' (\\x%.2x)",c,c);
      }
    }
  }
};

#endif /* DECODER_INCLUDED */

