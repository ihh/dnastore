#ifndef DECODER_INCLUDED
#define DECODER_INCLUDED

#include "trans.h"

template<class Writer>
struct Decoder {
  typedef map<State,string> StateString;
  typedef typename StateString::iterator StateStringIter;
  
  const Machine& machine;
  Writer& outs;
  StateString current;

  Decoder (const Machine& machine, Writer& outs)
    : machine(machine),
      outs(outs)
  {
    current[machine.startState()] = string();
    expand (false);
  }

  ~Decoder() {
    expand (true);
    close();
  }

  void close() {
    vguard<StateStringIter> ssIter;
    for (StateStringIter ss = current.begin(); ss != current.end(); ++ss)
      if (machine.state[ss->first].isEnd())
	ssIter.push_back (ss);
    if (ssIter.size() == 1)
      flush (ssIter.front());
    else if (ssIter.size() > 1) {
      Warn ("Decoder unresolved: %u possible end state(s)", ssIter.size());
      for (auto ss: ssIter)
	Warn ("State %s: input queue %s", machine.state[ss->first].name.c_str(), ss->second.c_str());
    } else if (current.size() > 1) {
      Warn ("Decoder unresolved: %u possible state(s)", ssIter.size());
      for (const auto& ss: current)
	Warn ("State %s: input queue %s", machine.state[ss.first].name.c_str(), ss.second.c_str());
    }
    current.clear();
  }
  
  void expand (bool includeEndStates) {
    StateString next;
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
	  else if (includeEndStates || !machine.state[t.dest].isEnd()) {
	    const string tStr = t.in ? (str + t.in) : str;
	    if (next.count (t.dest))
	      Assert (next.at(t.dest) == tStr, "Decoder error");
	    else {
	      LogThisAt(9,"Transition " << machine.state[state].name
			<< " -> " << machine.state[t.dest].name
			<< (tStr.empty() ? string() : (string(": input queue ") + tStr + ", "))
			<< endl);
	      next[t.dest] = tStr;
	      foundNew = true;
	    }
	  }
	if (hasOutput)
	  next[state] = str;
      }
      current.swap (next);
    } while (foundNew);
  }

  void flush (StateStringIter ss) {
    string& str = ss->second;
    if (str.size()) {
      LogThisAt(9,"Flushing input queue: " << str << endl);
      char* buf = (char*) malloc (sizeof(char) * (str.size() + 1));
      strcpy (buf, str.c_str());
      (void) outs.write (buf, str.size());
      free (buf);
      str.clear();
    }
  }
  
  void decodeBase (char base) {
    base = toupper(base);
    LogThisAt(8,"Decoding " << base << endl);
    StateString next;
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
	  LogThisAt(9,"Transition " << machine.state[state].name
		    << " -> " << machine.state[nextState].name
		    << ": "
		    << (nextStr.empty() ? string() : (string("input queue ") + nextStr + ", "))
		    << "output " << t.out
		    << endl);
	}
    }
    Assert (!next.empty(), "No inputs consistent with given output sequence");
    current.swap (next);
    expand (false);
    if (current.size() == 1) {
      auto iter = current.begin();
      const MachineState& ms = machine.state[iter->first];
      if (ms.isEnd() || ms.isInput())
	flush (iter);
    }
  }

  void decodeString (const string& seq) {
    for (char c: seq)
      decodeBase (c);
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
    if (!outbuf.empty()) {
      if (!msb0)
	reverse (outbuf.begin(), outbuf.end());
      Warn ("%u bits (%s) remain on output", outbuf.size(), to_string_join(outbuf,"").c_str());
    }
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

