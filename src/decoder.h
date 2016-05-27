#ifndef DECODER_INCLUDED
#define DECODER_INCLUDED

#include "trans.h"

template<class Writer>
struct Decoder {
  typedef map<State,deque<char> > StateString;
  typedef typename StateString::iterator StateStringIter;
  
  const Machine& machine;
  Writer& outs;
  StateString current;

  Decoder (const Machine& machine, Writer& outs)
    : machine(machine),
      outs(outs)
  {
    current[machine.startState()] = deque<char>();
    expand();
  }

  ~Decoder() {
    close();
  }

  void close() {
    if (current.size()) {
      expand();
      vguard<StateStringIter> ssIter;
      for (StateStringIter ss = current.begin(); ss != current.end(); ++ss)
	if (machine.state[ss->first].isEnd())
	  ssIter.push_back (ss);
      if (ssIter.size() == 1)
	flush (ssIter.front());
      else if (ssIter.size() > 1) {
	Warn ("Decoder unresolved: %u possible end state(s)", ssIter.size());
	for (auto ss: ssIter)
	  Warn ("State %s: input queue %s", machine.state[ss->first].name.c_str(), to_string_join(ss->second,"").c_str());
      } else if (current.size() > 1) {
	Warn ("Decoder unresolved: %u possible state(s)", ssIter.size());
	for (const auto& ss: current)
	  Warn ("State %s: input queue %s", machine.state[ss.first].name.c_str(), to_string_join(ss.second,"").c_str());
      }
      current.clear();
    }
  }
  
  void expand() {
    StateString next;
    bool foundNew;
    do {
      foundNew = false;
      for (const auto& ss: current) {
	const State state = ss.first;
	const auto& str = ss.second;
	const MachineState& ms = machine.state[state];
	if (ms.isEnd() || ms.emitsOutput())
	  next[state] = str;
      }
      for (const auto& ss: current) {
	const State state = ss.first;
	const auto& str = ss.second;
	const MachineState& ms = machine.state[state];
	for (const auto& t: ms.trans)
	  if (isUsable(t) && !t.out) {
	    auto nextStr = str;
	    if (!t.inputEmpty() && !t.isEOF())
	      nextStr.push_back (t.in);
	    if (next.count (t.dest))
	      Assert (next.at(t.dest) == nextStr, "Decoder error");
	    else {
	      next[t.dest] = nextStr;
	      LogThisAt(9,"Transition " << ms.name
			<< " -> " << machine.state[t.dest].name
			<< (nextStr.empty() ? string() : (string(": input queue ") + to_string_join(nextStr,"")))
			<< endl);
	      foundNew = true;
	    }
	  }
      }
      current.swap (next);
      next.clear();
    } while (foundNew);
  }

  void write (const char* s, size_t len) {
    char* buf = (char*) malloc (sizeof(char) * (len + 1));
    for (size_t n = 0; n < len; ++n)
      buf[n] = s[n];
    buf[len] = '\0';
    (void) outs.write (buf, len);
    free (buf);
  }
  
  void flush (StateStringIter ss) {
    const string str (ss->second.begin(), ss->second.end());
    if (str.size()) {
      LogThisAt(9,"Flushing input queue: " << str << endl);
      write (str.c_str(), str.size());
      ss->second.clear();
    }
  }

  static bool isUsable (const MachineTransition& t) {
    return t.in == MachineNull || t.in == MachineBit0 || t.in == MachineBit1 || t.in == MachineEOF;
  }
  
  void decodeBase (char base) {
    base = toupper(base);
    LogThisAt(8,"Decoding " << base << endl);
    StateString next;
    for (const auto& ss: current) {
      const State state = ss.first;
      const auto& str = ss.second;
      for (const auto& t: machine.state[state].trans)
	if (isUsable(t) && t.out == base) {
	  const State nextState = t.dest;
	  auto nextStr = str;
	  if (!t.inputEmpty())
	    nextStr.push_back (t.in);
	  Assert (!next.count(nextState) || next.at(nextState) == nextStr,
		  "Multiple outputs decode to single state");
	  next[nextState] = nextStr;
	  LogThisAt(9,"Transition " << machine.state[state].name
		    << " -> " << machine.state[nextState].name
		    << ": "
		    << (nextStr.empty() ? string() : (string("input queue ") + to_string_join(nextStr,"") + ", "))
		    << "output " << t.out
		    << endl);
	}
    }
    Assert (!next.empty(), "No inputs consistent with given output sequence");
    current.swap (next);
    expand();
    if (current.size() == 1) {
      auto iter = current.begin();
      const MachineState& ms = machine.state[iter->first];
      if (ms.exitsWithInput())
	flush (iter);
    } else
      shiftResolvedSymbols();
  }

  void shiftResolvedSymbols() {
    while (true) {
      bool foundQueue = false, queueNonempty, firstCharSame;
      char firstChar;
      for (const auto& ss: current) {
	if (!foundQueue) {
	  if ((queueNonempty = !ss.second.empty()))
	    firstChar = ss.second[0];
	  foundQueue = firstCharSame = true;
	} else if (queueNonempty
		   && (ss.second.empty()
		       || firstChar != ss.second[0])) {
	  firstCharSame = false;
	  break;
	}
      }
      if (foundQueue && queueNonempty && firstCharSame) {
	LogThisAt(9,"All input queues have '" << Machine::charToString(firstChar) << "' as first symbol; shifting" << endl);
	write (&firstChar, 1);
	for (auto& ss: current)
	  ss.second.erase (ss.second.begin());
      } else
	break;
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
      Warn ("%s (%s) remaining on output", plural(outbuf.size(),"bit").c_str(), to_string_join(outbuf,"").c_str());
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
      if (c == MachineBit0 || c == MachineBit1) {
	outbuf.push_back (c == MachineBit1);
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

