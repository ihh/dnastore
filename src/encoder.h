#ifndef ENCODER_INCLUDED
#define ENCODER_INCLUDED

#include "trans.h"

template<class Writer>
struct Encoder {
  typedef map<State,deque<InputSymbol> > StateString;
  typedef typename StateString::iterator StateStringIter;

  const Machine& machine;
  Writer& outs;
  StateString current;
  bool sentSOF, sentEOF;
  bool msb0;  // set this to encode MSB first, instead of LSB first

  Encoder (const Machine& machine, Writer& outs)
    : machine(machine),
      outs(outs),
      msb0(false),
      sentSOF(false),
      sentEOF(false)
  {
    current[machine.startState()] = deque<OutputSymbol>();
    expand();
  }

  ~Encoder() {
    close();
  }

  void close() {
    if (!sentEOF)
      encodeSymbol (MachineEOF);

    if (current.size()) {
      expand();
      vguard<StateStringIter> ssIter;
      for (StateStringIter ss = current.begin(); ss != current.end(); ++ss)
	if (machine.state[ss->first].isEnd())
	  ssIter.push_back (ss);
      if (ssIter.size() == 1)
	flush (ssIter.front());
      else if (ssIter.size() > 1) {
	Warn ("Encoder unresolved: %u possible end states", ssIter.size());
	for (auto ss: ssIter)
	  Warn ("State %s: output queue %s", machine.state[ss->first].name.c_str(), ss->second.empty() ? "empty" : to_string_join(ss->second,"").c_str());
      } else if (current.size() > 1) {
	Warn ("Encoder unresolved: %u possible states", current.size());
	showQueue();
      }
      current.clear();
    }
  }

  void showQueue() const {
    for (const auto& ss: current)
      Warn ("State %s: output queue %s", machine.state[ss.first].name.c_str(), ss.second.empty() ? "empty" : to_string_join(ss.second,"").c_str());
  }

  bool atEnd() const {
    return current.size() == 1 && machine.state[(*current.begin()).first].isEnd();
  }

  bool canEncodeSymbol (InputSymbol sym) const {
    for (const auto& ss: current)
      if (machine.state[ss.first].transFor(sym) != NULL)
	return true;
    return false;
  }
  
  void expand() {
    StateString next, seen;
    bool foundNew;
    do {
      foundNew = false;
      for (const auto& ss: current) {
	seen.insert (ss);
	const State state = ss.first;
	const auto& str = ss.second;
	const MachineState& ms = machine.state[state];
	LogThisAt(10,"Output queue for " << ms.name << " is " << (str.empty() ? string("empty") : string(str.begin(),str.end())) << endl);
	if (ms.isEnd() || ms.exitsWithInput())
	  next[state] = str;
      }
      for (const auto& ss: current) {
	const State state = ss.first;
	const auto& str = ss.second;
	const MachineState& ms = machine.state[state];
	for (const auto& t: ms.trans)
	  if (t.inputEmpty()) {
	    auto nextStr = str;
	    if (!t.outputEmpty())
	      nextStr.push_back (t.out);
	    if (seen.count (t.dest))
	      Assert (seen.at(t.dest) == nextStr,
		      "Encoder error: state %s has two possible output queues (%s, %s)",
		      machine.state[t.dest].name.c_str(),
		      to_string_join(seen.at(t.dest),"").c_str(),
		      to_string_join(nextStr,"").c_str());
	    else {
	      next[t.dest] = nextStr;
	      LogThisAt(9,"Transition " << ms.name
			<< " -> " << machine.state[t.dest].name
			<< (nextStr.empty() ? string() : (string(": output queue ") + to_string_join(nextStr,"")))
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
      LogThisAt(9,"Flushing output queue: " << str << endl);
      write (str.c_str(), str.size());
      ss->second.clear();
    }
  }

  void encodeSymbol (InputSymbol inSym) {
    if (!sentSOF && inSym != MachineSOF && canEncodeSymbol(MachineSOF))
      encodeSymbol (MachineSOF);
    if (inSym != MachineFlush && !canEncodeSymbol(inSym)) {
	Warn ("Sending FLUSH. Depending on the code, this may insert extra bits!");
	encodeSymbol (MachineFlush);
    }
    if (!canEncodeSymbol(inSym))
      showQueue();
    LogThisAt(8,"Encoding " << Machine::charToString(inSym) << endl);
    if (inSym == MachineSOF)
      sentSOF = true;
    else if (inSym == MachineEOF)
      sentEOF = true;

    StateString next;
    for (const auto& ss: current) {
      const State state = ss.first;
      const auto& str = ss.second;
      for (const auto& t: machine.state[state].trans)
	if (t.in == inSym) {
	  const State nextState = t.dest;
	  auto nextStr = str;
	  if (!t.outputEmpty())
	    nextStr.push_back (t.out);
	  Assert (!next.count(nextState) || next.at(nextState) == nextStr,
		  "Encoder error: state %s has two possible output queues (%s, %s)",
		  machine.state[nextState].name.c_str(),
		  to_string_join(next.at(nextState),"").c_str(),
		  to_string_join(nextStr,"").c_str());
	  next[nextState] = nextStr;
	  LogThisAt(9,"Transition " << machine.state[state].name
		    << " -> " << machine.state[nextState].name
		    << ": "
		    << (nextStr.empty() ? string() : (string("output queue ") + to_string_join(nextStr,"") + ", "))
		    << "input " << t.in
		    << endl);
	}
    }
    Assert (!next.empty(), "Can't encode symbol '%c'", inSym);
    current.swap (next);
    expand();
    if (current.size() == 1) {
      auto iter = current.begin();
      const MachineState& ms = machine.state[iter->first];
      if (ms.emitsOutput())
	flush (iter);
    } else
      shiftResolvedSymbols();
  }

  void shiftResolvedSymbols() {
    while (true) {
      bool foundQueue = false, queueNonempty, firstCharSame;
      OutputSymbol firstChar;
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
	LogThisAt(9,"All output queues have '" << firstChar << "' as first symbol; shifting" << endl);
	write (&firstChar, 1);
	for (auto& ss: current)
	  ss.second.erase (ss.second.begin());
      } else
	break;
    }
  }

  void encodeBit (bool bit) {
    encodeSymbol (bit ? MachineBit1 : MachineBit0);
  }

  template<typename CharType>
  void encodeByte (CharType byte) {
    LogThisAt(7,"Encoding '" << (char)byte << "' (\\x" << hex << (int)byte << ")" << endl);
    if (msb0)
      for (int n = 7; n >= 0; --n)
	encodeBit (byte & (1 << n));
    else
      for (int n = 0; n <= 7; ++n)
	encodeBit (byte & (1 << n));
  }

  void encodeStream (istream& in) {
    istreambuf_iterator<char> iter(in), iterEnd;
    while (iter != iterEnd) {
      encodeByte (*iter);
      ++iter;
    }
  }

  void encodeString (const string& s) {
    for (auto c: s)
      encodeByte (c);
  }

  void encodeSymbolString (const string& s) {
    for (auto c: s)
      encodeSymbol (c);
  }
};

struct FastaWriter {
  ostream& outs;
  size_t col, colsPerLine;
  
  FastaWriter (ostream& outs, const char* seqname = "SEQ")
    : outs(outs),
      col(0),
      colsPerLine(seqname == NULL ? 0 : 50)
  {
    if (seqname != NULL)
      outs << ">" << seqname << endl;
  }

  ~FastaWriter() {
    if (col > 0)
      outs << endl;
  }
  
  void write (char* buf, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      outs << buf[i];
      if (++col >= colsPerLine && colsPerLine > 0) {
	outs << endl;
	col = 0;
      }
    }
  }
};

#endif /* ENCODER_INCLUDED */
