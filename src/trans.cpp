#include <iomanip>
#include <fstream>
#include "trans.h"
#include "logger.h"
#include "jsonutil.h"

struct MachineTokenLookup {
  map<InputSymbol,InputToken> sym2tok;
  map<InputToken,InputSymbol> tok2sym;
  map<InputSymbol,string> sym2desc;
  void add (InputSymbol c, const char* s, const char* desc);
  string tokenDescriptionTable (const string& alphabet) const;
  MachineTokenLookup();
};
MachineTokenLookup machineTokenLookup;

void MachineTokenLookup::add (InputSymbol c, const char* s, const char* desc) {
  const string str (s);
  tok2sym[str] = c;
  sym2tok[c] = str;
  sym2desc[c] = string(desc);
}

MachineTokenLookup::MachineTokenLookup() {
  add (MachineNull, "NULL", "Null token");
  tok2sym[string()] = MachineNull;

  add (MachineBit0, "0", "Zero input bit (works in any context)");
  add (MachineBit1, "1", "One input bit (works in any context)");

  add (MachineFlush, "FLUSH", "Flush any queued input bits");

  add (MachineStrictBit0, "0%2", "Strict input bit 0 (works in radix-2 context)");
  add (MachineStrictBit1, "1%2", "Strict input bit 1 (works in radix-2 context)");

  add (MachineStrictTrit0, "0%3", "Strict input trit 0 (works in radix-3 context)");
  add (MachineStrictTrit1, "1%3", "Strict input trit 1 (works in radix-3 context)");
  add (MachineStrictTrit2, "2%3", "Strict input trit 2 (works in radix-3 context)");

  add (MachineStrictQuat0, "0%4", "Strict input quat 0 (works in radix-4 context)");
  add (MachineStrictQuat1, "1%4", "Strict input quat 1 (works in radix-4 context)");
  add (MachineStrictQuat2, "2%4", "Strict input quat 2 (works in radix-4 context)");
  add (MachineStrictQuat3, "3%4", "Strict input quat 3 (works in radix-4 context)");

  add (MachineSOF, "START", "Start-of-file control symbol");
  add (MachineEOF, "EOF", "End-of-file control symbol");

  for (InputSymbol c = MachineControlFirst; c <= MachineControlLast; ++c)
    add (c, (string("!") + (char) (c + 'a' - MachineControlFirst)).c_str(), "Control symbol");
}

string MachineTokenLookup::tokenDescriptionTable (const string& alphabet) const {
  size_t tw = 0;
  for (char c: alphabet)
    if (sym2tok.count(c))
      tw = max (tw, sym2tok.at(c).size());
  ostringstream out;
  for (char c: alphabet)
    if (sym2tok.count(c))
      out << c << ' ' << setw(tw) << sym2tok.at(c) << ' ' << sym2desc.at(c) << endl;
    else
      out << c << ' ' << setw(tw) << "?" << ' ' << "Unknown token" << endl;
  return out.str();
}

MachineTransition::MachineTransition()
{ }

MachineTransition::MachineTransition (InputSymbol in, char out, State dest)
  : in (in),
    out (out),
    dest (dest)
{ }

bool MachineTransition::inputEmpty() const {
  return in == MachineNull;
}

bool MachineTransition::outputEmpty() const {
  return out == MachineNull;
}

bool MachineTransition::isNull() const {
  return in == MachineNull && out == MachineNull;
}

bool MachineTransition::isEOF() const {
  return in == MachineEOF;
}

bool MachineTransition::isSOF() const {
  return in == MachineSOF;
}

MachineState::MachineState()
{ }

const MachineTransition* MachineState::transFor (InputSymbol in) const {
  for (const auto& t: trans)
    if (t.in == in)
      return &t;
  return NULL;
}

bool MachineState::isEnd() const {
  return trans.empty();
}

bool MachineState::exitsWithInput (const char* symbols) const {
  for (const auto& t: trans)
    if (t.in && strchr (symbols, t.in) != NULL)
      return true;
  return false;
}

bool MachineState::exitsWithInput() const {
  for (const auto& t: trans)
    if (t.in)
      return true;
  return false;
}

bool MachineState::exitsWithoutInput() const {
  for (const auto& t: trans)
    if (!t.in)
      return true;
  return false;
}

bool MachineState::isWait() const {
  return exitsWithInput() && !exitsWithoutInput();
}

bool MachineState::isNonWait() const {
  return !exitsWithInput() && exitsWithoutInput();
}

bool MachineState::emitsOutput() const {
  for (const auto& t: trans)
    if (t.out)
      return true;
  return false;
}

bool MachineState::isDeterministic() const {
  return trans.size() == 1 && trans.front().in == 0;
}

const MachineTransition& MachineState::next() const {
  Assert (isDeterministic(), "Called next() method on a non-deterministic state");
  return trans.front();
}

State Machine::nStates() const {
  return state.size();
}

State Machine::startState() const {
  Assert (nStates() > 0, "Machine has no states");
  return 0;
}

Machine::Machine()
{ }

void Machine::writeDot (ostream& out) const {
  out << "digraph G {\n";
  for (State s = 0; s < nStates(); ++s) {
    const MachineState& ms = state[s];
    out << " " << s << " [label=\"";
    out << ms.name << " " << ms.leftContext << " " << ms.rightContext;
    out << "\"];" << endl;
  }
  out << endl;
  for (State s = 0; s < nStates(); ++s) {
    const MachineState& ms = state[s];
    for (const auto& t: ms.trans) {
      out << " " << s << " -> " << t.dest << " [label=\"" << charToString(t.in) << "/";
      if (!t.outputEmpty())
	out << t.out;
      out << "\"];" << endl;
    }
    out << endl;
  }
  out << "}" << endl;
}

void Machine::write (ostream& out) const {
  const size_t iw = stateIndexWidth();
  const size_t nw = stateNameWidth();
  const size_t lw = maxLeftContext();
  const size_t rw = maxRightContext();
  for (State s = 0; s < nStates(); ++s) {
    const MachineState& ms = state[s];
    out << setw(iw+1) << left << stateIndex(s)
	<< setw(nw+1) << left << ms.name
	<< setw(lw) << right << ms.leftContext
	<< "."
	<< setw(rw) << left << ms.rightContext;
    for (const auto& t: ms.trans) {
      out << " ";
      if (t.in) out << charToString (t.in);
      out << "/";
      if (!t.outputEmpty())
	out << t.out;
      out << "->" << stateIndex(t.dest);
    }
    out << endl;
  }
}

string Machine::stateIndex (State s) {
  return string("#") + to_string(s);
}

InputSymbol Machine::controlChar (ControlIndex c) {
  InputSymbol ci = c + MachineControlFirst;
  Assert (ci <= MachineControlLast, "Ran out of control chars");
  return ci;
}

ControlIndex Machine::controlIndex (InputSymbol c) {
  return isControl(c) ? (c - MachineControlFirst) : -1;
}

bool Machine::isControl (InputSymbol c) {
  return c >= MachineControlFirst && c <= MachineControlLast;
}

bool Machine::isStrict (InputSymbol c) {
  return c == MachineStrictBit0 || c == MachineStrictBit1
    || c == MachineStrictTrit0 || c == MachineStrictTrit1 || c == MachineStrictTrit2
    || c == MachineStrictQuat0 || c == MachineStrictQuat1 || c == MachineStrictQuat2 || c == MachineStrictQuat3;
}

bool Machine::isRelaxed (InputSymbol c) {
  return c == MachineBit0 || c == MachineBit1;
}

InputToken Machine::charToString (InputSymbol c) {
  if (machineTokenLookup.sym2tok.count(c))
    return InputToken (machineTokenLookup.sym2tok.at(c));
  return InputToken("???");
}

InputSymbol Machine::stringToChar (const InputToken& in) {
  if (machineTokenLookup.tok2sym.count(in))
    return machineTokenLookup.tok2sym.at(in);
  return -1;
}

size_t Machine::maxLeftContext() const {
  size_t w = 0;
  for (const auto& ms: state)
    w = max (w, ms.leftContext.size());
  return w;
}

size_t Machine::maxRightContext() const {
  size_t w = 0;
  for (const auto& ms: state)
    w = max (w, ms.rightContext.size());
  return w;
}

size_t Machine::stateNameWidth() const {
  size_t w = 0;
  for (const auto& ms: state)
    w = max (w, ms.name.size());
  return w;
}

size_t Machine::stateIndexWidth() const {
  size_t w = 0;
  for (State s = 0; s < nStates(); ++s)
    w = max (w, stateIndex(s).size());
  return w;
}

string Machine::inputAlphabet (int inputFlags) const {
  set<char> alph;
  for (const auto& ms: state)
    for (const auto& t: ms.trans)
      if (!t.inputEmpty()
	  && (((t.isEOF() || t.isSOF()) && (inputFlags & MachineSEOFInputFlag))
	      || (isControl(t.in) && (inputFlags & MachineControlInputFlag))
	      || (t.in == MachineFlush && (inputFlags & MachineFlushInputFlag))
	      || (isRelaxed(t.in) && (inputFlags & MachineRelaxedInputFlag))
	      || (isStrict(t.in) && (inputFlags & MachineStrictInputFlag))))
	alph.insert (t.in);
  return string (alph.begin(), alph.end());
}

string Machine::outputAlphabet() const {
  set<char> alph;
  for (const auto& ms: state)
    for (const auto& t: ms.trans)
      if (!t.outputEmpty())
	alph.insert (t.out);
  return string (alph.begin(), alph.end());
}

string Machine::inputDescriptionTable() const {
  return machineTokenLookup.tokenDescriptionTable (inputAlphabet (MachineAllInputFlags));
}

map<InputSymbol,double> Machine::expectedBasesPerInputSymbol (const char* symbols) const {
  const size_t len = maxLeftContext() + maxRightContext();
  const size_t burnInSteps = len*4, simSteps = len*4;
  map<State,double> current;
  size_t nSources = 0;
  for (State s = 0; s < nStates(); ++s)
    if (state[s].exitsWithInput (symbols)) {
      current[s] = 1;
      ++nSources;
    }
  Assert (nSources > 0, "Couldn't find any input states");
  for (auto& ps: current)
    ps.second /= (double) nSources;
  const string alph (symbols);
  map<char,vguard<double> > eb;

  auto logCurrent = [&](int logLevel) -> void {
    double pTot = 0;
    for (auto& ps: current) {
      LogThisAt(logLevel,"P(" << state[ps.first].name << ") = " << ps.second << endl);
      pTot += ps.second;
    }
    LogThisAt(logLevel+1,"Total probability is " << pTot << endl);
  };
  
  auto evolve = [&]() -> void {
    logCurrent(5);
    map<State,double> next;
    map<char,double> bases;
    for (const auto& ps: current) {
      const MachineState& ms = state[ps.first];
      const double p = ps.second;

      map<char,const MachineTransition*> transFor;
      double nt = 0;
      for (char c: alph) {
	auto t = ms.transFor(c);
	if (t) {
	  transFor[c] = t;
	  if (c != MachineEOF)
	    ++nt;
	}
      }

      double ptot = 0;
      for (const auto& ct: transFor) {
	const char c = ct.first;
	auto t = ct.second;
	set<State> seen;
	State s;
	while (true) {
	  if (t->out)
	    bases[c] += p;
	  s = t->dest;
	  if (seen.count(s))  // guard against infinite loops
	    break;
	  seen.insert(s);
	  if (state[s].exitsWithInput() || state[s].isEnd())
	    break;
	  Assert (state[s].isDeterministic(), "Non-deterministic state without inputs: %s", state[s].name.c_str());
	  t = &state[s].next();
	}
	if (c != MachineEOF) {
	  LogThisAt(7,"P(" << ms.name << "->" << state[s].name << ")=" << (p/nt) << endl);
	  next[s] += p / nt;
	  ptot += p / nt;
	}
      }
      LogThisAt(8,"Total outgoing transition probability from state " << ms.name << " is " << ptot << "; state probability is " << p << endl);
    }
    for (char c: alph)
      eb[c].push_back (bases[c]);
    current.swap (next);
  };

  ProgressLog (plogSim, 1);
  plogSim.initProgress ("Estimating compression rate");
  for (size_t step = 0; step < burnInSteps; ++step) {
    plogSim.logProgress (step / (double) (burnInSteps + simSteps), "burn-in step %u/%u", step, simSteps + burnInSteps);
    evolve();
  }
  eb.clear();
  for (size_t step = 0; step < simSteps; ++step) {
    plogSim.logProgress ((step + burnInSteps) / (double) (burnInSteps + simSteps), "step %u/%u", step + burnInSteps, simSteps + burnInSteps);
    evolve();
  }

  logCurrent(3);

  map<char,double> bps;
  for (char c: alph)
    bps[c] = accumulate (eb[c].begin(), eb[c].end(), 0.) / (double) eb[c].size();
  return bps;
}

void Machine::writeJSON (ostream& out) const {
  out << "{\"state\": [" << endl;
  for (State s = 0; s < nStates(); ++s) {
    const MachineState& ms = state[s];
    out << " {\"n\":" << s << ",";
    if (ms.name.size())
	out << "\"id\":\"" << ms.name << "\",";
    if (ms.leftContext.size())
      out << "\"l\":\"" << ms.leftContext << "\",";
    if (ms.rightContext.size())
      out << "\"r\":\"" << ms.rightContext << "\",";
    out << "\"trans\":[";
    for (size_t nt = 0; nt < ms.trans.size(); ++nt) {
      const MachineTransition& t = ms.trans[nt];
      if (nt > 0) out << ",";
      out << "{";
      if (t.in) out << "\"in\":\"" << t.in << "\",";
      if (t.out) out << "\"out\":\"" << t.out << "\",";
      out << "\"to\":" << t.dest
	  << "}";
    }
    out << "]}";
    if (s < nStates() - 1)
      out << ",";
    out << endl;
  }
  out << "]}" << endl;
}

void Machine::readJSON (istream& in) {
  state.clear();
  ParsedJson pj (in);
  JsonValue jstate = pj.getType ("state", JSON_ARRAY);
  for (JsonIterator iter = begin(jstate); iter != end(jstate); ++iter) {
    const JsonMap jsmap (iter->value);
    MachineState ms;
    if (jsmap.contains("n")) {
      const size_t n = jsmap.getNumber("n");
      Require (state.size() == n, "State n=%u out of sequence", n);
    }
    if (jsmap.contains("id"))
      ms.name = jsmap.getString("id");
    if (jsmap.contains("l"))
      ms.leftContext = jsmap.getString("l");
    if (jsmap.contains("r"))
      ms.rightContext = jsmap.getString("r");
    JsonValue jtrans = jsmap.getType ("trans", JSON_ARRAY);
    for (JsonIterator transIter = begin(jtrans); transIter != end(jtrans); ++transIter) {
      const JsonMap jtmap (transIter->value);
      MachineTransition t;
      t.in = t.out = 0;
      t.dest = (State) jtmap.getNumber("to");
      if (jtmap.contains("in")) {
	const string& tin = jtmap.getString("in");
	Assert (tin.size() == 1, "Invalid input character: %s", tin.c_str());
	t.in = tin[0];
      }
      if (jtmap.contains("out")) {
	const string& tout = jtmap.getString("out");
	Assert (tout.size() == 1, "Invalid output character: %s", tout.c_str());
	t.out = tout[0];
      }
      ms.trans.push_back (t);
    }
    state.push_back (ms);
  }
  verifyContexts();
  Assert (isWaitingMachine(), "Not a waiting machine");
}

Machine Machine::fromJSON (istream& in) {
  Machine machine;
  machine.readJSON (in);
  return machine;
}

Machine Machine::fromFile (const char* filename) {
  ifstream infile (filename);
  if (!infile)
    Fail ("File not found: %s", filename);
  return fromJSON (infile);
}

void Machine::verifyContexts() const {
  for (const auto& ms: state) {
    for (const auto& t: ms.trans) {
      const auto& md = state[t.dest];
      if (t.out) {
	if (ms.rightContext.size())
	  Assert (t.out == ms.rightContext[0], "In transition from %s to %s: emitted character (%c) does not match source's right context (%s)", ms.name.c_str(), md.name.c_str(), t.out, ms.rightContext.c_str());
	if (md.leftContext.size())
	  Assert (t.out == md.leftContext.back(), "In transition from %s to %s: emitted character (%c) does not match destination's left context (%s)", ms.name.c_str(), md.name.c_str(), t.out, md.leftContext.c_str());
      }
    }
  }
}

bool Machine::isWaitingMachine() const {
  for (const auto& ms: state)
    if (!ms.isWait() && !ms.isNonWait() && !ms.isEnd())
      return false;
  return true;
}

Machine Machine::compose (const Machine& first, const Machine& second) {
  LogThisAt(3,"Composing " << first.nStates() << "-state transducer with " << second.nStates() << "-state transducer" << endl);
  Assert (second.isWaitingMachine(), "Attempt to compose transducers A*B where B is not a waiting machine");
  vguard<MachineState> comp (first.nStates() * second.nStates());
  auto compState = [&](State i,State j) -> State {
    return i * second.nStates() + j;
  };
  auto compStateName = [&](State i,State j) -> string {
    return string("(") + first.state[i].name + "," + second.state[j].name + ")";
  };
  for (State i = 0; i < first.nStates(); ++i)
    for (State j = 0; j < second.nStates(); ++j) {
      MachineState& ms = comp[compState(i,j)];
      const MachineState& msi = first.state[i];
      const MachineState& msj = second.state[j];
      ms.name = compStateName(i,j);
      ms.leftContext = msj.leftContext;
      ms.rightContext = msj.rightContext;
      if (msj.isWait()) {
	for (const auto& it: msi.trans)
	  if (it.out == MachineNull) {
	    ms.trans.push_back (MachineTransition (it.in, MachineNull, compState(it.dest,j)));
	    LogThisAt(6,"Adding transition from " << ms.name << " to " << compStateName(it.dest,j) << endl);
	  } else
	    for (const auto& jt: msj.trans)
	      if (it.out == jt.in) {
		ms.trans.push_back (MachineTransition (it.in, jt.out, compState(it.dest,jt.dest)));
		LogThisAt(6,"Adding transition from " << ms.name << " to " << compStateName(it.dest,jt.dest) << endl);
	      }
      } else
	for (const auto& jt: msj.trans) {
	  ms.trans.push_back (MachineTransition (MachineNull, jt.out, compState(i,jt.dest)));
	  LogThisAt(6,"Adding transition from " << ms.name << " to " << compStateName(i,jt.dest) << endl);
	}
    }
  vguard<bool> seen (comp.size(), false);
  deque<State> queue;
  queue.push_back (compState(first.startState(),second.startState()));
  while (queue.size()) {
    const State c = queue.front();
    queue.pop_front();
    seen[c] = true;
    for (const auto& t: comp[c].trans)
      if (!seen[t.dest])
	queue.push_back (t.dest);
  }
  map<State,State> nullEquiv;
  for (State s = 0; s < comp.size(); ++s)
    if (seen[s]) {
      State d = s;
      while (comp[d].trans.size() == 1 && comp[d].trans.front().isNull())
	d = comp[d].trans.front().dest;
      if (d != s)
	nullEquiv[s] = d;
    }
  vguard<State> old2new (comp.size());
  State nStates = 0;
  for (State oldIdx = 0; oldIdx < comp.size(); ++oldIdx)
    if (seen[oldIdx] && !nullEquiv.count(oldIdx))
      old2new[oldIdx] = nStates++;
  for (State oldIdx = 0; oldIdx < comp.size(); ++oldIdx)
    if (seen[oldIdx] && nullEquiv.count(oldIdx))
      old2new[oldIdx] = old2new[nullEquiv.at(oldIdx)];
  for (State oldIdx = 0; oldIdx < comp.size(); ++oldIdx)
    if (seen[oldIdx])
      for (auto& t: comp[oldIdx].trans)
	t.dest = old2new[t.dest];
  LogThisAt(3,"Transducer composition yielded " << nStates << "-state machine; " << plural (comp.size() - nStates, "more state was", "more states were") << " unreachable" << endl);
  Machine compMachine;
  compMachine.state.reserve (nStates);
  for (State oldIdx = 0; oldIdx < comp.size(); ++oldIdx)
    if (seen[oldIdx] && !nullEquiv.count(oldIdx))
      compMachine.state.push_back (comp[oldIdx]);
  return compMachine;
}

vguard<State> Machine::decoderToposort (const string& inputAlphabet) const {
  LogThisAt(5,"Toposorting transducer for decoder" << endl);
  deque<State> S;
  vguard<State> L;
  vguard<int> nParents (nStates());
  vguard<vguard<State> > children (nStates());
  int edges = 0;
  for (State s = 0; s < nStates(); ++s)
    for (const auto& t: state[s].trans)
      if (t.outputEmpty() && (t.inputEmpty() || inputAlphabet.find(t.in) != string::npos)) {
	++nParents[t.dest];
	++edges;
	children[s].push_back (t.dest);
      }
  for (State s = 0; s < nStates(); ++s)
    if (nParents[s] == 0)
      S.push_back (s);
  while (S.size()) {
    const State n = S.front();
    S.pop_front();
    L.push_back (n);
    for (auto m : children[n]) {
      --edges;
      if (--nParents[m] == 0)
	S.push_back (m);
    }
  }
  if (edges > 0)
    throw std::domain_error ("Transducer is cyclic, can't toposort");
  return L;
}
