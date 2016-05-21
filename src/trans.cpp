#include <iomanip>
#include <fstream>
#include "trans.h"
#include "logger.h"
#include "jsonutil.h"

struct MachineTokenLookup {
  map<InputSymbol,InputToken> sym2tok;
  map<InputToken,InputSymbol> tok2sym;
  void add (InputSymbol c, const char* s);
  MachineTokenLookup();
};
MachineTokenLookup machineTokenLookup;

void MachineTokenLookup::add (InputSymbol c, const char* s) {
  const string str (s);
  tok2sym[str] = c;
  sym2tok[c] = str;
}

MachineTokenLookup::MachineTokenLookup() {
  add (MachineNull, "NULL");
  tok2sym[string()] = MachineNull;

  add (MachineBit0, "0");
  add (MachineBit1, "1");

  add (MachineFlushedBit0, "0.");
  add (MachineFlushedBit1, "1.");

  add (MachineEscapedA, "A");
  add (MachineEscapedC, "C");
  add (MachineEscapedG, "G");
  add (MachineEscapedT, "T");

  add (MachineStrictBit0, "0%2");
  add (MachineStrictBit1, "1%2");

  add (MachineStrictTrit0, "0%3");
  add (MachineStrictTrit1, "1%3");
  add (MachineStrictTrit2, "2%3");

  add (MachineStrictQuat0, "0%4");
  add (MachineStrictQuat1, "1%4");
  add (MachineStrictQuat2, "2%4");
  add (MachineStrictQuat3, "3%4");

  add (MachineEOF, "EOF");

  for (InputSymbol c = MachineControlFirst; c <= MachineControlLast; ++c)
    add (c, (string("^") + (char) (c + 'a' - MachineControlFirst)).c_str());
}

MachineTransition::MachineTransition()
{ }

MachineTransition::MachineTransition (InputSymbol in, char out, State dest)
  : in (in),
    out (out),
    dest (dest)
{ }

bool MachineTransition::inputPrintable() const {
  return in != MachineNull && in != MachineEOF;
}

bool MachineTransition::outputNonempty() const {
  return out != MachineNull;
}

bool MachineTransition::isEof() const {
  return in == MachineEOF;
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

bool MachineState::acceptsInputOrEof() const {
  for (const auto& t: trans)
    if (t.in)
      return true;
  return false;
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
      if (t.outputNonempty()) out << t.out;
      out << "\"];" << endl;
    }
    out << endl;
  }
  out << "}" << endl;
}

void Machine::write (ostream& out) const {
  const size_t iw = stateIndexWidth();
  const size_t nw = stateNameWidth();
  const size_t lw = leftContextWidth();
  const size_t rw = rightContextWidth();
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
      if (t.outputNonempty()) out << t.out;
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
  return (c >= MachineControlFirst && c <= MachineControlLast) ? (c - MachineControlFirst) : -1;
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

size_t Machine::leftContextWidth() const {
  size_t w = 0;
  for (const auto& ms: state)
    w = max (w, ms.leftContext.size());
  return w;
}

size_t Machine::rightContextWidth() const {
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

string Machine::inputAlphabet() const {
  set<char> alph;
  for (const auto& ms: state)
    for (const auto& t: ms.trans)
      if (t.inputPrintable())
	alph.insert (t.in);
  return string (alph.begin(), alph.end());
}

string Machine::outputAlphabet() const {
  set<char> alph;
  for (const auto& ms: state)
    for (const auto& t: ms.trans)
      if (t.outputNonempty())
	alph.insert (t.out);
  return string (alph.begin(), alph.end());
}

map<InputSymbol,double> Machine::expectedBasesPerInputSymbol (const char* symbols) const {
  const size_t len = leftContextWidth() + rightContextWidth();
  const size_t burnInSteps = len*4, simSteps = len*4;
  map<State,double> current;
  size_t nSources = 0;
  for (State s = 0; s < nStates(); ++s)
    if (state[s].acceptsInputOrEof()) {
      current[s] = 1;
      ++nSources;
    }
  Assert (nSources > 0, "Couldn't find any input states");
  for (auto& ps: current)
    ps.second /= (double) nSources;
  const string alph (symbols);
  map<char,vguard<double> > eb;
  auto evolve = [&]() -> void {
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
	  if (state[s].acceptsInputOrEof() || state[s].isEnd())
	    break;
	  Assert (state[s].isDeterministic(), "Non-deterministic state without inputs: %s", state[s].name.c_str());
	  t = &state[s].next();
	}
	if (c != MachineEOF)
	  next[s] += p / nt;
      }
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
  double pTot = 0;
  for (auto& ps: current) {
    LogThisAt(3,"P(" << state[ps.first].name << ") = " << ps.second << endl);
    pTot += ps.second;
  }
  LogThisAt(4,"Total probability is " << pTot << endl);
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
      if (t.in) out << "\"in\":\"" << charToString(t.in) << "\",";
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
      if (jtmap.contains("in"))
	t.in = stringToChar (jtmap.getString("in"));
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
