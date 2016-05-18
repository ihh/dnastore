#include <iomanip>
#include "trans.h"
#include "logger.h"

MachineTransition::MachineTransition (char in, char out, State dest)
  : in (in),
    out (out),
    dest (dest)
{ }

MachineState::MachineState()
{ }

const MachineTransition* MachineState::transFor (char in) const {
  for (const auto& t: trans)
    if (t.in == in)
      return &t;
  return NULL;
}

bool MachineState::isEnd() const {
  return trans.empty();
}

bool MachineState::hasInput() const {
  for (const auto& t: trans)
    if (t.in)
      return true;
  return false;
}

bool MachineState::hasOutput() const {
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
	<< setw(rw+1) << left << ms.rightContext;
    for (const auto& t: ms.trans) {
      out << " ";
      if (t.in) out << t.in;
      out << "/";
      if (t.out) out << t.out;
      out << "->" << stateIndex(t.dest);
    }
    out << endl;
  }
}

string Machine::stateIndex (State s) {
  return string("#") + to_string(s+1);
}

string controlChars ("23456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
char Machine::controlChar (ControlIndex c) {
  Assert (c < controlChars.size(), "Ran out of control chars");
  return controlChars[c];
}

ControlIndex Machine::controlIndex (char c) {
  const char* cc = controlChars.c_str();
  const char* s = strchr (cc, c);
  return s == NULL ? -1 : s - cc;
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
      if (t.in)
	alph.insert (t.in);
  return string (alph.begin(), alph.end());
}

string Machine::outputAlphabet() const {
  set<char> alph;
  for (const auto& ms: state)
    for (const auto& t: ms.trans)
      if (t.out)
	alph.insert (t.out);
  return string (alph.begin(), alph.end());
}

map<char,double> Machine::expectedBasesPerInputSymbol() const {
  const size_t len = leftContextWidth() + rightContextWidth();
  const size_t burnInSteps = len*4, simSteps = len*4;
  map<State,double> current;
  size_t nSources = 0;
  for (State s = 0; s < nStates(); ++s)
    if (state[s].hasInput()) {
      current[s] = 1;
      ++nSources;
    }
  Assert (nSources > 0, "Couldn't find any input states");
  for (auto& ps: current)
    ps.second /= (double) nSources;
  const string alph = inputAlphabet();
  map<char,vguard<double> > eb;
  auto evolve = [&]() -> void {
    map<State,double> next;
    map<char,double> bases;
    for (const auto& ps: current) {
      const MachineState& ms = state[ps.first];
      const double p = ps.second;

      map<char,const MachineTransition*> transFor;
      for (char c: alph) {
	auto t = ms.transFor(c);
	if (t)
	  transFor[c] = t;
      }

      const double nt = transFor.count('0') + transFor.count('1');

      for (const auto& ct: transFor) {
	const char c = ct.first;
	auto t = ct.second;
	if (t->out)
	  bases[c] += p;
	State s;
	while (true) {
	  if (t->out)
	    bases[c] += p;
	  s = t->dest;
	  if (state[s].hasInput())
	    break;
	  Assert (state[s].isDeterministic(), "Non-deterministic state without inputs");
	  t = &state[s].next();
	}
	if (c == '0' || c == '1')
	  next[s] += p / nt;
      }
    }
    for (const auto& cb: bases)
      eb[cb.first].push_back (cb.second);
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
  for (const auto& cb: eb)
    bps[cb.first] = accumulate (cb.second.begin(), cb.second.end(), 0.) / (double) cb.second.size();
  return bps;
}
