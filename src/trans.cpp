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

bool MachineState::isInput() const {
  return transFor('0') != NULL && transFor('1') != NULL;
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

double Machine::expectedBasesPerBit() const {
  const size_t len = leftContextWidth() + rightContextWidth();
  const size_t burnInSteps = len*4, simSteps = len*4;
  map<State,double> current;
  size_t nSources = 0;
  for (State s = 0; s < nStates(); ++s)
    if (state[s].isInput()) {
      current[s] = 1;
      ++nSources;
    }
  Assert (nSources > 0, "Couldn't find any input states");
  for (auto& ps: current)
    ps.second /= (double) nSources;
  vguard<double> bpb;
  auto evolve = [&]() -> void {
    map<State,double> next;
    double basesPerBit = 0;
    for (const auto& ps: current) {
      const MachineState& ms = state[ps.first];
      const double p = ps.second;

      auto t0 = ms.transFor('0');
      State s0;
      while (true) {
	if (t0->out)
	  basesPerBit += p / 2;
	s0 = t0->dest;
	if (state[s0].isInput())
	  break;
	t0 = &state[s0].next();
      }
      next[s0] += p / 2;

      auto t1 = ms.transFor('1');
      State s1;
      while (true) {
	if (t1->out)
	  basesPerBit += p / 2;
	s1 = t1->dest;
	if (state[s1].isInput())
	  break;
	t1 = &state[s1].next();
      }
      next[s1] += p / 2;
    }
    bpb.push_back (basesPerBit);
    current.swap (next);
  };
  ProgressLog (plogSim, 1);
  plogSim.initProgress ("Estimating compression rate");
  for (size_t step = 0; step < burnInSteps; ++step) {
    plogSim.logProgress (step / (double) (burnInSteps + simSteps), "burn-in step %u/%u", step, simSteps + burnInSteps);
    evolve();
  }
  bpb.clear();
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
  return accumulate (bpb.begin(), bpb.end(), 0.) / (double) bpb.size();
}
