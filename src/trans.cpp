#include <iomanip>
#include "trans.h"
#include "logger.h"

MachineTransition::MachineTransition (char in, char out, State dest)
  : in (in),
    out (out),
    dest (dest)
{ }

MachineState::MachineState()
  : type (UndefinedState)
{ }

string MachineState::typeString() const {
  switch (type) {
  case SourceState:
    return string("Source");
  case ControlState:
    return string("Meta(") + Machine::controlChar(control) + ")";
  case CodeState:
    return string("Code");
  case SplitState:
    return string("Split");
  case PadState:
    return string("Pad(") + Machine::controlChar(control) + ")";
  default:
    break;
  }
  return string("Undefined");
}

const MachineTransition* MachineState::transFor (char in) const {
  for (const auto& t: trans)
    if (t.in == in)
      return &t;
  return NULL;
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

Machine::Machine (Pos len)
  : len(len)
{ }

void Machine::write (ostream& out) const {
  const size_t tw = typeStringWidth();
  const size_t sw = stateNameWidth();
  for (State s = 0; s < nStates(); ++s) {
    out << setw(sw+1) << left << stateName(s)
	<< setw(tw+1) << left << state[s].typeString()
	<< kmerString(state[s].context,len);
    for (const auto& t: state[s].trans) {
      out << " ";
      if (t.in) out << t.in;
      out << "/";
      if (t.out) out << t.out;
      out << "->" << stateName(t.dest);
    }
    out << endl;
  }
}

string Machine::stateName (State s) {
  return string("#") + to_string(s+1);
}

string controlChars ("XYZWVPQRSDEFIJKLM23456789<>[]!?abcdefghijklmnopqrstuvwxyz");
char Machine::controlChar (ControlIndex c) {
  Assert (c < controlChars.size(), "Ran out of control cahrs");
  return controlChars[c];
}

size_t Machine::stateNameWidth() const {
  size_t w = 0;
  for (State s = 0; s < nStates(); ++s)
    w = max (w, stateName(s).size());
  return w;
}

size_t Machine::typeStringWidth() const {
  size_t w = 0;
  for (const auto& ms: state)
    w = max (w, ms.typeString().size());
  return w;
}

double Machine::expectedBasesPerBit() const {
  const size_t burnInSteps = len*4, simSteps = len*4;
  map<State,double> current;
  size_t nSources = 0;
  for (State s = 0; s < nStates(); ++s)
    if (state[s].type == SourceState || state[s].type == ControlState) {
      State sIn = s;
      while (!state[sIn].isInput())
	sIn = state[sIn].next().dest;
      current[sIn] = 1;
      ++nSources;
    }
  if (nSources == 0)
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
    LogThisAt(3,"P(" << kmerString(state[ps.first].context,len) << ") = " << ps.second << endl);
    pTot += ps.second;
  }
  LogThisAt(4,"Total probability is " << pTot << endl);
  return accumulate (bpb.begin(), bpb.end(), 0.) / (double) bpb.size();
}
