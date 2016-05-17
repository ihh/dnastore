#include <iomanip>
#include "trans.h"

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

string controlChars ("XYPQVWKLEFIJLM23456789<>[]!?abcdefghijklmnopqrstuvwxyz");
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
