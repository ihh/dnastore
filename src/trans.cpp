#include <iomanip>
#include "trans.h"

MachineTransition::MachineTransition (char in, char out, State dest)
  : in (in),
    out (out),
    dest (dest)
{ }

MachineState::MachineState()
{ }

MachineState::MachineState (Kmer context, StateType type)
  : context (context),
    type (type)
{ }

State Machine::nStates() const {
  return state.size();
}

Machine::Machine (Pos len)
  : len(len)
{ }

void Machine::write (ostream& out) const {
  for (State s = 0; s < nStates(); ++s) {
    out << setw(6) << stateName(s)
	<< " " << setw(3) << typeName(state[s].type)
	<< " " << kmerString(state[s].context,len);
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

string Machine::typeName (StateType t) {
  return t == 0 ? string(".") : string(1,t);
}
