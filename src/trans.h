#ifndef TRANSDUCER_INCLUDED
#define TRANSDUCER_INCLUDED

#include <string>
#include "vguard.h"

using namespace std;

typedef unsigned long long State;

typedef int ControlIndex;

struct MachineTransition {
  char in, out;
  State dest;
  MachineTransition (char, char, State);
};

struct MachineState {
  string name, leftContext, rightContext;
  vguard<MachineTransition> trans;
  MachineState();
  const MachineTransition* transFor (char in) const;
  bool isEnd() const;  // true if this has no outgoing transitions
  bool isInput() const;  // true if this has '0' and '1' inputs
  bool isDeterministic() const;  // true if this has only one non-absorbing transition
  const MachineTransition& next() const;
};

struct Machine {
  vguard<MachineState> state;
  
  Machine();
  State nStates() const;

  State startState() const;
  
  void write (ostream& out) const;

  static char controlChar (ControlIndex c);
  static ControlIndex controlIndex (char c);

  static string stateIndex (State s);

  size_t leftContextWidth() const;
  size_t rightContextWidth() const;
  size_t stateNameWidth() const;
  size_t stateIndexWidth() const;

  double expectedBasesPerBit() const;
};

#endif /* TRANSDUCER_INCLUDED */
