#ifndef TRANSDUCER_INCLUDED
#define TRANSDUCER_INCLUDED

#include "kmer.h"

typedef unsigned long long State;

enum StateType { SourceState,
		 ControlState,
		 CodeState,
		 SplitState,
		 PadState,
		 UndefinedState };

typedef size_t ControlIndex;

struct MachineTransition {
  char in, out;
  State dest;
  MachineTransition (char, char, State);
};

struct MachineState {
  Kmer context;
  StateType type;
  ControlIndex control;
  vguard<MachineTransition> trans;
  MachineState();
  MachineState (Kmer);
  string typeString() const;
  const MachineTransition* transFor (char in) const;
  bool isInput() const;
  bool isDeterministic() const;
  const MachineTransition& next() const;
};

struct Machine {
  const Pos len;
  vguard<MachineState> state;
  
  Machine (Pos len);
  State nStates() const;
  State startState() const;
  
  void write (ostream& out) const;

  static string stateName (State s);
  static char controlChar (ControlIndex c);

  size_t stateNameWidth() const;
  size_t typeStringWidth() const;

  double expectedBasesPerBit() const;
};

#endif /* TRANSDUCER_INCLUDED */
