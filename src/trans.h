#ifndef TRANSDUCER_INCLUDED
#define TRANSDUCER_INCLUDED

#include "kmer.h"

typedef unsigned long long State;
typedef char StateType;

struct MachineTransition {
  char in, out;
  State dest;
  MachineTransition (char, char, State);
};

struct MachineState {
  Kmer context;
  StateType type;  // 0 for regular coding states
  vguard<MachineTransition> trans;
  MachineState();
  MachineState (Kmer, StateType);
};

struct Machine {
  const Pos len;
  map<char,Kmer> control;
  vguard<MachineState> state;

  Machine (Pos len);
  State nStates() const;

  void write (ostream& out) const;
  static string stateName (State s);
  static string typeName (StateType t);
};

#endif /* TRANSDUCER_INCLUDED */
