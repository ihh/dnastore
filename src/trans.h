#ifndef TRANSDUCER_INCLUDED
#define TRANSDUCER_INCLUDED

#include <string>
#include <map>
#include "vguard.h"

using namespace std;

typedef unsigned long long State;

typedef int ControlIndex;

struct MachineTransition {
  char in, out;
  State dest;
  MachineTransition();
  MachineTransition (char, char, State);
  bool isInput() const;
  bool isOutput() const;
  bool isEof() const;
};

struct MachineState {
  string name, leftContext, rightContext;
  vguard<MachineTransition> trans;
  MachineState();
  const MachineTransition* transFor (char in) const;
  bool isEnd() const;  // true if this has no outgoing transitions
  bool acceptsInputOrEof() const;
  bool emitsOutput() const;
  bool isDeterministic() const;  // true if this has only one non-absorbing transition
  const MachineTransition& next() const;
};

struct Machine {
  vguard<MachineState> state;
  
  Machine();
  State nStates() const;
  State startState() const;
  
  void verifyContexts() const;

  void write (ostream& out) const;
  void writeJSON (ostream& out) const;
  void readJSON (istream& in);
  static Machine fromJSON (istream& in);
  static Machine fromFile (const char* filename);
  
  static char controlChar (ControlIndex c);
  static ControlIndex controlIndex (char c);
  static string stateIndex (State s);

  static char eofChar;
  static string charToString (char in);
  static char stringToChar (const string& in);
  
  size_t leftContextWidth() const;
  size_t rightContextWidth() const;
  size_t stateNameWidth() const;
  size_t stateIndexWidth() const;

  string inputAlphabet() const;
  string outputAlphabet() const;

  map<char,double> expectedBasesPerInputSymbol (bool includeEof = false) const;
};

#endif /* TRANSDUCER_INCLUDED */
