#ifndef TRANSDUCER_INCLUDED
#define TRANSDUCER_INCLUDED

#include <string>
#include <map>
#include "vguard.h"

using namespace std;

typedef unsigned long long State;

typedef int ControlIndex;

#define MachineNull         '\0'

#define MachineBit0         '0'
#define MachineBit1         '1'

#define MachineFlush        '.'

#define MachineStrictBit0   'i'
#define MachineStrictBit1   'j'

#define MachineStrictTrit0  'x'
#define MachineStrictTrit1  'y'
#define MachineStrictTrit2  'z'

#define MachineStrictQuat0  'p'
#define MachineStrictQuat1  'q'
#define MachineStrictQuat2  'r'
#define MachineStrictQuat3  's'

#define MachineSOF          '^'
#define MachineEOF          '$'

#define MachineControlFirst 'A'
#define MachineControlLast  'Z'

#define MachineWildContext  '*'

#define MachineStrictInputFlag   1
#define MachineRelaxedInputFlag  2
#define MachineFlushInputFlag    4
#define MachineControlInputFlag  8
#define MachineSEOFInputFlag     16
#define MachineDefaultInputFlags (MachineRelaxedInputFlag | MachineControlInputFlag)
#define MachineAllInputFlags     (MachineStrictInputFlag | MachineRelaxedInputFlag | MachineFlushInputFlag | MachineControlInputFlag |  MachineSEOFInputFlag)

typedef char OutputSymbol;
typedef char InputSymbol;
typedef string InputToken;

struct MachineTransition {
  InputSymbol in;
  OutputSymbol out;
  State dest;
  MachineTransition();
  MachineTransition (InputSymbol, OutputSymbol, State);
  bool inputEmpty() const;
  bool outputEmpty() const;
  bool isNull() const;
  bool isSOF() const;
  bool isEOF() const;
};

struct MachineState {
  string name, leftContext, rightContext;
  vguard<MachineTransition> trans;
  MachineState();
  const MachineTransition* transFor (InputSymbol in) const;
  bool isEnd() const;  // true if this has no outgoing transitions
  bool exitsWithInput (const char* symbols) const;  // true if this has an input transition for the specified symbols
  bool exitsWithInput() const;  // true if this has an input transition
  bool exitsWithoutInput() const;  // true if this has a non-input transition
  bool emitsOutput() const;  // true if this has an output transition
  bool isDeterministic() const;  // true if this has only one transition and it is non-input
  bool isWait() const;  // exitsWithInput() && !exitsWithoutInput()
  bool isNonWait() const;  // !exitsWithInput() && exitsWithoutInput()
  const MachineTransition& next() const;
};

struct Machine {
  vguard<MachineState> state;
  
  Machine();
  State nStates() const;
  State startState() const;
  
  void verifyContexts() const;
  bool isWaitingMachine() const;

  static Machine compose (const Machine& first, const Machine& second);
  
  void write (ostream& out) const;
  void writeDot (ostream& out) const;
  void writeJSON (ostream& out) const;
  void readJSON (istream& in);
  static Machine fromJSON (istream& in);
  static Machine fromFile (const char* filename);
  
  static InputSymbol controlChar (ControlIndex c);
  static ControlIndex controlIndex (InputSymbol c);

  static bool isControl (InputSymbol c);
  static bool isStrict (InputSymbol c);
  static bool isRelaxed (InputSymbol c);
  
  static InputToken charToString (InputSymbol in);
  static InputSymbol stringToChar (const InputToken& in);
  
  static string stateIndex (State s);

  size_t maxLeftContext() const;
  size_t maxRightContext() const;
  size_t stateNameWidth() const;
  size_t stateIndexWidth() const;

  string inputAlphabet (int inputFlags = MachineDefaultInputFlags) const;
  string outputAlphabet() const;
  string inputDescriptionTable() const;

  map<InputSymbol,double> expectedBasesPerInputSymbol (const char* symbols = "01") const;

  Machine waitingMachine() const;  // convert to waiting machine
  vguard<State> decoderToposort (const string& inputAlphabet) const;  // topological sort by non-output transitions
};

#endif /* TRANSDUCER_INCLUDED */
