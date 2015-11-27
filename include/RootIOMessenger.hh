/* ------------------------------------------------------------------------
            |\___/|       
            )     (    
           =\     /=
             )===(
            /     \         CaTS: Calorimeter and Tracker Simulation
            |     |         Author: Hans Wenzel (Fermilab)
           /       \
           \       /
            \__  _/
              ( (
               ) )
              (_(
-------------------------------------------------------------------------*/
#ifndef RootIOMessenger_h
#define RootIOMessenger_h 1

#include "RootIO.hh"

#include "G4UIdirectory.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"

class RootIO;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;



class RootIOMessenger : public G4UImessenger
{
public:
  RootIOMessenger(RootIO*);
  ~RootIOMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  
  RootIO* pRootIO;
  G4UIdirectory* RootIODir;
  G4UIcmdWithAString* pFilenameCmd;
  G4UIcmdWithAnInteger* pBranchStatus1Cmd;
  G4UIcmdWithAnInteger* pBranchStatus2Cmd;
  G4UIcmdWithAnInteger* pBranchStatus3Cmd;
  G4UIcmdWithAnInteger* pBranchStatus4Cmd;
  G4UIcmdWithAnInteger* pBranchStatus5Cmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
