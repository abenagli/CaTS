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

#include "G4UImessenger.hh"
#include "globals.hh"

class RootIO;
class G4UIdirectory;
class G4UIcmdWithAString;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RootIOMessenger : public G4UImessenger {
public:
    RootIOMessenger(RootIO*);
    ~RootIOMessenger();

    void SetNewValue(G4UIcommand*, G4String);

private:
    
    RootIO* pRootIO;
    G4UIdirectory*      RootIODir;
    G4UIcmdWithAString*   pFilenameCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

