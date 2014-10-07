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
#ifndef RunAction_h
#define RunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Timer;
class G4Run;
class PrimaryGeneratorAction;

class RunAction : public G4UserRunAction {
public:
    RunAction(G4String fname, G4String pname, G4bool eo, G4bool es);
    ~RunAction();
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

private:
    G4Timer* timer;
    G4String gdmlFile;
    G4String PhysicsList;
    G4bool enableoptics;
    G4bool enablescint;
  
    static PrimaryGeneratorAction* pgA; // pointer to the particle source 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*RunAction_h*/
