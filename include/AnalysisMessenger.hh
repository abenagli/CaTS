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
#ifndef AnalysisMessenger_h
#define AnalysisMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class Analysis;
class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class AnalysisMessenger : public G4UImessenger {
public:
    AnalysisMessenger(Analysis*);
    ~AnalysisMessenger();

    void SetNewValue(G4UIcommand*, G4String);
  
private:
    
    Analysis* pAnalysis;
    G4UIdirectory* AnalysisDir;
    G4UIcmdWithAString* pFilenameCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

