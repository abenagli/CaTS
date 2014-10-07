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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "AnalysisMessenger.hh"
#include "Analysis.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AnalysisMessenger::AnalysisMessenger(Analysis* pA):pAnalysis(pA) {
    
    AnalysisDir = new G4UIdirectory("/CaTS/Analysis/");
    AnalysisDir->SetGuidance("CaTS commands to control analysis.");
    
    pFilenameCmd = new G4UIcmdWithAString("/CaTS/Analysis/Filename", this);
    pFilenameCmd->SetGuidance("change name of root file for outputting the hits");
    pFilenameCmd->SetParameterName("choice", false);
    pFilenameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AnalysisMessenger::~AnalysisMessenger() {
    delete pFilenameCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
    if (command == pFilenameCmd) {
        pAnalysis->SetFileName(newValue);
    }
}
