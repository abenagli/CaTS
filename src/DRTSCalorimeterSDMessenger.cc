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

#include "DRTSCalorimeterSDMessenger.hh"
#include "DRTSCalorimeterSD.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTSCalorimeterSDMessenger::DRTSCalorimeterSDMessenger(DRTSCalorimeterSD* pA):
  pDRTSCalorimeterSD(pA)
{
  G4String mess_directory = "/CaTS/SensitiveDetectors/" + pA->GetName() + "/";
  AnalysisDir = new G4UIdirectory(mess_directory);
  AnalysisDir->SetGuidance("CaTS commands to control setting of sensitive Detectors.");
  
  G4String c1dir = mess_directory + "SetBirksC1";
  G4String c2dir = mess_directory + "SetBirksC2";
  psetc1Cmd = new G4UIcmdWithADouble(c1dir, this);
  psetc1Cmd->SetGuidance("change value of C1 Birks const");
  psetc1Cmd->SetParameterName("choice", false);
  psetc1Cmd->SetDefaultValue(1.29e-2);
  psetc1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  psetc2Cmd = new G4UIcmdWithADouble(c2dir, this);
  psetc2Cmd->SetGuidance("change value of C2 Birks const");
  psetc2Cmd->SetParameterName("choice", false);
  psetc2Cmd->SetDefaultValue(9.59e-6);
  psetc2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTSCalorimeterSDMessenger::~DRTSCalorimeterSDMessenger() {
  delete psetc1Cmd;
  delete psetc2Cmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRTSCalorimeterSDMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == psetc1Cmd) {
    pDRTSCalorimeterSD->SetBirksc1(psetc1Cmd->GetNewDoubleValue(newValue));
  }
  else if (command == psetc2Cmd) {
    pDRTSCalorimeterSD->SetBirksc2(psetc1Cmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String DRTSCalorimeterSDMessenger::GetCurrentValue(G4UIcommand* command) {
  G4String cv;
  if (command == psetc1Cmd) {
    cv = psetc1Cmd->ConvertToString(pDRTSCalorimeterSD->GetBirksc1());
  }
  else if (command == psetc2Cmd) {
    cv = psetc2Cmd->ConvertToString(pDRTSCalorimeterSD->GetBirksc2());
  }
  return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
