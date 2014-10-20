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
  G4String mess_directory = "/CaTS/SensitiveDetectors/" + pA->GetName()+"/";    
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
  
  
  mess_directory = "/CaTS/SensitiveDetectors/";    
  
  G4String c3dir = mess_directory + "TimeSliceSizeLow";
  psetc3Cmd = new G4UIcmdWithADoubleAndUnit(c3dir, this);
  psetc3Cmd->SetGuidance("change size of the time slice bin");
  psetc3Cmd->SetParameterName("choice", false);
  psetc3Cmd->SetDefaultValue(0.01*CLHEP::ns);
  psetc3Cmd->SetUnitCategory("Time");
  psetc3Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  G4String c3dirMin = mess_directory + "MinTimeLow";
  psetc3MinCmd = new G4UIcmdWithADoubleAndUnit(c3dirMin, this);
  psetc3MinCmd->SetGuidance("change the minimum time");
  psetc3MinCmd->SetParameterName("choice", false);
  psetc3MinCmd->SetDefaultValue(0.*CLHEP::ns);
  psetc3MinCmd->SetUnitCategory("Time");
  psetc3MinCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  G4String c3dirMax = mess_directory + "MaxTimeLow";
  psetc3MaxCmd = new G4UIcmdWithADoubleAndUnit(c3dirMax, this);
  psetc3MaxCmd->SetGuidance("change the maximum time");
  psetc3MaxCmd->SetParameterName("choice", false);
  psetc3MaxCmd->SetDefaultValue(0.*CLHEP::ns);
  psetc3MaxCmd->SetUnitCategory("Time");
  psetc3MaxCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  G4String c4dir = mess_directory + "TimeSliceSizeMed";
  psetc4Cmd = new G4UIcmdWithADoubleAndUnit(c4dir, this);
  psetc4Cmd->SetGuidance("change size of the time slice bin");
  psetc4Cmd->SetParameterName("choice", false);
  psetc4Cmd->SetDefaultValue(1.*CLHEP::ns);
  psetc4Cmd->SetUnitCategory("Time");
  psetc4Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  G4String c4dirMin = mess_directory + "MinTimeMed";
  psetc4MinCmd = new G4UIcmdWithADoubleAndUnit(c4dirMin, this);
  psetc4MinCmd->SetGuidance("change the minimum time");
  psetc4MinCmd->SetParameterName("choice", false);
  psetc4MinCmd->SetDefaultValue(0.*CLHEP::ns);
  psetc4MinCmd->SetUnitCategory("Time");
  psetc4MinCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  G4String c4dirMax = mess_directory + "MaxTimeMed";
  psetc4MaxCmd = new G4UIcmdWithADoubleAndUnit(c4dirMax, this);
  psetc4MaxCmd->SetGuidance("change the maximum time");
  psetc4MaxCmd->SetParameterName("choice", false);
  psetc4MaxCmd->SetDefaultValue(0.*CLHEP::ns);
  psetc4MaxCmd->SetUnitCategory("Time");
  psetc4MaxCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  G4String c5dir = mess_directory + "TimeSliceSizeHig";
  psetc5Cmd = new G4UIcmdWithADoubleAndUnit(c5dir, this);
  psetc5Cmd->SetGuidance("change size of the time slice bin");
  psetc5Cmd->SetParameterName("choice", false);
  psetc5Cmd->SetDefaultValue(100.*CLHEP::ns);
  psetc5Cmd->SetUnitCategory("Time");
  psetc5Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  G4String c5dirMin = mess_directory + "MinTimeHig";
  psetc5MinCmd = new G4UIcmdWithADoubleAndUnit(c5dirMin, this);
  psetc5MinCmd->SetGuidance("change the minimum time");
  psetc5MinCmd->SetParameterName("choice", false);
  psetc5MinCmd->SetDefaultValue(0.*CLHEP::ns);
  psetc5MinCmd->SetUnitCategory("Time");
  psetc5MinCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  G4String c5dirMax = mess_directory + "MaxTimeHig";
  psetc5MaxCmd = new G4UIcmdWithADoubleAndUnit(c5dirMax, this);
  psetc5MaxCmd->SetGuidance("change the maximum time");
  psetc5MaxCmd->SetParameterName("choice", false);
  psetc5MaxCmd->SetDefaultValue(0.*CLHEP::ns);
  psetc5MaxCmd->SetUnitCategory("Time");
  psetc5MaxCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTSCalorimeterSDMessenger::~DRTSCalorimeterSDMessenger() {
  delete psetc1Cmd;
  delete psetc2Cmd;
  delete psetc3Cmd;
  delete psetc3MinCmd;
  delete psetc3MaxCmd;
  delete psetc4Cmd;
  delete psetc4MinCmd;
  delete psetc4MaxCmd;
  delete psetc5Cmd;
  delete psetc5MinCmd;
  delete psetc5MaxCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRTSCalorimeterSDMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == psetc1Cmd) {
    pDRTSCalorimeterSD->SetBirksc1(psetc1Cmd->GetNewDoubleValue(newValue));
  }
  else if (command == psetc2Cmd) {
    pDRTSCalorimeterSD->SetBirksc2(psetc1Cmd->GetNewDoubleValue(newValue));
  }
  else if (command == psetc3Cmd) {
    pDRTSCalorimeterSD->SetTimeSliceSizeLow(psetc3Cmd->GetNewDoubleValue(newValue));
  }
  else if (command == psetc3MinCmd) {
    pDRTSCalorimeterSD->SetMinTimeLow(psetc3MinCmd->GetNewDoubleValue(newValue));
  }
  else if (command == psetc3MaxCmd) {
    pDRTSCalorimeterSD->SetMaxTimeLow(psetc3MaxCmd->GetNewDoubleValue(newValue));
  }
  else if (command == psetc4Cmd) {
    pDRTSCalorimeterSD->SetTimeSliceSizeMed(psetc4Cmd->GetNewDoubleValue(newValue));
  }
  else if (command == psetc4MinCmd) {
    pDRTSCalorimeterSD->SetMinTimeMed(psetc4MinCmd->GetNewDoubleValue(newValue));
  }
  else if (command == psetc4MaxCmd) {
    pDRTSCalorimeterSD->SetMaxTimeMed(psetc4MaxCmd->GetNewDoubleValue(newValue));
  }  
  else if (command == psetc5Cmd) {
    pDRTSCalorimeterSD->SetTimeSliceSizeHig(psetc5Cmd->GetNewDoubleValue(newValue));
  }
  else if (command == psetc5MinCmd) {
    pDRTSCalorimeterSD->SetMinTimeHig(psetc5MinCmd->GetNewDoubleValue(newValue));
  }
  else if (command == psetc5MaxCmd) {
    pDRTSCalorimeterSD->SetMaxTimeHig(psetc5MaxCmd->GetNewDoubleValue(newValue));
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
  else if (command == psetc3Cmd) {
    cv = psetc3Cmd->ConvertToString(pDRTSCalorimeterSD->GetTimeSliceSizeLow());
  }
  else if (command == psetc3MinCmd) {
    cv = psetc3MinCmd->ConvertToString(pDRTSCalorimeterSD->GetMinTimeLow());
  }
  else if (command == psetc3MaxCmd) {
    cv = psetc3MaxCmd->ConvertToString(pDRTSCalorimeterSD->GetMaxTimeLow());
  }
  else if (command == psetc4Cmd) {
    cv = psetc4Cmd->ConvertToString(pDRTSCalorimeterSD->GetTimeSliceSizeMed());
  }
  else if (command == psetc4MinCmd) {
    cv = psetc4MinCmd->ConvertToString(pDRTSCalorimeterSD->GetMinTimeMed());
  }
  else if (command == psetc4MaxCmd) {
    cv = psetc4MaxCmd->ConvertToString(pDRTSCalorimeterSD->GetMaxTimeMed());
  }
  else if (command == psetc5Cmd) {
    cv = psetc5Cmd->ConvertToString(pDRTSCalorimeterSD->GetTimeSliceSizeHig());
  }
  else if (command == psetc5MinCmd) {
    cv = psetc5MinCmd->ConvertToString(pDRTSCalorimeterSD->GetMinTimeHig());
  }
  else if (command == psetc5MaxCmd) {
    cv = psetc5MaxCmd->ConvertToString(pDRTSCalorimeterSD->GetMaxTimeHig());
  }
  return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
