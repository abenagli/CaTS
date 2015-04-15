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

#include "RootIOMessenger.hh"



RootIOMessenger::RootIOMessenger(RootIO* pR):
  pRootIO(pR)
{
  RootIODir = new G4UIdirectory("/CaTS/RootIO/");
  RootIODir->SetGuidance("CaTS RootIO control commands.");
  
  pFilenameCmd = new G4UIcmdWithAString("/CaTS/RootIO/Filename", this);
  pFilenameCmd->SetGuidance("change name of root file for outputting the hits");
  pFilenameCmd->SetParameterName("choice", false);
  pFilenameCmd->SetDefaultValue("hits.root");
  pFilenameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  pBranchStatus1Cmd = new G4UIcmdWithAnInteger("/CaTS/RootIO/MapByPosBranch", this);
  pBranchStatus1Cmd->SetGuidance("activate / deactivate MapByPos branch");
  pBranchStatus1Cmd->SetParameterName("status", false);
  pBranchStatus1Cmd->SetDefaultValue(1);
  pBranchStatus1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  pBranchStatus2Cmd = new G4UIcmdWithAnInteger("/CaTS/RootIO/HCMapBranch", this);
  pBranchStatus2Cmd->SetGuidance("activate / deactivate HCMap branch");
  pBranchStatus2Cmd->SetParameterName("status", false);
  pBranchStatus2Cmd->SetDefaultValue(1);
  pBranchStatus2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



RootIOMessenger::~RootIOMessenger()
{
  delete pFilenameCmd;
  delete pBranchStatus1Cmd;
  delete pBranchStatus2Cmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void RootIOMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == pFilenameCmd ){
    pRootIO->SetFileName(newValue);
  }
  if( command == pBranchStatus1Cmd ){
    pRootIO -> SetBranchStatus("m_Edep_byPos",pBranchStatus1Cmd->GetNewIntValue(newValue));
    pRootIO -> SetBranchStatus("m_Eobs_byPos",pBranchStatus1Cmd->GetNewIntValue(newValue));
    pRootIO -> SetBranchStatus("m_NCeren_byPos",pBranchStatus1Cmd->GetNewIntValue(newValue));
  }
  if( command == pBranchStatus2Cmd ){
    pRootIO -> SetBranchStatus("HCMap",pBranchStatus1Cmd->GetNewIntValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
