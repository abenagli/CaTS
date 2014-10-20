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

#include "EventActionMessenger.hh"

#include "EventAction.hh"



EventActionMessenger::EventActionMessenger(EventAction* evAct):
  eventAction(evAct)
{
  saveSeed = new G4UIcmdWithABool("/CaTS/random/saveSeed", this);
  saveSeed->SetGuidance("Save random seed");
  saveSeed->SetParameterName("saveSeed", true);
  saveSeed->SetDefaultValue(false);
  saveSeed->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



EventActionMessenger::~EventActionMessenger()
{
  delete saveSeed;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void EventActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == saveSeed )
  {
    eventAction -> SetSeedSave(saveSeed->GetNewBoolValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
