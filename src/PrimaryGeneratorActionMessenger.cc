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

#include "PrimaryGeneratorActionMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
// clhep
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/Randomize.h"

using namespace CLHEP;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorActionMessenger::PrimaryGeneratorActionMessenger(PrimaryGeneratorAction* CaTSGun)
  : CaTSAction(CaTSGun)
{
  genCmd = new G4UIcmdWithAString("/CaTS/generator", this);
  genCmd->SetGuidance("Select primary generator.");
  genCmd->SetGuidance("Available generators : particleGun,GPS,HEPEvt");
  genCmd->SetParameterName("generator", true);
  genCmd->SetDefaultValue("particleGun");
  genCmd->SetCandidates("particleGun GPS HEPEvt");
  
  randomSeed = new G4UIcmdWithAnInteger("/CaTS/random/randomSeed", this);
  randomSeed->SetGuidance("Set random seed");
  randomSeed->SetParameterName("randomSeed", false);
  randomSeed->SetDefaultValue(0);
  randomSeed->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



PrimaryGeneratorActionMessenger::~PrimaryGeneratorActionMessenger()
{
  delete genCmd;
  delete randomSeed;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void PrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == randomSeed )
  {
    G4int seed = randomSeed->GetNewIntValue(newValue);
    
    // seed it
    HepRandom::setTheSeed(seed);
    
    G4cout << "set random seed: " << seed << G4endl;
  }
  
  else if (command == genCmd)
  {
    CaTSAction->SetGenerator(newValue);
    G4cout << "setting new primary Event Generator to: " << newValue << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

