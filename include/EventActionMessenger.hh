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
#ifndef EventActionMessenger_h
#define EventActionMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithABool.hh"

class EventAction;



class EventActionMessenger : public G4UImessenger
{
public:
  EventActionMessenger(EventAction*);
  ~EventActionMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  EventAction* eventAction;
  G4UIcmdWithABool* saveSeed;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
