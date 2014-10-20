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

#ifndef EventAction_h
#define EventAction_h 1

#include "EventActionMessenger.hh"

#include "G4UserEventAction.hh"

#include <vector>

class Event;
class G4Event;
class G4string;



class EventAction : public G4UserEventAction
{
private:
  Event* CaTSEvt;
  EventActionMessenger* eventMessenger;
  G4bool seedSave;
  
  static EventAction* instance;
  
  
public:
  EventAction();
  ~EventAction();
  
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  
  void SetSeedSave(const G4bool& val) { seedSave = val; };
  
  Event* GetEvent() { return CaTSEvt; };
  
  static EventAction* GetInstance();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
