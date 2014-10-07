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

#include "G4UserEventAction.hh"
#include <vector>
class Event;
class G4Event;
class G4string;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction {
private:
    Event* CaTSEvt;
    static EventAction* instance;
public:
    EventAction();
    ~EventAction();

    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    Event* GetEvent(){return CaTSEvt;}
    static EventAction* GetInstance();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


