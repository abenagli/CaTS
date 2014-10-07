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

#ifndef PhotonSD_h
#define PhotonSD_h 1

#include "G4VSensitiveDetector.hh"
#include "PhotonHit.hh"
#include "globals.hh"

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhotonSD : public G4VSensitiveDetector {
public:
    PhotonSD(G4String);
    ~PhotonSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);

private:
    PhotonHitsCollection* photonCollection;
    G4int HCID;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

