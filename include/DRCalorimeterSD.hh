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

#ifndef DRCalorimeterSD_h
#define DRCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "DRCalorimeterHit.hh"
#include "Cerenkov.hh"
#include <vector>
//#include "G4PhysicsOrderedFreeVector.hh"
#include "G4Material.hh" 
//#include "G4MaterialPropertyVector.hh"
class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DRCalorimeterSD : public G4VSensitiveDetector {
public:
    DRCalorimeterSD(G4String);
    ~DRCalorimeterSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);

private:
    DRCalorimeterHitsCollection* drcalorimeterCollection;
    Cerenkov* CerenGenerator;
    //std::vector<G4String> CAI;
    //std::vector<G4PhysicsOrderedFreeVector*> CerenkovAngleIntegrals;
    //std::vector<G4MaterialPropertyVector*> RefractionIndeces;
    //G4double GetAverageNumberOfPhotons(const G4double, const G4double, const G4int) const;
    G4int HCID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

