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

#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"

#include "RootIO.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(G4String name)
: G4VSensitiveDetector(name) {
    G4String HCname;
    collectionName.insert(HCname = "trackerCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::~TrackerSD() {
    RootIO::GetInstance()->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent* HCE) {
    trackerCollection = new TrackerHitsCollection
            (SensitiveDetectorName, collectionName[0]);
    static G4int HCID = -1;
    if (HCID < 0) {
        G4cout << "TrackerSD::Initialize:  " << SensitiveDetectorName << "   " << collectionName[0] << G4endl;
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
    HCE->AddHitsCollection(HCID, trackerCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    G4double edep = aStep->GetTotalEnergyDeposit();
    if (edep == 0.) return false;

    TrackerHit* newHit = new TrackerHit();
    newHit->SetEdep(edep);
    newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
    trackerCollection->insert(newHit);
    return true;
}

