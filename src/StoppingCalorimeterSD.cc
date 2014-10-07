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

#include "StoppingCalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"

#include "RootIO.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StoppingCalorimeterSD::StoppingCalorimeterSD(G4String name)
: G4VSensitiveDetector(name) {
    G4String HCname;
    collectionName.insert(HCname = "stpcalorimeterCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StoppingCalorimeterSD::~StoppingCalorimeterSD() {
    RootIO::GetInstance()->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StoppingCalorimeterSD::Initialize(G4HCofThisEvent* HCE) {
    calorimeterCollection = new CalorimeterHitsCollection
            (SensitiveDetectorName, collectionName[0]);
    static G4int HCID = -1;
    if (HCID < 0) {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
        HCE->AddHitsCollection(HCID, calorimeterCollection);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool StoppingCalorimeterSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    G4double edep = aStep->GetTotalEnergyDeposit();
    G4Track* theTrack = aStep->GetTrack();
    G4double Energy = theTrack->GetTotalEnergy() + edep;
    const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
    const G4ThreeVector cellpos = touch->GetTranslation();
    for (G4int j = 0; j < calorimeterCollection->entries(); j++) {
        CalorimeterHit* aPreviousHit = (*calorimeterCollection)[j];
        if (cellpos == aPreviousHit->GetPos()) {
            aPreviousHit->SetEdep(Energy + aPreviousHit->GetEdep());
            if (theTrack->GetParentID() != 0) {
                theTrack->SetTrackStatus(fStopAndKill);
            }
            return true;
        }
    }
    CalorimeterHit* newHit = new CalorimeterHit();
    newHit->SetEdep(edep);
    newHit->SetPos(cellpos);
    calorimeterCollection->insert(newHit);
    if (theTrack->GetParentID() != 0) {
        theTrack->SetTrackStatus(fStopAndKill);
    }
    //    newHit->Print();
    newHit->Draw();
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


