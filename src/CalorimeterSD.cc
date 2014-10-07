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

#include "CalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::CalorimeterSD(G4String name)
: G4VSensitiveDetector(name) {
    G4String HCname = name + "_HC";
    collectionName.insert(HCname);
    G4cout << collectionName.size() << "   CalorimeterSD name:  " << name << " collection Name: " << HCname << G4endl;
    HCID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::~CalorimeterSD() {

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::Initialize(G4HCofThisEvent* HCE) {

    calorimeterCollection = new CalorimeterHitsCollection
            (SensitiveDetectorName, collectionName[0]);
    if (HCID < 0) {
        G4cout << "CalorimeterSD::Initialize:  " << SensitiveDetectorName << "   " << collectionName[0] << G4endl;
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

    }
    HCE->AddHitsCollection(HCID, calorimeterCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CalorimeterSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    G4double edep = aStep->GetTotalEnergyDeposit() / CLHEP::MeV;
    if (edep == 0.) return false;
    const G4double time = aStep->GetPreStepPoint()->GetGlobalTime() / CLHEP::ns;
    const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
    const G4ThreeVector cellpos = touch->GetTranslation();
    G4Track* theTrack = aStep->GetTrack();
    G4String particleType = theTrack->GetDefinition()->GetParticleName();
    for (G4int j = 0; j < calorimeterCollection->entries(); j++) {
        CalorimeterHit* aPreviousHit = (*calorimeterCollection)[j];
        if (cellpos == aPreviousHit->GetPos()) {
            aPreviousHit->SetEdep(aStep->GetTotalEnergyDeposit() + aPreviousHit->GetEdep());
            if ((particleType == "e+") || (particleType == "gamma") || (particleType == "e-")) {
                aPreviousHit->SetEdepEM(edep + aPreviousHit->GetEdepEM());
            } else {
                aPreviousHit->SetEdepnonEM(edep + aPreviousHit->GetEdepnonEM());
            }
            return true;
        }
    }
    CalorimeterHit* newHit = new CalorimeterHit();
    newHit->SetEdep(edep);
    newHit->SetPos(cellpos);
    newHit->SetTime(time);
    if ((particleType == "e+") || (particleType == "gamma") || (particleType == "e-")) {
        newHit->SetEdepEM(edep);
        newHit->SetEdepnonEM(0.0);
    } else {
        newHit->SetEdepnonEM(edep);
        newHit->SetEdepEM(0.0);
    }

    calorimeterCollection->insert(newHit);
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

