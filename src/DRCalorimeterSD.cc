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
#include "DRCalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Poisson.hh"
#include "G4VVisManager.hh"

#include "RootIO.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRCalorimeterSD::DRCalorimeterSD(G4String name)
: G4VSensitiveDetector(name) {
    CerenGenerator = new(Cerenkov);
    G4String HCname = name + "_HC";
    collectionName.insert(HCname);
    G4cout << collectionName.size() << "   DRCalorimeterSD name:  " << name << " collection Name: " << HCname << G4endl;
    HCID = -1;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRCalorimeterSD::~DRCalorimeterSD() {

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRCalorimeterSD::Initialize(G4HCofThisEvent* HCE) {
    drcalorimeterCollection = new DRCalorimeterHitsCollection
            (SensitiveDetectorName, collectionName[0]);
    if (HCID < 0) {
        G4cout << "DRCalorimeterSD::Initialize:  " << SensitiveDetectorName << "   " << collectionName[0] << G4endl;
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

    }
    HCE->AddHitsCollection(HCID, drcalorimeterCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DRCalorimeterSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    G4double edep = aStep->GetTotalEnergyDeposit() / CLHEP::MeV;
    if (edep == 0.) return false;
    const G4double time = aStep->GetPreStepPoint()->GetGlobalTime() / CLHEP::ns;
    const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
    G4String thematerial = touch->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName();
    G4int NCerenPhotons = 0;
    G4Track* theTrack = aStep->GetTrack();
    const G4double charge = theTrack->GetDefinition()->GetPDGCharge();
    G4String particleType = theTrack->GetDefinition()->GetParticleName();
    G4StepPoint* pPreStepPoint = aStep->GetPreStepPoint();
    G4StepPoint* pPostStepPoint = aStep->GetPostStepPoint();
    G4double beta = 0.5 * (pPreStepPoint ->GetBeta() + pPostStepPoint->GetBeta());
    G4double MeanNumberOfPhotons = CerenGenerator->GetAverageNumberOfPhotons(charge, beta, thematerial);
    if (MeanNumberOfPhotons > 0.0) {
        G4double step_length = aStep->GetStepLength();
        MeanNumberOfPhotons = MeanNumberOfPhotons * step_length;
        NCerenPhotons = (G4int) G4Poisson(MeanNumberOfPhotons);
    } else {
        NCerenPhotons = 0;
    }

    const G4ThreeVector cellpos = touch->GetTranslation();
    for (G4int j = 0; j < drcalorimeterCollection->entries(); j++) {
        DRCalorimeterHit* aPreviousHit = (*drcalorimeterCollection)[j];
        if (cellpos == aPreviousHit->GetPos()) {
            aPreviousHit->SetEdep(edep + aPreviousHit->GetEdep());
            aPreviousHit->SetNCeren(NCerenPhotons + aPreviousHit->GetNCeren());
            if ((particleType == "e+") || (particleType == "gamma") || (particleType == "e-")) {
                aPreviousHit->SetEdepEM(edep + aPreviousHit->GetEdepEM());
            } else {
                aPreviousHit->SetEdepnonEM(edep + aPreviousHit->GetEdepnonEM());
            }
            return true;
        }
    }

    DRCalorimeterHit* newHit = new DRCalorimeterHit();
    newHit->SetEdep(edep);
    newHit->SetPos(cellpos);
    newHit->SetNCeren(NCerenPhotons);
    newHit->SetTime(time);
    if ((particleType == "e+") || (particleType == "gamma") || (particleType == "e-")) {
        newHit->SetEdepEM(edep);
        newHit->SetEdepnonEM(0.0);
    } else {
        newHit->SetEdepnonEM(edep);
        newHit->SetEdepEM(0.0);
    }
    drcalorimeterCollection->insert(newHit);

    return true;
}
