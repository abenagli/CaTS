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
#include "DRCalorimeterHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<DRCalorimeterHit> DRCalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRCalorimeterHit::DRCalorimeterHit() {
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DRCalorimeterHit::DRCalorimeterHit(G4double ed, G4double eem, G4double enonem, G4int nc, G4ThreeVector p, G4double t) {
    edep = ed;
    edepem = eem;
    edepnonem = enonem;
    nceren = nc;
    pos = p;
    time = t;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRCalorimeterHit::~DRCalorimeterHit() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRCalorimeterHit::DRCalorimeterHit(const DRCalorimeterHit& right)
: G4VHit() {

    edep = right.edep;
    edepem = right.edepem;
    edepnonem = right.edepnonem;
    nceren = right.nceren;
    pos = right.pos;
    time = right.time;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const DRCalorimeterHit& DRCalorimeterHit::operator=(const DRCalorimeterHit& right) {

    edep = right.edep;
    edepem = right.edepem;
    edepnonem = right.edepnonem;
    nceren = right.nceren;
    pos = right.pos;
    time = right.time;

    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DRCalorimeterHit::operator==(const DRCalorimeterHit& right) const {
    return (this == &right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRCalorimeterHit::Draw() {
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if (pVVisManager) {
        G4Circle circle(pos);
        circle.SetScreenSize(4.);
        circle.SetFillStyle(G4Circle::filled);
        if (nceren > 0) {
            G4Colour colour(0., 0., 1.);
            G4VisAttributes attribs(colour);
            circle.SetVisAttributes(attribs);
            pVVisManager->Draw(circle);
        } else {
            G4Colour colour(1., 0., 0.);
            G4VisAttributes attribs(colour);
            circle.SetVisAttributes(attribs);
            pVVisManager->Draw(circle);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRCalorimeterHit::Print() {
    G4cout << "  energy deposit[MeV]: " << edep
            << "  em energy deposit[MeV]: " << edepem
            << "  non em energy deposit[MeV]: " << edepnonem
            << "  Nr of cerenkov Photons " << nceren
            << "  position[mm]: " << pos 
            << " Time of first energy deposit: "<< time<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

