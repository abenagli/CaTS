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

#include "CalorimeterHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<CalorimeterHit> CalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterHit::CalorimeterHit() {
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
CalorimeterHit::CalorimeterHit(G4double ed, G4double eem, G4double enonem, G4ThreeVector p, G4double t) {
    edep = ed;
    edepem = eem;
    edepnonem = enonem;
    pos = p;
    time = t;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterHit::~CalorimeterHit() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterHit::CalorimeterHit(const CalorimeterHit& right)
: G4VHit() {

    edep = right.edep;
    edepem = right.edepem;
    edepnonem = right.edepnonem;
    pos = right.pos;
    time = right.time;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const CalorimeterHit& CalorimeterHit::operator=(const CalorimeterHit& right) {

    edep = right.edep;
    edepem = right.edepem;
    edepnonem = right.edepnonem;
    pos = right.pos;
    time = right.time;
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int CalorimeterHit::operator==(const CalorimeterHit& right) const {
    return (this == &right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterHit::Draw() {
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if (pVVisManager) {
        G4Circle circle(pos);
        circle.SetScreenSize(2.);
        circle.SetFillStyle(G4Circle::filled);
        G4Colour colour(1., 0., 0.);
        G4VisAttributes attribs(colour);
        circle.SetVisAttributes(attribs);
        pVVisManager->Draw(circle);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterHit::Print() {

    G4cout << "  energy deposit[MeV]: " << edep
            << "  em energy deposit[MeV]: " << edepem
            << "  non em energy deposit[MeV]: " << edepnonem
            << "  position[mm]: " << pos 
            << " Time of first energy deposit: "<< time<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

