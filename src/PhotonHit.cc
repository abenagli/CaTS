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

#include "PhotonHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<PhotonHit> PhotonHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhotonHit::PhotonHit() {
}

PhotonHit::PhotonHit(G4int id, G4double e, G4ThreeVector p, G4double t) {

    processID = id;
    edep = e;
    pos = p;
    time = t;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhotonHit::~PhotonHit() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhotonHit::PhotonHit(const PhotonHit& right)
: G4VHit() {
    processID = right.processID;
    edep = right.edep;
    pos = right.pos;
    time = right.time;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const PhotonHit& PhotonHit::operator=(const PhotonHit& right) {
    processID = right.processID;
    edep = right.edep;
    pos = right.pos;
    time = right.time;
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int PhotonHit::operator==(const PhotonHit& right) const {
    return (this == &right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhotonHit::Draw() {
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

void PhotonHit::Print() {
    G4cout << "  processID: " << processID
            << "  energy deposit[MeV]: " << edep
            << "  position[mm]: " << pos
            << "  time[nsec]: " << time
            << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

