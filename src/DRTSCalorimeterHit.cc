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
#include "DRTSCalorimeterHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<DRTSCalorimeterHit> DRTSCalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTSCalorimeterHit::DRTSCalorimeterHit() {
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTSCalorimeterHit::DRTSCalorimeterHit(G4double ed, G4double eem, G4double enonem, G4double eobsb, G4int nc, G4ThreeVector p,
                                       G4double gt, G4double lt1, G4double lt2, G4double ts,
                                       G4double el, G4double v) {
    edep = ed;
    edepem = eem;
    edepnonem = enonem;
    eobsbirks = eobsb;
    nceren = nc;
    pos = p;
    globalTime = gt;
    localTime1 = lt1;
    localTime2 = lt2;
    dedx = el;
    vel = v;
    
    G4int slice = int(globalTime/ts);
    DRTimeSliceHit tslht = DRTimeSliceHit(slice,edep,edepem,edepnonem,nceren);
    fDRGlobalTimeSliceHitCollection.push_back(tslht);
    
    slice = int(localTime1/ts);
    tslht = DRTimeSliceHit(slice,edep,edepem,edepnonem,nceren);
    fDRLocalTime1SliceHitCollection.push_back(tslht);
    
    slice = int(localTime2/ts);
    tslht = DRTimeSliceHit(slice,edep,edepem,edepnonem,nceren);
    fDRLocalTime2SliceHitCollection.push_back(tslht);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTSCalorimeterHit::~DRTSCalorimeterHit() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTSCalorimeterHit::DRTSCalorimeterHit(const DRTSCalorimeterHit& right)
: G4VHit() {

    edep = right.edep;
    edepem = right.edepem;
    edepnonem = right.edepnonem;
    eobsbirks = right.eobsbirks;
    nceren = right.nceren;
    pos = right.pos;
    globalTime = right.globalTime;
    localTime1 = right.localTime1;
    localTime2 = right.localTime2;
    dedx = right.dedx;
    vel = right.vel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const DRTSCalorimeterHit& DRTSCalorimeterHit::operator=(const DRTSCalorimeterHit& right) {

    edep = right.edep;
    edepem = right.edepem;
    edepnonem = right.edepnonem;
    eobsbirks = right.eobsbirks;
    nceren = right.nceren;
    pos = right.pos;
    globalTime = right.globalTime;
    localTime1 = right.localTime1;
    localTime2 = right.localTime2;
    dedx = right.dedx;
    vel = right.vel;
    
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DRTSCalorimeterHit::operator==(const DRTSCalorimeterHit& right) const {
    return (this == &right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRTSCalorimeterHit::Draw() {
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

void DRTSCalorimeterHit::Print() {
    G4cout << "  energy deposit[MeV]: " << edep
            << "  em energy deposit[MeV]: " << edepem
            << "  non em energy deposit[MeV]: " << edepnonem
            << "  birks suppressed energy deposit[MeV]: " <<eobsbirks
            << "  Nr of cerenkov Photons " << nceren
            << "  position[mm]: " << pos
            << "  Global time of first energy deposit: " << globalTime 
            << "  Energy loss per distance traveled: " << dedx 
            << "  Velocity of protons: " << vel << G4endl;
    G4cout << "================================================================" << G4endl;
    for (std::vector<DRTimeSliceHit>::iterator it = fDRGlobalTimeSliceHitCollection.begin(); it != fDRGlobalTimeSliceHitCollection.end(); ++it) {
        (*it).Print();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRTSCalorimeterHit::AddDRTimeSliceHit(G4int type, G4int slc, G4double edp, G4double edpEM, G4double edpnonEM, G4int NCer) {
    DRTimeSliceHit TSHit = DRTimeSliceHit(slc, edp, edpEM, edpnonEM, NCer);
    if( type == 0 ) fDRGlobalTimeSliceHitCollection.push_back(TSHit);
    if( type == 1 ) fDRLocalTime1SliceHitCollection.push_back(TSHit);
    if( type == 2 ) fDRLocalTime2SliceHitCollection.push_back(TSHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRTSCalorimeterHit::AddDRTimeSliceHit(DRTimeSliceHit TSlHit, G4int type) {
    if( type == 0 ) fDRGlobalTimeSliceHitCollection.push_back(TSlHit);
    if( type == 1 ) fDRLocalTime1SliceHitCollection.push_back(TSlHit);
    if( type == 2 ) fDRLocalTime2SliceHitCollection.push_back(TSlHit);
}
