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
#include "DRTimeSliceHit.hh"
//#include "G4UnitsTable.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTimeSliceHit::DRTimeSliceHit() {
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTimeSliceHit::DRTimeSliceHit(G4int sl, G4double ed, G4double eem, G4double enonem, G4int nc) {
    slice = sl;
    edep = ed;
    edepem = eem;
    edepnonem = enonem;
    nceren = nc;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTimeSliceHit::~DRTimeSliceHit() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTimeSliceHit::DRTimeSliceHit(const DRTimeSliceHit& right)
 {
    slice = right.slice;
    edep = right.edep;
    edepem = right.edepem;
    edepnonem = right.edepnonem;
    nceren = right.nceren;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const DRTimeSliceHit& DRTimeSliceHit::operator=(const DRTimeSliceHit& right) {
    slice = right.slice;
    edep = right.edep;
    edepem = right.edepem;
    edepnonem = right.edepnonem;
    nceren = right.nceren;
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DRTimeSliceHit::operator==(const DRTimeSliceHit& right) const {
    return (this == &right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRTimeSliceHit::Print() {
    G4cout  << " Time slice Nr:  " << slice
            << "  energy deposit[MeV]: " << edep
            << "  em energy deposit[MeV]: " << edepem
            << "  non em energy deposit[MeV]: " << edepnonem
            << "  Nr of Cerenkov Photons " << nceren << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

