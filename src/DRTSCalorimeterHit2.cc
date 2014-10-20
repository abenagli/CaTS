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
#include "DRTSCalorimeterHit2.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<DRTSCalorimeterHit2> DRTSCalorimeterHit2Allocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTSCalorimeterHit2::DRTSCalorimeterHit2() {
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTSCalorimeterHit2::DRTSCalorimeterHit2(G4String pn, G4double de, G4double eobsb, G4int nc, G4ThreeVector p,
                                         G4int ts)
{
  particleName = pn;
  edep = de;
  eobsbirks = eobsb;
  nceren = nc;
  pos = p;
  timeSlice = ts;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTSCalorimeterHit2::~DRTSCalorimeterHit2() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DRTSCalorimeterHit2::DRTSCalorimeterHit2(const DRTSCalorimeterHit2& right):
  G4VHit()
{
  particleName = right.particleName;
  edep = right.edep;
  eobsbirks = right.eobsbirks;
  nceren = right.nceren;
  pos = right.pos;
  timeSlice = right.timeSlice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const DRTSCalorimeterHit2& DRTSCalorimeterHit2::operator=(const DRTSCalorimeterHit2& right)
{
  particleName = right.particleName;
  edep = right.edep;
  eobsbirks = right.eobsbirks;
  nceren = right.nceren;
  pos = right.pos;
  timeSlice = right.timeSlice;
  
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DRTSCalorimeterHit2::operator==(const DRTSCalorimeterHit2& right) const
{
  return (this == &right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DRTSCalorimeterHit2::Draw()
{
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

void DRTSCalorimeterHit2::Print()
{
  G4cout << "   particle name: " << particleName
         << "  energy deposit[MeV]: " << edep/CLHEP::MeV
         << "  birks suppressed energy deposit[MeV]: " <<eobsbirks/CLHEP::MeV
         << "  Nr of cerenkov Photons " << nceren
         << "  position[mm]: " << pos
         << "  Time slice of energy deposit: " << timeSlice
         << G4endl;
  G4cout << "================================================================" << G4endl;
}
