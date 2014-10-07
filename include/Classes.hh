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
#include "G4VHit.hh"
#include "CalorimeterHit.hh"
#include "DRCalorimeterHit.hh"
#include "DRTimeSliceHit.hh"
#include "DRTSCalorimeterHit.hh"
#include "DRTSCalorimeterHit2.hh"
#include "products.hh"
#include "Event.hh"
#include "PhotonHit.hh"
#include "RunHeader.hh"
#include "TrackerHit.hh"

#include "TH1F.h"

DRTimeSliceHit tslh;
products p;
Event e;
RunHeader rh;
std::map<G4String, G4int > pmult;
std::map<G4String, G4double> Emap; // Energy deposited by particle type
std::map<G4String, std::map<G4String, G4int> > Vpmult;
std::map<G4String, std::map<G4String, G4double> > VEmap;
std::map<G4int, std::map<G4String, G4int> > Vpmultt;
std::map<G4int, std::map<G4String, G4double> > VEmapt;
std::map<G4String, std::map<G4int, std::map<G4String, G4int> > > Vpmultt2;
std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > VEmapt2;
std::map<G4String, std::vector<CalorimeterHit*> > cmap;
std::map<G4String, std::vector<DRCalorimeterHit*> > drmap;
std::map<G4String, std::vector<DRTSCalorimeterHit*> > drtsmap;
std::map<G4String, std::vector<DRTSCalorimeterHit2*> > drts2map;
std::map<G4String, std::vector<G4VHit*> > hcmap; // map of Hit Collections
std::map<G4String, std::vector<PhotonHit*> > pmap;
std::map<G4String, std::vector<TrackerHit*> > tmap;
std::vector<G4String> vs;
std::vector<G4VHit*> vh;
std::vector<CalorimeterHit*> c;
std::vector<DRCalorimeterHit*> d;
std::vector<DRTimeSliceHit*> tslhvec;
std::vector<DRTSCalorimeterHit*> drts;
std::vector<DRTSCalorimeterHit2*> drts2;
std::vector<PhotonHit*> b;
std::vector<TrackerHit*> a;
std::map<G4String, products> mp;
std::map<G4String, std::map<G4String, products> > pm;
std::map<G4String, std::map<G4String, std::map<G4String, products> > >ppm;
std::map<G4int, std::map<G4String, std::map<G4String, products> > > pmt;
std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > ppmt;
#undef __G4String
