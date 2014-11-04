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
#include "MyMainFrame.hh"

Event e;
RunHeader rh;
products p;
G4VHit h;
DRTSCalorimeterHit2* drtsch;

MyMainFrame* mmf;
GroupBox* gb;

// energy by detector
std::map<G4String,G4float> sfMap;
// energy by detector and time
std::map<G4int,G4float> ifMap;
std::map<G4String,std::map<G4int,G4float> > sifMap;
std::map<G4String,std::map<G4String,std::map<G4int,G4float> > > ssifMap;

// energy by detector and particle
std::map<G4String,std::map<G4String,G4float> > ssfMap;
// energy by detector and particle and time
std::map<G4int,std::map<G4String,G4float> > isfMap;
std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > sisfMap;
std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > > ssisfMap;

// process multiplicity
std::map<G4String,G4int> siMap;
std::map<G4String,std::map<G4String,G4int> > ssiMap;
// process multiplicity and time
std::map<G4int,std::map<G4String,G4int> > isiMap;
std::map<G4String,std::map<G4int,std::map<G4String,G4int> > > sisiMap;
std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4int> > > > ssisiMap;

// processes
std::map<G4String,products> sp;
std::map<G4String,std::map<G4String,products> > ssp;
std::map<G4String,std::map<G4String,std::map<G4String,products> > > sssp;
// processes and time
std::map<G4int,std::map<G4String,std::map<G4String,products> > > issp;
std::map<G4String,std::map<G4int,std::map<G4String,std::map<G4String,products> > > > sissp;
std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,std::map<G4String,products> > > > > ssissp;

// hits
std::vector<G4VHit*> vh;
std::map<G4ThreeVector,std::vector<G4VHit*> > vvh;
std::map<G4String,std::map<G4ThreeVector,std::vector<G4VHit*> > > svvh;
std::map<G4String,std::map<G4String,std::map<G4ThreeVector,std::vector<G4VHit*> > > > ssvvh;
std::vector<DRTSCalorimeterHit2*> vdrtsch;
std::map<G4ThreeVector,std::vector<DRTSCalorimeterHit2*> > vvdrtsch;
std::map<G4String,std::map<G4ThreeVector,std::vector<DRTSCalorimeterHit2*> > > svvdrtsch;
std::map<G4String,std::map<G4String,std::map<G4ThreeVector,std::vector<DRTSCalorimeterHit2*> > > > ssvvdrtsch;

#undef __G4String
