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
#include <vector>
#include <map>
#include <sstream>
#include <utility>

#include "G4RunManager.hh"
#include "EventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4HCofThisEvent.hh"
#ifdef G4ANALYSIS_USE
#include "Analysis.hh"
#endif
#include "RootIO.hh"
#include "Event.hh"
#include "TrackerHit.hh"
#include "CalorimeterHit.hh"
#include "DRCalorimeterHit.hh"
#include "DRTSCalorimeterHit.hh"
#include "DRTSCalorimeterHit2.hh"
#include "PhotonHit.hh"
#include "TrackerHit.hh"

EventAction* EventAction::instance = 0;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

EventAction* EventAction::GetInstance() {
    if (instance == 0) instance = new EventAction();
    return instance;
}

EventAction::EventAction() {
    if (instance == 0) {
        CaTSEvt = new Event();
    }
    instance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt) {
#ifdef G4ANALYSIS_USE
    Analysis* analysis = Analysis::getInstance();
//    analysis->BeginOfEvent(evt->GetEventID());
#endif  
    CaTSEvt->SetEventNr(evt->GetEventID());
    long seed = CLHEP::HepRandom::getTheSeed();
    G4cout << "seed: " << seed << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*evt)
{
  G4RunManager* runManager = G4RunManager::GetRunManager();
  runManager -> rndmSaveThisEvent();
  
  evt->Print();
#ifdef G4ANALYSIS_USE
  Analysis* analysis = Analysis::getInstance();
#endif
  //    Event* EventContainer = new Event();
  //    EventContainer->Set(1, evt->GetEventID());
#ifdef G4ANALYSIS_USE
  //    EventContainer->SetPi0Energy(analysis->Pi0Energy);
#else
  //      EventContainer->SetPi0Energy(0);
#endif
  std::map<G4String, std::vector<G4VHit* > >* hcmap = CaTSEvt->GetHCMap();
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  std::vector<G4VHit*> hitsVector;
  //G4cout << "Number of collections:  " << HCE->GetNumberOfCollections()<<G4endl;
  for (int i = 0; i < HCE->GetNumberOfCollections(); i++)
  {
    hitsVector.clear();
    G4VHitsCollection* hc = HCE->GetHC(i);
    G4String hcname = hc->GetName();
    std::vector<std::string> y = split(hcname, '_');
    std::string Classname = y[1];
    std::string Colltype = y[2];
    //G4cout << "Classname: " << Classname << "   Colltype: " << Colltype << G4endl;
    if (Classname == "DRCalorimeter") {
      G4int NbHits = hc->GetSize();
      for (G4int ii = 0; ii < NbHits; ii++) {
        G4VHit* hit = hc->GetHit(ii);
        DRCalorimeterHit* DRHit = dynamic_cast<DRCalorimeterHit*> (hit);
        hitsVector.push_back(DRHit);
      }
      hcmap->insert(std::make_pair(hcname, hitsVector));
    } else if (Classname == "DRTSCalorimeter" && Colltype == "HC") {
      G4int NbHits = hc->GetSize();
      for (G4int ii = 0; ii < NbHits; ii++) {
        G4VHit* hit = hc->GetHit(ii);
        DRTSCalorimeterHit* DRTSHit = dynamic_cast<DRTSCalorimeterHit*> (hit);
        hitsVector.push_back(DRTSHit);
      }
      hcmap->insert(std::make_pair(hcname, hitsVector));
    } else if (Classname == "DRTSCalorimeter" && Colltype == "HC2") {
      G4int NbHits = hc->GetSize();
      for (G4int ii = 0; ii < NbHits; ii++) {
        G4VHit* hit = hc->GetHit(ii);
        DRTSCalorimeterHit2* DRTSHit = dynamic_cast<DRTSCalorimeterHit2*> (hit);
        hitsVector.push_back(DRTSHit);
      }
      hcmap->insert(std::make_pair(hcname, hitsVector));
    } else if (Classname == "Calorimeter") {
      G4int NbHits = hc->GetSize();
      for (G4int ii = 0; ii < NbHits; ii++) {
        G4VHit* hit = hc->GetHit(ii);
        CalorimeterHit* Hit = dynamic_cast<CalorimeterHit*> (hit);
        hitsVector.push_back(Hit);
      }
      hcmap->insert(std::make_pair(hcname, hitsVector));
    } else if (Classname == "Tracker") {
      G4int NbHits = hc->GetSize();
      for (G4int ii = 0; ii < NbHits; ii++) {
        G4VHit* hit = hc->GetHit(ii);
        TrackerHit* Hit = dynamic_cast<TrackerHit*> (hit);
        hitsVector.push_back(Hit);
      }
      hcmap->insert(std::make_pair(hcname, hitsVector));
    } else if (Classname == "PhotonDetector") {
      G4int NbHits = hc->GetSize();
      for (G4int ii = 0; ii < NbHits; ii++) {
        G4VHit* hit = hc->GetHit(ii);
        PhotonHit* Hit = dynamic_cast<PhotonHit*> (hit);
        hitsVector.push_back(Hit);
      }
      hcmap->insert(std::make_pair(hcname, hitsVector));
    } else {
      G4cout << "SD type: " << Classname << " unknown" << G4endl;
    }
  }
  
  //   G4cout << "Size of hcmap:  " << hcmap->size() << "  " << EventContainer->GetHCMap()->size() << G4endl;
  RootIO::GetInstance()->Write(CaTSEvt);
  CaTSEvt->Reset();
#ifdef G4ANALYSIS_USE
  analysis->EndOfEvent();
#endif
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
