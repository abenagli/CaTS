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

#ifndef DRTSCalorimeterSD_h
#define DRTSCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "DRTSCalorimeterHit.hh"
#include "DRTSCalorimeterHit2.hh"
#include "RunHeader.hh"
#include "RunAction.hh"

#include "G4Material.hh" 
#ifdef G4ANALYSIS_USE
#include "Analysis.hh"
#include "EventAction.hh"
#include "TH1F.h"
#endif

class G4Step;
class G4HCofThisEvent;
class Cerenkov;
class DRTSCalorimeterSDMessenger;

G4String GetParticleName(G4Track* aTrack);



class DRTSCalorimeterSD : public G4VSensitiveDetector
{
private:
  G4float birksc1;
  G4float birksc2;
  std::vector<G4String> timeSliceTypes;
  std::map<G4String,G4float> timeSliceSizes;
  std::map<G4String,G4float> minTimes;
  std::map<G4String,G4float> maxTimes;
  std::map<G4String,G4float> times;
  
  
public:
  DRTSCalorimeterSD(G4String);
  ~DRTSCalorimeterSD();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  
  void SetBirksc1(G4float c1) {
    birksc1 = c1;
  };
  void SetBirksc2(G4float c2) {
    birksc2 = c2;
  };
  
  G4float GetBirksc1() {
    return birksc1;
  };
  G4float GetBirksc2() {
    return birksc2;
  };
  
  static DRTSCalorimeterSD* getInstance() { return instance; };
  
private:
  static DRTSCalorimeterSD* instance;
  static EventAction* EvtAction;
  static RunHeader* RunHead;
  Cerenkov* CerenGenerator;
  DRTSCalorimeterSDMessenger* pMessenger;
  
  std::vector<G4String>* particleList;
  std::vector<G4String>* processList;
  
#ifdef G4ANALYSIS_USE
    static Analysis* analysis;
    TDirectory *mydir;
    TH1F* dEdxweighted;
    TH1F* dEdxunweighted;
    TH1F* Birksweighted;
    TH1F* Birksunweighted;
    TH1F* globalPositionX;
    TH1F* globalPositionY;
    TH1F* globalPositionZ;
    TH1F* localPositionX;
    TH1F* localPositionY;
    TH1F* localPositionZ;
#endif
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
