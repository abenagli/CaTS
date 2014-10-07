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

#include "G4Material.hh" 
#ifdef G4ANALYSIS_USE
#include "Analysis.hh"
#include "EventAction.hh"
#include "TH1.h"
#endif

class G4Step;
class G4HCofThisEvent;
class Cerenkov;
class DRTSCalorimeterSDMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DRTSCalorimeterSD : public G4VSensitiveDetector {

private:
  G4double birksc1;
  G4double birksc2;
  G4double timeslicelow;
  G4double mintimelow;
  G4double maxtimelow;
  G4double timeslicemed;
  G4double mintimemed;
  G4double maxtimemed;
  G4double timeslicehig;
  G4double mintimehig;
  G4double maxtimehig;

public:
  DRTSCalorimeterSD(G4String);
  ~DRTSCalorimeterSD();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  
  void SetBirksc1(G4double c1) {
    birksc1 = c1;
  };
  void SetBirksc2(G4double c2) {
    birksc2 = c2;
  };
  void SetTimeSliceLow(G4double tslice) {
    timeslicelow = tslice;
  };
  void SetMinTimeLow(G4double time) {
    mintimelow = time;
  };
  void SetMaxTimeLow(G4double time) {
    maxtimelow = time;
  };
  void SetTimeSliceMed(G4double tslice) {
    timeslicemed = tslice;
  };
  void SetMinTimeMed(G4double time) {
    mintimemed = time;
  };
  void SetMaxTimeMed(G4double time) {
    maxtimemed = time;
  };
  void SetTimeSliceHig(G4double tslice) {
    timeslicehig = tslice;
  };
  void SetMinTimeHig(G4double time) {
    mintimehig = time;
  };
  void SetMaxTimeHig(G4double time) {
    maxtimehig = time;
  };
  
  G4double GetBirksc1() {
    return birksc1;
  };
  G4double GetBirksc2() {
    return birksc2;
  };
  G4double GetTimeSliceLow() {
    return timeslicelow;
  };
  G4double GetMinTimeLow() {
    return mintimelow;
  };
  G4double GetMaxTimeLow() {
    return maxtimelow;
  };
  G4double GetTimeSliceMed() {
    return timeslicemed;
  };
  G4double GetMinTimeMed() {
    return mintimemed;
  };
  G4double GetMaxTimeMed() {
    return maxtimemed;
  };
  G4double GetTimeSliceHig() {
    return timeslicehig;
  };
  G4double GetMinTimeHig() {
    return mintimehig;
  };
  G4double GetMaxTimeHig() {
    return maxtimehig;
  };
  
  std::vector<G4String>* GetParticleList() { return particleList; };
  std::vector<G4String>* GetParticleListCeren() { return particleListCeren; };
  std::vector<G4String>* GetParticleListShort() { return particleListShort; };
  
  static DRTSCalorimeterSD* getInstance() { return instance; };
  
private:
  DRTSCalorimeterHitsCollection* drtscalorimeterCollection;
  DRTSCalorimeterHits2Collection* drtscalorimeterCollection2;
  G4int HCID;
  static EventAction* EvtAction;
  static RunHeader* RunHead;
  Cerenkov* CerenGenerator;
  DRTSCalorimeterSDMessenger* pMessenger;
  std::vector<G4String>* particleList;
  std::vector<G4String>* particleListCeren;
  std::vector<G4String>* particleListShort;
  
  static DRTSCalorimeterSD* instance;
  
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
