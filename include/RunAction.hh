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
#ifndef RunAction_h
#define RunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"
#include "G4Timer.hh"

#include "RunActionMessenger.hh"
#include "G4Run.hh"
#include "G4ParticleGun.hh"
#ifdef G4ANALYSIS_USE
#include "Analysis.hh"
#endif
#include "RootIO.hh"
#include "DetectorConstruction.hh"
#include "StackingAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DRTSCalorimeterSD.hh"
#include "RunHeader.hh"

class RunActionMessenger;
class G4Run;
class PrimaryGeneratorAction;
class G4Timer;



class RunAction : public G4UserRunAction
{
public:
  RunAction(G4String fname, G4String pname, G4bool eo, G4bool es);
  ~RunAction();
  void BeginOfRunAction(const G4Run* aRun);
  void EndOfRunAction(const G4Run* aRun);

  static RunAction* getInstance() { return instance; };
  
private:
  G4Timer* timer;
  G4String gdmlFile;
  G4String PhysicsList;
  G4bool enableoptics;
  G4bool enablescint;
  RunActionMessenger* pMessenger;
  
  std::map<G4String,G4float> timeSliceSizes;
  std::map<G4String,G4float> minTimes;
  std::map<G4String,G4float> maxTimes;
  
  std::vector<G4String>* particleList;
  std::vector<G4String>* processList;
  
  static PrimaryGeneratorAction* pgA; // pointer to the particle source 
  
  static RunAction* instance;
  
public:
  void SetTimeSliceSize(const G4float& tslice, const G4String& type) {
    timeSliceSizes[type] = tslice;
  };
  void SetMinTime(const G4float& time, const G4String& type) {
    minTimes[type] = time;
  };
  void SetMaxTime(const G4float& time, const G4String& type) {
    maxTimes[type] = time;
  };
  
  G4float GetTimeSliceSize(const G4String& type) {
    return timeSliceSizes[type];
  };
  G4float GetMinTime(const G4String& type) {
    return minTimes[type];
  };
  G4float GetMaxTime(const G4String& type) {
    return maxTimes[type];
  };
  
  std::map<G4String,G4float> GetTimeSliceSizes() { return timeSliceSizes; };
  std::map<G4String,G4float> GetMinTimes()       { return       minTimes; };
  std::map<G4String,G4float> GetMaxTimes()       { return       maxTimes; };
  
  std::vector<G4String>* GetParticleList() { return particleList; };
  std::vector<G4String>* GetProcessList()  { return processList;  };
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*RunAction_h*/
