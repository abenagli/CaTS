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
 Created on May 2, 2012, 8:44 AM
-------------------------------------------------------------------------*/

#ifndef RunHeader_HH
#define	RunHeader_HH

#include <vector>

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"



class RunHeader
{
private:
  G4int RunNumber;
  G4int RunNumber2;
  G4String ParticleSource;
  G4String ParticleName;
  G4double ParticleEnergy;
  G4double ParticleTime;
  G4ThreeVector ParticlePosition;
  G4ParticleMomentum ParticleMomentum;
  G4String GDMLFile;
  G4String PhysicsList;
  G4bool enableoptics;
  G4bool enablescint;
  
  std::vector<G4String> Volumes;
  std::vector<G4String> particleList;
  std::vector<G4String> particleListCeren;
  std::vector<G4String> particleListShort;
  
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
  
  RunHeader();
  
  virtual ~RunHeader() {};
  
  G4int GetRunNumber() const {
    return RunNumber;
  }
  
  G4String GetParticleSource() {
    return ParticleSource;
  }
  
  G4String GetParticleName() {
    return ParticleName;
  }
  
  G4double GetParticleEnergy() {
    return ParticleEnergy;
  }
  
  G4double GetParticleTime() {
    return ParticleTime;
  }
  
  G4ThreeVector GetParticlePosition() {
    return ParticlePosition;
  }
  
  G4ParticleMomentum GetParticleMomentum() {
    return ParticleMomentum;
  }
  
  G4String GetGDMLFile() {
    return GDMLFile;
  }
  
  G4String GetPhysicsList() {
    return PhysicsList;
  }
  
  G4bool Getenableoptics() {
    return enableoptics;
  }
  
  G4bool Getenablescint() {
    return enablescint;
  }
  
  std::vector<G4String> GetVolumes() {
    return Volumes;
  }
  
  std::vector<G4String> GetParticleList() {
    return particleList;
  }
  
  std::vector<G4String> GetParticleListCeren() {
    return particleListCeren;
  }
  
  std::vector<G4String> GetParticleListShort() {
    return particleListShort;
  }
  
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
  
  
  void SetRunNumber(G4int r) {
    RunNumber = r;
  }
  
  void SetParticleSource(G4String ps) {
    ParticleSource = ps;
  }
  
  void SetParticleName(G4String pn) {
    ParticleName = pn;
  }

  void SetParticleEnergy(G4double pe) {
    ParticleEnergy = pe;
  }
  
  void SetParticleTime(G4double t) {
    ParticleTime = t;
  }
  
  void SetParticlePosition(G4ThreeVector pp) {
    ParticlePosition = pp;
  }
  
  void SetParticleMomentum(G4ParticleMomentum pm) {
    ParticleMomentum = pm;
  }
  
  void SetGDMLFile(G4String f) {
    GDMLFile = f;
  }
  
  void SetPhysicsList(G4String pl) {
    PhysicsList = pl;
  }
  
  void Setenableoptics(G4bool eo) {
    enableoptics = eo;
  }
  
  void Setenablescint(G4bool es) {
    enablescint = es;
  }
  
  void SetVolumes(std::vector<G4String>* vec) {
    for(unsigned int it = 0; it < vec->size(); ++it) 
      Volumes.push_back(vec->at(it));
  }
  
  void SetParticleList(std::vector<G4String>* vec) {
    for(unsigned int it = 0; it < vec->size(); ++it) 
      particleList.push_back(vec->at(it));
  }
  
  void SetParticleListCeren(std::vector<G4String>* vec) {
    for(unsigned int it = 0; it < vec->size(); ++it) 
      particleListCeren.push_back(vec->at(it));
  }
  
  void SetParticleListShort(std::vector<G4String>* vec) {
    for(unsigned int it = 0; it < vec->size(); ++it) 
      particleListShort.push_back(vec->at(it));
  }
  
  void SetTimeSliceLow(G4double val) {
    timeslicelow = val;
  }
  
  void SetMinTimeLow(G4double val) {
    mintimelow = val;
  }
  
  void SetMaxTimeLow(G4double val) {
    maxtimelow = val;
  }
  
  void SetTimeSliceMed(G4double val) {
    timeslicemed = val;
  }
  
  void SetMinTimeMed(G4double val) {
    mintimemed = val;
  }
  
  void SetMaxTimeMed(G4double val) {
    maxtimemed = val;
  }
  
  void SetTimeSliceHig(G4double val) {
    timeslicehig = val;
  }
  
  void SetMinTimeHig(G4double val) {
    mintimehig = val;
  }
  
  void SetMaxTimeHig(G4double val) {
    maxtimehig = val;
  }
  
  void Print();
};

#endif	/* RunHeader_HH */
