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
  G4float ParticleEnergy;
  G4float ParticleTime;
  G4ThreeVector ParticlePosition;
  G4ParticleMomentum ParticleMomentum;
  G4String GDMLFile;
  G4String PhysicsList;
  G4bool enableoptics;
  G4bool enablescint;
  
  std::vector<G4String>* Solids;
  std::vector<G4float>* SolidsXHalfLength;
  std::vector<G4float>* SolidsYHalfLength;
  std::vector<G4float>* SolidsZHalfLength;
  std::vector<G4String>* Volumes;
  G4int numXCells;
  G4int numYCells;
  G4int numZLayers;
  
  std::vector<G4String>* particleList;
  std::vector<G4String>* particleTypeList;
  std::vector<G4String>* processList;
  
  G4float timeslicesizelow;
  G4float mintimelow;
  G4float maxtimelow;
  G4float timeslicesizemed;
  G4float mintimemed;
  G4float maxtimemed;
  G4float timeslicesizehig;
  G4float mintimehig;
  G4float maxtimehig;
  
  static RunHeader* instance;
  
public:
  
  RunHeader();
  
  virtual ~RunHeader() {};
  
  static RunHeader* getInstance() { return instance; };
  
  G4int GetRunNumber() const {
    return RunNumber;
  }
  
  G4String GetParticleSource() {
    return ParticleSource;
  }
  
  G4String GetParticleName() {
    return ParticleName;
  }
  
  G4float GetParticleEnergy() {
    return ParticleEnergy;
  }
  
  G4float GetParticleTime() {
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
  
  std::vector<G4String>* GetSolids() {
    return Solids;
  }
  
  std::vector<G4float>* GetSolidsXHalfLength() {
    return SolidsXHalfLength;
  }
  
  std::vector<G4float>* GetSolidsYHalfLength() {
    return SolidsYHalfLength;
  }
  
  std::vector<G4float>* GetSolidsZHalfLength() {
    return SolidsZHalfLength;
  }
  
  std::vector<G4String>* GetVolumes() {
    return Volumes;
  }
  
  G4int GetXCellNum() {
    return numXCells;
  }
  
  G4int GetYCellNum() {
    return numYCells;
  }
  
  G4int GetZLayerNum() {
    return numZLayers;
  }
  
  std::vector<G4String>* GetParticleList() {
    return particleList;
  }
  
  std::vector<G4String>* GetParticleTypeList() {
    return particleTypeList;
  }
  
  std::vector<G4String>* GetProcessList() {
    return processList;
  }
  
  G4float GetTimeSliceSizeLow() {
    return timeslicesizelow;
  };
  G4float GetMinTimeLow() {
    return mintimelow;
  };
  G4float GetMaxTimeLow() {
    return maxtimelow;
  };
  G4float GetTimeSliceSizeMed() {
    return timeslicesizemed;
  };
  G4float GetMinTimeMed() {
    return mintimemed;
  };
  G4float GetMaxTimeMed() {
    return maxtimemed;
  };
  G4float GetTimeSliceSizeHig() {
    return timeslicesizehig;
  };
  G4float GetMinTimeHig() {
    return mintimehig;
  };
  G4float GetMaxTimeHig() {
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
  
  void SetParticleEnergy(G4float pe) {
    ParticleEnergy = pe;
  }
  
  void SetParticleTime(G4float t) {
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
  
  void SetSolids(std::vector<G4String>* vec) {
    for(unsigned int it = 0; it < vec->size(); ++it) 
      Solids->push_back(vec->at(it));
  }
  
  void SetSolidsXHalfLength(std::vector<G4float>* vec) {
    for(unsigned int it = 0; it < vec->size(); ++it) 
      SolidsXHalfLength->push_back(vec->at(it));
  }
  
  void SetSolidsYHalfLength(std::vector<G4float>* vec) {
    for(unsigned int it = 0; it < vec->size(); ++it) 
      SolidsYHalfLength->push_back(vec->at(it));
  }
  
  void SetSolidsZHalfLength(std::vector<G4float>* vec) {
    for(unsigned int it = 0; it < vec->size(); ++it) 
      SolidsZHalfLength->push_back(vec->at(it));
  }
  
  void SetVolumes(std::vector<G4String>* vec) {
    for(unsigned int it = 0; it < vec->size(); ++it) 
      Volumes->push_back(vec->at(it));
  }
  
  void SetXCellNum(const G4int& val) {
    numXCells = val;
  }
  
  void SetYCellNum(const G4int& val) {
    numYCells = val;
  }
  
  void SetZLayerNum(const G4int& val) {
    numZLayers = val;
  }
  
  void SetParticleList(std::vector<G4String>* vec) {
    for(unsigned int it = 0; it < vec->size(); ++it) 
      particleList->push_back(vec->at(it));
    std::sort(particleList->begin(),particleList->end());
  }
  
  void SetParticleTypeList(std::vector<G4String>* vec) {
    for(unsigned int it = 0; it < vec->size(); ++it) 
      particleTypeList->push_back(vec->at(it));
    std::sort(particleTypeList->begin(),particleTypeList->end());
  }
  
  void SetProcessList(std::vector<G4String>* vec) {
    for(unsigned int it = 0; it < vec->size(); ++it) 
      processList->push_back(vec->at(it));
    std::sort(processList->begin(),processList->end());
  }
  
  void SetTimeSliceSizeLow(G4float val) {
    timeslicesizelow = val;
  }
  
  void SetMinTimeLow(G4float val) {
    mintimelow = val;
  }
  
  void SetMaxTimeLow(G4float val) {
    maxtimelow = val;
  }
  
  void SetTimeSliceSizeMed(G4float val) {
    timeslicesizemed = val;
  }
  
  void SetMinTimeMed(G4float val) {
    mintimemed = val;
  }
  
  void SetMaxTimeMed(G4float val) {
    maxtimemed = val;
  }
  
  void SetTimeSliceSizeHig(G4float val) {
    timeslicesizehig = val;
  }
  
  void SetMinTimeHig(G4float val) {
    mintimehig = val;
  }
  
  void SetMaxTimeHig(G4float val) {
    maxtimehig = val;
  }
  
  void Print();
};

#endif	/* RunHeader_HH */
