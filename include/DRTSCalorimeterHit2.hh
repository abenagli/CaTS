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

#ifndef DRTSCalorimeterHit2_h
#define DRTSCalorimeterHit2_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"



class DRTSCalorimeterHit2 : public G4VHit
{
private:
  
  G4String particleName;
  G4double edep; // total energy deposit 
  G4double eobsbirks; //Birks supressed energy deposition
  G4int nceren; // nr of produced cerenkov photons
  G4ThreeVector pos; // position of calorimeter cell
  G4int timeSlice; // time of first energy deposit
  
public:
  
  DRTSCalorimeterHit2();
  DRTSCalorimeterHit2(G4String pn, G4double de, G4double eobsb, G4int nc, G4ThreeVector p,
                      G4int ts);
  ~DRTSCalorimeterHit2();
  DRTSCalorimeterHit2(const DRTSCalorimeterHit2&);
  const DRTSCalorimeterHit2& operator=(const DRTSCalorimeterHit2&);
  G4int operator==(const DRTSCalorimeterHit2&) const;
  
  inline void* operator new(size_t);
  inline void operator delete(void*);
  
  void Draw();
  void Print();
  
  
  inline void SetEdep(G4double de) {
    edep = de;
  };
  
  inline void SetParticleName(G4String pn) {
    particleName = pn;
  };
  
  inline void SetEobsbirks(G4double de) {
    eobsbirks = de;
  };
  
  inline void SetNCeren(G4int nc) {
    nceren = nc;
  };
  
  inline void SetPos(G4ThreeVector xyz) {
    pos = xyz;
  };
  
  inline void SetTimeSlice(G4int ts) {
    timeSlice = ts;
  };
  
  
  inline G4String GetParticleName() {
    return particleName;
  };
  
  inline G4double GetEdep() {
    return edep;
  };
  
  inline G4double GetEobsbirks() {
    return eobsbirks;
  };
  
  inline G4int GetNCeren() {
    return nceren;
  };
  
  inline G4ThreeVector GetPos() {
    return pos;
  };
  
  inline G4int GetTimeSlice() {
    return timeSlice;
  };
  
};


typedef G4THitsCollection<DRTSCalorimeterHit2> DRTSCalorimeterHits2Collection;

extern G4Allocator<DRTSCalorimeterHit2> DRTSCalorimeterHit2Allocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* DRTSCalorimeterHit2::operator new(size_t) {
  void *aHit;
  aHit = (void *) DRTSCalorimeterHit2Allocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void DRTSCalorimeterHit2::operator delete(void *aHit) {
  DRTSCalorimeterHit2Allocator.FreeSingle((DRTSCalorimeterHit2*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
