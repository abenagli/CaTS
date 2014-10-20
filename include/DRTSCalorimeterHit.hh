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

#ifndef DRTSCalorimeterHit_h
#define DRTSCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "DRTimeSliceHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DRTSCalorimeterHit : public G4VHit {
private:

    G4double edep; // total energy deposit 
    G4double edepem; // energy deposit by e+/-, gamma 
    G4double edepnonem; // energy deposit by non e+/-, gamma
    G4double eobsbirks; //Birks supressed energy deposition
    G4int nceren; // nr of produced cerenkov photons
    G4ThreeVector pos; // position of calorimeter cell
    G4double globalTime; // time of first energy deposit
    G4double localTime1; // time of first energy deposit
    G4double localTime2; // time of first energy deposit
    G4double dedx; // energy loss per distance traveled
    G4double vel;// momentum of
    std::vector<DRTimeSliceHit> fDRGlobalTimeSliceHitCollection;
    std::vector<DRTimeSliceHit> fDRLocalTime1SliceHitCollection;
    std::vector<DRTimeSliceHit> fDRLocalTime2SliceHitCollection;
 
   
public:

    DRTSCalorimeterHit();
    //
    // The vector of DRTimeSliceHit is created automatically when the constructor is called 
    // and therefore we don't need to have it as an input argument
    //
    DRTSCalorimeterHit(G4double e, G4double eem, G4double enonem, G4double eobsb, G4int nc, G4ThreeVector p,
                       G4double gt, G4double lt1, G4double lt2, G4double ts,
                       G4double el, G4double v);
    ~DRTSCalorimeterHit();
    DRTSCalorimeterHit(const DRTSCalorimeterHit&);
    const DRTSCalorimeterHit& operator=(const DRTSCalorimeterHit&);
    G4int operator==(const DRTSCalorimeterHit&) const;

    inline void* operator new(size_t);
    inline void operator delete(void*);

    void Draw();
    void Print();

    void AddDRTimeSliceHit(G4int,G4int,G4double,G4double,G4double,G4int);
    void AddDRTimeSliceHit(DRTimeSliceHit,G4int);

    inline std::vector<DRTimeSliceHit>* GetDRGlobalTimeSliceHitCol() {
        return &fDRGlobalTimeSliceHitCollection;
    };

    inline std::vector<DRTimeSliceHit>* GetDRLocalTime1SliceHitCol() {
        return &fDRLocalTime1SliceHitCollection;
    };

    inline std::vector<DRTimeSliceHit>* GetDRLocalTime2SliceHitCol() {
        return &fDRLocalTime2SliceHitCollection;
    };

    inline G4int GetNrofTimeSliceHits(G4int type) {
      if( type == 0 ) return fDRGlobalTimeSliceHitCollection.size();
      if( type == 1 ) return fDRLocalTime1SliceHitCollection.size();
      if( type == 2 ) return fDRLocalTime2SliceHitCollection.size();
      return -1;
    };

    inline void SetEdep(G4double de) {
        edep = de;
    };

    inline void SetEdepEM(G4double de) {
        edepem = de;
    };

    inline void SetEdepnonEM(G4double de) {
        edepnonem = de;
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

    inline void SetGlobalTime(G4double de) {
        globalTime = de;
    };

    inline void SetLocalTime1(G4double de) {
        localTime1 = de;
    };

    inline void SetLocalTime2(G4double de) {
        localTime2 = de;
    };

    inline void SetDedx(G4double de) {
        dedx = de;
    };
    inline void SetVel(G4double v){
        vel = v;
    }

    inline G4double GetEdep() {
        return edep;
    };

    inline G4double GetEdepEM() {
        return edepem;
    };

    inline G4double GetEdepnonEM() {
        return edepnonem;
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

    inline G4double GetGlobalTime() {
        return globalTime;
    };

    inline G4double GetLocalTime1() {
        return localTime1;
    };

    inline G4double GetLocalTime2() {
        return localTime2;
    };
  
    inline G4double GetDedx() {
        return dedx;
    };
    inline G4double GetV(){
        return vel;
    };
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<DRTSCalorimeterHit> DRTSCalorimeterHitsCollection;

extern G4Allocator<DRTSCalorimeterHit> DRTSCalorimeterHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* DRTSCalorimeterHit::operator new(size_t) {
    void *aHit;
    aHit = (void *) DRTSCalorimeterHitAllocator.MallocSingle();
    return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void DRTSCalorimeterHit::operator delete(void *aHit) {
    DRTSCalorimeterHitAllocator.FreeSingle((DRTSCalorimeterHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
