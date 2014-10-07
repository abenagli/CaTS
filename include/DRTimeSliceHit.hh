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

#ifndef DRTimeSliceHit_h
#define DRTimeSliceHit_h 1
#include "G4UnitsTable.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DRTimeSliceHit 
{
      private:
      G4int         slice;    // time slice number
      G4double      edep;     // total energy deposit 
      G4double      edepem;   // energy deposit by e+/-, gamma 
      G4double      edepnonem;// energy deposit by non e+/-, gamma 
      G4int         nceren;   // nr of produced cerenkov photons

  public:

      DRTimeSliceHit();
      DRTimeSliceHit(G4int s,G4double e, G4double eem, G4double enonem, G4int nc);
     ~DRTimeSliceHit();
      DRTimeSliceHit(const DRTimeSliceHit&);
      const DRTimeSliceHit& operator=(const DRTimeSliceHit&);
      G4int operator==(const DRTimeSliceHit&) const;

 //     inline void* operator new(size_t);
//      inline void  operator delete(void*);

      void Print();

      void SetEdep(G4int sl)           {slice = sl;};   
      void SetEdep(G4double de)        {edep = de;};
      void SetEdepEM(G4double de)      {edepem = de;};
      void SetEdepnonEM(G4double de)   {edepnonem = de;};
      void SetNCeren(G4int nc)         {nceren = nc;};
      
      G4int    GetSlice()     {return slice;}; 
      G4double GetEdep()      {return edep;}; 
      G4double GetEdepEM()    {return edepem;}; 
      G4double GetEdepnonEM() {return edepnonem;}; 
      G4int    GetNCeren()    {return nceren;};       

};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
