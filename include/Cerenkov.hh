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

#ifndef Cerenkov_h
#define Cerenkov_h 1


#include <vector>
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4Material.hh" 
#include "G4MaterialPropertyVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Cerenkov {
public:
    Cerenkov();
    ~Cerenkov();
    void Initialize();
    G4double GetAverageNumberOfPhotons(const G4double, const G4double, const G4String) const;
private:
    std::vector<G4String> CAI;
    std::vector<G4PhysicsOrderedFreeVector*> CerenkovAngleIntegrals;
    std::vector<G4MaterialPropertyVector*> RefractionIndeces;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

