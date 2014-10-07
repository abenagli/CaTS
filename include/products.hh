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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....../* 
#ifndef PRODUCTS_HH
#define	PRODUCTS_HH
#include "G4Types.hh"

struct products {
public:
    G4int NParticles;
    G4double kinE;
    G4double totE;

    products() {
    }

    virtual ~products() {
    }
    
};


#endif	/* PRODUCTS_HH */

