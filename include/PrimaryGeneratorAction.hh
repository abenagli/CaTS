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

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <map>

//class G4ParticleGun;
//class G4GeneralParticleSource;
//class G4HEPEvtInterface;
class G4Event;
class PrimaryGeneratorActionMessenger;
class G4VPrimaryGenerator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
private:

    G4VPrimaryGenerator* currentGenerator;
    G4String currentGeneratorName;
    std::map<G4String, G4VPrimaryGenerator*> gentypeMap;
    PrimaryGeneratorActionMessenger* gunMessenger;
    static PrimaryGeneratorAction* instance;


public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction();
    virtual void GeneratePrimaries(G4Event* anEvent);
    void SetGenerator(G4VPrimaryGenerator* gen);
    void SetGenerator(G4String genname);

    G4VPrimaryGenerator* GetGenerator() const;
    G4String GetGeneratorName() const;
    static PrimaryGeneratorAction* getInstance();
};

// ====================================================================
// inline functions
// ====================================================================

inline void PrimaryGeneratorAction::SetGenerator(G4VPrimaryGenerator* gen) {
    currentGenerator = gen;
}

inline void PrimaryGeneratorAction::SetGenerator(G4String genname) {
    std::map<G4String, G4VPrimaryGenerator*>::iterator pos = gentypeMap.find(genname);
    if (pos != gentypeMap.end()) {
        currentGenerator = pos->second;
        currentGeneratorName = genname;
    }
}

inline G4VPrimaryGenerator* PrimaryGeneratorAction::GetGenerator() const {
    return currentGenerator;
}

inline G4String PrimaryGeneratorAction::GetGeneratorName() const {
    return currentGeneratorName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*PrimaryGeneratorAction_h*/
