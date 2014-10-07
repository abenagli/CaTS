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

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorActionMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GeneralParticleSource.hh"
#include "G4HEPEvtInterface.hh"
PrimaryGeneratorAction* PrimaryGeneratorAction::instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction* PrimaryGeneratorAction::getInstance() {
    if (instance == 0) instance = new PrimaryGeneratorAction();
    return instance;
}

PrimaryGeneratorAction::PrimaryGeneratorAction() {

    //    const char* filename = "pythia_event.data";
    //    gentypeMap["HEPEvt"] = new G4HEPEvtInterface(filename);
    G4int n_particle = 1;
    G4ParticleGun* fParticleGun = new G4ParticleGun(n_particle);
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "mu+");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    fParticleGun->SetParticleEnergy(10. * CLHEP::GeV);
    fParticleGun->SetParticlePosition(G4ThreeVector(0. * CLHEP::cm, 0. * CLHEP::cm, -130. * CLHEP::cm));
    gentypeMap["particleGun"] = fParticleGun;
    gentypeMap["GPS"] = new G4GeneralParticleSource;
    //create a messenger for this class
    gunMessenger = new PrimaryGeneratorActionMessenger(this);
    currentGenerator = gentypeMap["particleGun"];
    currentGeneratorName = "particleGun";
    instance=this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() {

    //    delete gunMessenger;
    //    for (std::map<G4String, G4VPrimaryGenerator*>::iterator ii = gentypeMap.begin(); ii != gentypeMap.end(); ++ii) {
    //        delete (*ii).second;
    //        gentypeMap.erase((*ii).first);
    //    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    currentGenerator->GeneratePrimaryVertex(anEvent);
}
