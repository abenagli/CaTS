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
#include "G4Timer.hh"
#include "RunAction.hh"
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

PrimaryGeneratorAction* RunAction::pgA;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(G4String gf, G4String pl, G4bool eo, G4bool es) {
    gdmlFile = gf;
    PhysicsList = pl;
    enableoptics = eo;
    enablescint = es;
    timer = new G4Timer;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction() {
    delete timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun) {
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
#ifdef G4ANALYSIS_USE
    Analysis* analysis = Analysis::getInstance();
    analysis->BeginOfRun(aRun->GetRunID());
#endif
    RootIO::GetInstance()->BeginOfRun();
    StackingAction::getInstance()->BeginRun();
    timer->Start();
    pgA = PrimaryGeneratorAction::getInstance();
    G4cout << "particle generator:  " << pgA->GetGeneratorName() << G4endl;
    G4ParticleGun* pvG = (G4ParticleGun*) pgA->GetGenerator();
    G4cout << "Particle energy: " << pvG->GetParticleEnergy() << G4endl;
    G4cout << "Particle name: " << pvG->GetParticleDefinition()->GetParticleName() << G4endl;
    RunHeader* rh = new RunHeader();
    rh->SetRunNumber(aRun->GetRunID());
    rh->SetParticleSource(pgA->GetGeneratorName());
    rh->SetParticleName(pvG->GetParticleDefinition()->GetParticleName());
    rh->SetParticleEnergy(pvG->GetParticleEnergy()/CLHEP::GeV);
    rh->SetParticleTime(pvG->GetParticleTime()/CLHEP::ns);
    rh->SetParticlePosition(pvG->GetParticlePosition()/CLHEP::mm);
    rh->SetParticleMomentum(pvG->GetParticleMomentumDirection());
    rh->SetGDMLFile(gdmlFile);
    rh->SetPhysicsList(PhysicsList);
    rh->Setenableoptics(enableoptics);
    rh->Setenablescint(enablescint);
    DetectorConstruction* DC = DetectorConstruction::getInstance();
    rh->SetVolumes(DC->GetVolumes());
    DRTSCalorimeterSD* cSD = DRTSCalorimeterSD::getInstance();
    rh->SetParticleList(cSD->GetParticleList());
    rh->SetTimeSliceSizeLow(cSD->GetTimeSliceSizeLow());
    rh->SetMinTimeLow(cSD->GetMinTimeLow());
    rh->SetMaxTimeLow(cSD->GetMaxTimeLow());
    rh->SetTimeSliceSizeMed(cSD->GetTimeSliceSizeMed());
    rh->SetMinTimeMed(cSD->GetMinTimeMed());
    rh->SetMaxTimeMed(cSD->GetMaxTimeMed());
    rh->SetTimeSliceSizeHig(cSD->GetTimeSliceSizeHig());
    rh->SetMinTimeHig(cSD->GetMinTimeHig());
    rh->SetMaxTimeHig(cSD->GetMaxTimeHig());
    RootIO::GetInstance()->Write(rh);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun) {
    timer->Stop();
    G4cout << "number of events = " << aRun->GetNumberOfEvent()
            << " " << *timer << G4endl;
#ifdef G4ANALYSIS_USE
    Analysis* analysis = Analysis::getInstance();
    analysis->EndOfRun(aRun->GetRunID());
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
