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
#include "RunAction.hh"
#include "RunActionMessenger.hh"
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
#include "G4Timer.hh"

RunAction* RunAction::instance = 0;

PrimaryGeneratorAction* RunAction::pgA;



RunAction::RunAction(G4String gf, G4String pl, G4bool eo, G4bool es)
{
  pMessenger = new RunActionMessenger(this);
  particleList = new std::vector<G4String>;
  processList = new std::vector<G4String>;
  gdmlFile = gf;
  PhysicsList = pl;
  enableoptics = eo;
  enablescint = es;
  timer = new G4Timer;
  
  instance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



RunAction::~RunAction()
{
  delete timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
#ifdef G4ANALYSIS_USE
  Analysis* analysis = Analysis::getInstance();
  analysis->BeginOfRun(aRun->GetRunID());
#endif
  
  RunHeader* RH = new RunHeader();
  
  RootIO::GetInstance()->BeginOfRun();
  StackingAction::getInstance()->BeginRun();
  timer->Start();
  pgA = PrimaryGeneratorAction::getInstance();
  
  G4ParticleGun* pvG = (G4ParticleGun*) pgA->GetGenerator();
  RH->SetRunNumber(aRun->GetRunID());
  RH->SetParticleSource(pgA->GetGeneratorName());
  RH->SetParticleName(pvG->GetParticleDefinition()->GetParticleName());
  RH->SetParticleEnergy(pvG->GetParticleEnergy()/CLHEP::GeV);
  RH->SetParticleTime(pvG->GetParticleTime()/CLHEP::ns);
  RH->SetParticlePosition(pvG->GetParticlePosition()/CLHEP::mm);
  RH->SetParticleMomentum(pvG->GetParticleMomentumDirection());
  RH->SetGDMLFile(gdmlFile);
  RH->SetPhysicsList(PhysicsList);
  RH->Setenableoptics(enableoptics);
  RH->Setenablescint(enablescint);
  
  DetectorConstruction* DC = DetectorConstruction::getInstance();
  RH->SetSolids(DC->GetSolids());
  RH->SetSolidsXHalfLength(DC->GetSolidsXHalfLength());
  RH->SetSolidsYHalfLength(DC->GetSolidsYHalfLength());
  RH->SetSolidsZHalfLength(DC->GetSolidsZHalfLength());
  RH->SetVolumes(DC->GetVolumes());
  RH->SetXCellNum(DC->GetXCellNum());
  RH->SetYCellNum(DC->GetYCellNum());
  RH->SetZLayerNum(DC->GetZLayerNum());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void RunAction::EndOfRunAction(const G4Run* aRun)
{
  RunHeader* RH = RunHeader::getInstance();
  
  RH->SetParticleList(particleList);
  RH->SetProcessList(processList);
  RH->SetTimeSliceSizeLow(timeSliceSizes["Low"]);
  RH->SetMinTimeLow(minTimes["Low"]);
  RH->SetMaxTimeLow(maxTimes["Low"]);
  RH->SetTimeSliceSizeMed(timeSliceSizes["Med"]);
  RH->SetMinTimeMed(minTimes["Med"]);
  RH->SetMaxTimeMed(maxTimes["Med"]);
  RH->SetTimeSliceSizeHig(timeSliceSizes["Hig"]);
  RH->SetMinTimeHig(minTimes["Hig"]);
  RH->SetMaxTimeHig(maxTimes["Hig"]);
  
  RootIO::GetInstance()->Write(RH);
  
  timer->Stop();
  G4cout << "number of events = " << aRun->GetNumberOfEvent()
         << " " << *timer
         << G4endl;
  
#ifdef G4ANALYSIS_USE
  Analysis* analysis = Analysis::getInstance();
  analysis->EndOfRun(aRun->GetRunID());
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
