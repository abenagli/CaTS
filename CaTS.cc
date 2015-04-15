/* ------------------------------------------------------------------------
            |\___/|       
            )     (    
           =\     /=
             )===(
            /     \         Cats: Calorimeter and Tracker Simulation
            |     |         Author: Hans Wenzel (Fermilab)
           /       \
           \       /
            \__  _/
              ( (
               ) )
              (_(
-------------------------------------------------------------------------*/
//-----------------
// Geant4 includes:
//-----------------
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4PhysListFactory.hh"
#include "G4OpticalPhysics.hh"
#include "G4VModularPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4String.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif
#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif
#include "G4GDMLParser.hh"
//--------------
// CaTS include:
//--------------
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "StackingAction.hh"
#include "TrackingAction.hh"
#include "RootIO.hh"
#ifdef G4ANALYSIS_USE
#include "Analysis.hh"
#endif



int main(int argc, char **argv)
{
  if (argc < 2)
  {
    G4cout << "Error! Mandatory input file is not specified!" << G4endl;
    G4cout << G4endl;
    G4cout << G4endl;
    G4cout << "Usage: CATS <intput_gdml_file:mandatory>"
           << G4endl;
    G4cout << G4endl;
    return -1;
  }
  
  // #ifdef G4ANALYSIS_USE
  //   Analysis* pAnalysis = Analysis::getInstance();
  // #endif
  
    
  //Construct the default run manager
  G4RunManager * runManager = new G4RunManager();
  
  
  // random seeds
  runManager -> SetRandomNumberStore(true);
  
  
  // PhysicsList
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = NULL;
  G4String physName = "";
  //-----------------------------------------------------
  // Physics List name defined via environmental variable
  // The following 19 physics lists are available:
  //  CHIPS
  //  FTFP_BERT
  //  FTFP_BERT_TRV
  //  FTFP_BERT_HP
  //  FTF_BIC 
  //  LBE
  //  LHEP
  //  QBBC
  //  QGSC_BERT
  //  QGSP
  //  QGSP_BERT
  //  QGSP_BERT_CHIPS
  //  QGSP_BERT_HP
  //  QGSP_BIC
  //  QGSP_BIC_HP
  //  QGSP_FTFP_BERT
  //  QGS_BIC
  //  QGSP_INCLXX
  //  Shielding
  //-----------------------------------------------------
  char* path = getenv("PHYSLIST");
  if( path ) physName = G4String(path);
  else       physName = "FTFP_BERT";     // default
  //else       physName = "QGSP_BERT_LIV"; // livermore
  
  // reference PhysicsList via its name
  if( factory.IsReferencePhysList(physName) ) {
    phys = factory.GetReferencePhysList(physName);
  }
  
  
  // Optical Physics
  // 
  // since the cerenkov contribution is calculated in the sensitive detector (DRCAL/DRTSCal) 
  // it's not necessary to enable optical physics for a dual readout calorimeter
  // (only if you are interested in actually tracing the photons) 
  //
  G4bool enable_scint = false;
  G4bool enable_optical = false;
  char* opticalphys = getenv("ENABLEOPTICAL");
  
  if( opticalphys )
  {
    enable_optical = true;
    G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
    phys->RegisterPhysics(opticalPhysics);
    //
    //  opticalPhysics->Configure(kNoProcess, true);
    //  The following optical processes (indeces) are available
    //    
    //   kCerenkov,      ///< Cerenkov process index
    //   kScintillation, ///< Scintillation process index
    //   kAbsorption,    ///< Absorption process index
    //   kRayleigh,      ///< Rayleigh scattering process index
    //   kMieHG,         ///< Mie scattering process index
    //   kBoundary,      ///< Boundary process index
    //   kWLS,           ///< Wave Length Shifting process index
    //   kNoProcess      ///< Number of processes, no selected process
    //
    opticalPhysics->Configure(kCerenkov, true);
    char* enablescintphys = getenv("ENABLESCINTILLATION");
    if (enablescintphys) {
      enable_scint = true;
      opticalPhysics->Configure(kScintillation, true);
      G4cout << "Scintillation and Cerenkov enabled" << G4endl;
    } else {
      G4cout << "Scintillation disabled and Cerenkov enabled" << G4endl;
    }
  }
  else
  {
    G4cout << "Optical physics processes (scintillation/Cerenkov) disabled" << G4endl;
  }
  phys->DumpList();
  phys->DumpCutValuesTable();
  runManager->SetUserInitialization(phys);
  
  
  runManager->SetUserInitialization(new DetectorConstruction(argv[1]));
  runManager->SetUserAction(new RunAction(argv[1], physName, enable_optical, enable_scint));
  runManager->SetUserAction(new PrimaryGeneratorAction);
  EventAction *EvtAct = EventAction::GetInstance();
  runManager->SetUserAction(EvtAct);
  runManager->SetUserAction(new TrackingAction());
  // the stacking action is basically only useful when doing analysis
  // on the particles produced
#ifdef G4ANALYSIS_USE
  runManager->SetUserAction(StackingAction::getInstance());
#endif
  
  // #ifdef G4ANALYSIS_USE
  //   Analysis* pAnalysis = Analysis::getInstance();
  // #endif
  runManager->Initialize();
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  RootIO* pRootIO = RootIO::GetInstance();
  
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  if (argc == 3) // batch mode  
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[2];
    UImanager->ApplyCommand(command + fileName);
  } else // interactive mode
  {
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute vis.mac");
#endif
    ui->SessionStart();
    delete ui;
#endif
  }
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  
  
  return 0;
}
