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
#include <sstream>

#include "RootIO.hh"
#include "RootIOMessenger.hh"
#include "PrimaryGeneratorAction.hh"
//
#include "Cintex/Cintex.h"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"



static RootIO* instance = 0;

RootIO::RootIO()
{
  // initialize ROOT
  fevttree=0;
  fruntree=0;
  //TSystem ts;
  gSystem->Load("libClassesDict");
  
  ROOT::Cintex::Cintex::Enable();
  // ROOT::Cintex::Cintex::SetDebug(2);
  gDebug = 0;
  //gDebug = 1;
  FileName = "hits.root";
  evtTreeInitialized = false;
  runTreeInitialized = false;
  evtInitialized = false;
  runInitialized = false;
  pMessenger = new RootIOMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



RootIO::~RootIO()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



RootIO* RootIO::GetInstance()
{
  if( instance == 0 )
  {
    instance = new RootIO();
  }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void RootIO::BeginOfRun()
{
  G4cout << "Output File name  " << FileName << G4endl;
  fo = new TFile(FileName.c_str(), "RECREATE");
  TTree::SetMaxTreeSize(1000 * Long64_t(2000000000));
  
  // Create a ROOT Tree and one superbranch
  fevttree = new TTree("EventTree", "ROOT tree containing events");
  G4cout << "fevttree: " << fevttree << G4endl;
  fevttree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
  if( !evtTreeInitialized )
  {
    TTree::SetBranchStyle(1);
    evtTreeInitialized = true;
  }
  
  fruntree = new TTree("RunTree", "ROOT tree containing run header");
  G4cout << "fruntree: " << fruntree << G4endl;
  fruntree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
  if( !runTreeInitialized )
  {
    TTree::SetBranchStyle(1);
    runTreeInitialized = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void RootIO::Write(Event* fevent)
{
  if( !evtInitialized )
  {
    Int_t bufsize = 64000;
    fevtbranch = fevttree->Branch("Event", &fevent, bufsize, 2);
    fevtbranch->SetAutoDelete(kFALSE);
    
    for(std::map<G4String,G4int>::const_iterator mapIt = branchStatuses.begin(); mapIt != branchStatuses.end(); ++mapIt)
      fevttree->SetBranchStatus(mapIt->first,mapIt->second);
    
    evtInitialized = true;
  }
  
  fo = fevttree->GetCurrentFile(); //just in case we switched to a new file
  fevttree->Fill();
  fo->Write("", TObject::kOverwrite);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void RootIO::Write(RunHeader* fRunHeader)
{
  if( !runInitialized )
  {
    Int_t bufsize = 64000;
    frunbranch = fruntree->Branch("RunHeader", &fRunHeader, bufsize, 2);
    frunbranch->SetAutoDelete(kFALSE);
    runInitialized = true;
  }
  
  fo = fruntree->GetCurrentFile(); //just in case we switched to a new file
  fruntree->Fill();
  fo->Write("", TObject::kOverwrite);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void RootIO::Close()
{
  fo->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void RootIO::SetFileName(G4String fname)
{
  FileName = fname;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void RootIO::SetBranchStatus(G4String bname, G4int bstatus)
{
  branchStatuses[bname] = bstatus;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
