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

//

static RootIO* instance = 0;

RootIO::RootIO() {
    // initialize ROOT
  ftree=0;
  fhtree=0;
    TSystem ts;
    gSystem->Load("libClassesDict");

    //    ROOT::Cintex::Cintex::SetDebug(2);
    ROOT::Cintex::Cintex::Enable();
    //gDebug = 1;
    gDebug = 0;
    FileName = "hits.root";
    treeinitialized = false;
    runinitialized = false;
    evtinitialized = false;
    pMessenger = new RootIOMessenger(this);
}

RootIO::~RootIO() {
}

RootIO* RootIO::GetInstance() {
    if (instance == 0) {
        instance = new RootIO();
    }
    return instance;
}

void RootIO::BeginOfRun() {

    G4cout << "Output File name  " << FileName << G4endl;
    fo = new TFile(FileName.c_str(), "RECREATE");
    TTree::SetMaxTreeSize(1000 * Long64_t(2000000000));
    // Create a ROOT Tree and one superbranch
    ftree = new TTree("Events", "ROOT tree containing Hit collections");
    G4cout << "ftree: " << ftree << G4endl;
    ftree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
    if (!treeinitialized) {
        Int_t branchStyle = 1;
        TTree::SetBranchStyle(branchStyle);
        treeinitialized = true;
    }
  fhtree = new TTree("Runheader", "ROOT tree containing Hit collections");
    G4cout << "fhtree: " << fhtree << G4endl;
    fhtree->SetAutoSave(1000000000); // autosave when 1 Gbyte written
    //if (!treeinitialized) {
      //Int_t branchStyle = 1;
	//  TTree::SetBranchStyle(branchStyle);
        //treeinitialized = true;
    //}
}

void RootIO::Write(Event* fevent) {

    if (!evtinitialized) {
        Int_t bufsize = 64000;
        fevtbranch = ftree->Branch("event.", &fevent, bufsize, 0);
        fevtbranch->SetAutoDelete(kFALSE);
        evtinitialized = true;
        /*
        std::map<G4String, std::map<G4String, products> >* pmap = fevent->Getprocessmap();
        std::map<G4String, std::map<G4String, products> >::iterator pmap_iter;
        for (pmap_iter = pmap->begin(); pmap_iter != pmap->end(); pmap_iter++) {
            G4cout << "(*pmap_iter).first: " << (*pmap_iter).first << "(*pmap_iter).second.size():    " << (*pmap_iter).second.size() << G4endl;
            std::map<G4String, products>::iterator proc_iter;
            for (proc_iter = (*pmap_iter).second.begin(); proc_iter != (*pmap_iter).second.end(); proc_iter++) {
                G4cout << "pname:  " << (*proc_iter).first
                        << "  nr.:    " << (*proc_iter).second.NParticles
                        << "  kinE.:    " << (*proc_iter).second.kinE
                        << "  totE.:    " << (*proc_iter).second.totE
                        << G4endl;
            }
        }
        */
    }


    fo = ftree->GetCurrentFile(); //just in case we switched to a new file
    // std::map<G4String, std::map<G4int, std::map<G4String, G4int> > >* pmap = fevent->GetProcessAndGlobalTimeMult();
    // std::map<G4String, std::map<G4int, std::map<G4String, G4int> > >::const_iterator mapIt;
    // for(mapIt = pmap->begin(); mapIt != pmap->end(); ++mapIt)
    // {
    //   std::map<G4int, std::map<G4String, G4int> > pmap2 = mapIt -> second;
    //   std::map<G4int, std::map<G4String, G4int> >::const_iterator mapIt2;
      
    //   for(mapIt2 = pmap2.begin(); mapIt2 != pmap2.end(); ++mapIt2)
    //   {
    //     std::map<G4String, G4int> pmap3 = mapIt2 -> second;
    //     std::map<G4String, G4int>::const_iterator mapIt3;
        
    //     for(mapIt3 = pmap3.begin(); mapIt3 != pmap3.end(); ++mapIt3)
    //     {
    //       G4cout << "processMult[" << mapIt->first << "][" << mapIt2->first << "][" << mapIt3->first << "] = " << mapIt3->second << G4endl;
    //     }
    //   }
    // }
    // G4cout << "done" << G4endl;
      
    fnb += ftree->Fill();
    fo->Write("", TObject::kOverwrite);
    //    fevts++;
    //   G4cout << "fnb: " << fnb << "   " << fevts << G4endl;
}

void RootIO::Write(RunHeader* fRunHeader) {

    if (!runinitialized) {
        Int_t bufsize = 64000;
        frunbranch = fhtree->Branch("RunHeader.", &fRunHeader, bufsize, 0);
        frunbranch->SetAutoDelete(kFALSE);
        runinitialized = true;
    }

    fo = fhtree->GetCurrentFile(); //just in case we switched to a new file
    fnb += fhtree->Fill();
    fo->Write("", TObject::kOverwrite);
    //    fevts++;
    //   G4cout << "fnb: " << fnb << "   " << fevts << G4endl;
}

void RootIO::Close() {
    fo->Close();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootIO::SetFileName(G4String fname) {
    FileName = fname;
}
