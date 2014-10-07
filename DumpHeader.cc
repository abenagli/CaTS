// Include files
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TBranch.h"
#include "TBranchElement.h"
//
#include "Cintex/Cintex.h"
//
#include "Event.hh"
#include "RunHeader.hh"

using namespace std;


int main(int argc, char** argv) {
    // initialize ROOT
    TSystem ts;
    gSystem->Load("libCintex");
    ROOT::Cintex::Cintex::Enable();
    gSystem->Load("libClassesDict");
    //  ROOT::Cintex::Cintex::SetDebug(2);

    if (argc < 2) {
        G4cout << "Program requires 1 arguments: name of input file" << G4endl;
        exit(1);
    }
   
    TFile fo(argv[1]);
    
    RunHeader *runh = new RunHeader();
    Event *event = new Event();
    TTree *Tevt = (TTree*) fo.Get("Events");
    Tevt->SetBranchAddress("event.", &event);
    TBranch* fevtbranch = Tevt->GetBranch("event.");
    TTree *Trh = (TTree*) fo.Get("Runheader");
    Trh->SetBranchAddress("RunHeader.", &runh);
    TBranch* frunhbranch = Trh->GetBranch("RunHeader.");
    frunhbranch->GetEntry(0);
    runh->Print();
  
}

