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
// Include files
// Root:
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TSystem.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH1D.h"
//
#include "Cintex/Cintex.h"
//
//Cats:
#include "Classes.hh"
#include "Event.hh"
#include "RunHeader.hh"
#include "TrackerHit.hh"
#include "CalorimeterHit.hh"
#include "DRCalorimeterHit.hh"
#include "DRTSCalorimeterHit.hh"
#include "PhotonHit.hh"
#include "TrackerHit.hh"
//
//STL:
#include <vector>
#include<iomanip>
#include <fstream>
using namespace std;



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    G4cout << "Program requires 1 argument: input file list, label" << G4endl;
    exit(1);
  }
  
  
  //----------------
  // initialize ROOT
  
  TSystem ts;
  gSystem->Load("libCintex");
  ROOT::Cintex::Cintex::Enable();
  gSystem->Load("libClassesDict");
  //  ROOT::Cintex::Cintex::SetDebug(2);
  
  
  
  //---------------------
  // read input file list
  
  std::ifstream inputFileList(argv[1],std::ios::in);
  std::string buffer;
  if( !inputFileList.is_open() )
  {
    std::cerr << "** ERROR: Can't open '" << argv[1] << "' for input file list" << std::endl;
    return false;
  }
  
  
  //----------------------
  // loop over input files
  
  int fileIt = 0;
  while(1)
  {
    std::getline(inputFileList,buffer);
    if( !inputFileList.good() ) break;
    if( buffer.at(0) == '#' ) continue;
    
    TFile* inFile = TFile::Open(buffer.c_str(),"READ");
    if( !inFile ) continue;
    else std::cout << ">>> file " << buffer << " opened" << std::endl;
    
    TTree* Trh  = (TTree*)( inFile->Get("RunTree") );
    
    RunHeader* runHeader = new RunHeader();
    TBranch* b_runHeader = Trh -> GetBranch("RunHeader");
    
    if( !Trh->GetListOfBranches()->FindObject("RunHeader") )
    {
      std::cout << ">>> RunHeader not found <<<" << std::endl;
      continue;
    }
    
    b_runHeader -> SetAddress(&runHeader);
    
    std::vector<G4String>* processList;
    
    Trh->GetEntry(0);
    //runHeader->Print();
    
    processList  = runHeader -> GetProcessList();
    
    std::ofstream outFile(Form("processList_%s_%d.txt",argv[2],fileIt),std::ios::out);
    for(unsigned int vecIt = 0; vecIt < processList->size(); ++vecIt)
    {
      outFile << processList->at(vecIt) << "\n";
    }
    outFile.close();
    
    ++fileIt;
  }
  
  G4cout << G4endl;
}
