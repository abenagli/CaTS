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

bool FillChain(TChain* chain, const std::string& inputFileList);



int main(int argc, char** argv)
{
  if( argc < 3 )
  {
    G4cout << "Program requires 2 arguments: input file list, base name of output file" << G4endl;
    exit(1);
  }
  
  
  //----------------
  // initialize ROOT
  
  TSystem ts;
  gSystem->Load("libCintex");
  ROOT::Cintex::Cintex::Enable();
  gSystem->Load("libClassesDict");
  //  ROOT::Cintex::Cintex::SetDebug(2);
  
  
  
  //------------------------
  // define global variables
  std::vector<std::pair<int,int> > inputFrequencies;
  
  inputFrequencies.push_back(std::pair<int,int>(1,0));
  
  inputFrequencies.push_back(std::pair<int,int>(  2,  1));
  inputFrequencies.push_back(std::pair<int,int>(  4,  2));
  inputFrequencies.push_back(std::pair<int,int>(  6,  3));
  inputFrequencies.push_back(std::pair<int,int>(  8,  4));
  inputFrequencies.push_back(std::pair<int,int>( 10,  5));
  inputFrequencies.push_back(std::pair<int,int>( 12,  6));
  inputFrequencies.push_back(std::pair<int,int>( 14,  7));
  inputFrequencies.push_back(std::pair<int,int>( 16,  8));
  inputFrequencies.push_back(std::pair<int,int>( 18,  9));
  inputFrequencies.push_back(std::pair<int,int>( 20, 10));
  
  inputFrequencies.push_back(std::pair<int,int>(  1,  1));
  inputFrequencies.push_back(std::pair<int,int>(  2,  2));
  inputFrequencies.push_back(std::pair<int,int>(  3,  3));
  inputFrequencies.push_back(std::pair<int,int>(  4,  4));
  inputFrequencies.push_back(std::pair<int,int>(  5,  5));
  inputFrequencies.push_back(std::pair<int,int>(  6,  6));
  inputFrequencies.push_back(std::pair<int,int>(  7,  7));
  inputFrequencies.push_back(std::pair<int,int>(  8,  8));
  inputFrequencies.push_back(std::pair<int,int>(  9,  9));
  inputFrequencies.push_back(std::pair<int,int>( 10, 10));
  
  inputFrequencies.push_back(std::pair<int,int>(  1,  2));
  inputFrequencies.push_back(std::pair<int,int>(  2,  4));
  inputFrequencies.push_back(std::pair<int,int>(  3,  6));
  inputFrequencies.push_back(std::pair<int,int>(  4,  8));
  inputFrequencies.push_back(std::pair<int,int>(  5, 10));
  inputFrequencies.push_back(std::pair<int,int>(  6, 12));
  inputFrequencies.push_back(std::pair<int,int>(  7, 14));
  inputFrequencies.push_back(std::pair<int,int>(  8, 16));
  inputFrequencies.push_back(std::pair<int,int>(  9, 18));
  inputFrequencies.push_back(std::pair<int,int>( 10, 20));
  
  inputFrequencies.push_back(std::pair<int,int>(  1,  3));
  inputFrequencies.push_back(std::pair<int,int>(  2,  6));
  inputFrequencies.push_back(std::pair<int,int>(  3,  9));
  inputFrequencies.push_back(std::pair<int,int>(  4, 12));
  inputFrequencies.push_back(std::pair<int,int>(  5, 15));
  inputFrequencies.push_back(std::pair<int,int>(  6, 18));
  inputFrequencies.push_back(std::pair<int,int>(  7, 21));
  inputFrequencies.push_back(std::pair<int,int>(  8, 24));
  inputFrequencies.push_back(std::pair<int,int>(  9, 27));
  inputFrequencies.push_back(std::pair<int,int>( 10, 30));
 
  inputFrequencies.push_back(std::pair<int,int>(  1,  4));
  inputFrequencies.push_back(std::pair<int,int>(  2,  8));
  inputFrequencies.push_back(std::pair<int,int>(  3, 12));
  inputFrequencies.push_back(std::pair<int,int>(  4, 16));
  inputFrequencies.push_back(std::pair<int,int>(  5, 20));
  inputFrequencies.push_back(std::pair<int,int>(  6, 24));
  inputFrequencies.push_back(std::pair<int,int>(  7, 28));
  inputFrequencies.push_back(std::pair<int,int>(  8, 32));
  inputFrequencies.push_back(std::pair<int,int>(  9, 36));
  inputFrequencies.push_back(std::pair<int,int>( 10, 40));
  
  inputFrequencies.push_back(std::pair<int,int>(  1,  9));
  inputFrequencies.push_back(std::pair<int,int>(  2, 18));
  inputFrequencies.push_back(std::pair<int,int>(  3, 27));
  inputFrequencies.push_back(std::pair<int,int>(  4, 36));
  inputFrequencies.push_back(std::pair<int,int>(  5, 45));
  inputFrequencies.push_back(std::pair<int,int>(  6, 54));
  inputFrequencies.push_back(std::pair<int,int>(  7, 63));
  inputFrequencies.push_back(std::pair<int,int>(  8, 72));
  inputFrequencies.push_back(std::pair<int,int>(  9, 81));
  inputFrequencies.push_back(std::pair<int,int>( 10, 90));
  
  inputFrequencies.push_back(std::pair<int,int>(  1, 19));
  inputFrequencies.push_back(std::pair<int,int>(  2, 38));
  inputFrequencies.push_back(std::pair<int,int>(  3, 57));
  inputFrequencies.push_back(std::pair<int,int>(  4, 76));
  inputFrequencies.push_back(std::pair<int,int>(  5, 95));
  inputFrequencies.push_back(std::pair<int,int>(  6,114));
  inputFrequencies.push_back(std::pair<int,int>(  7,133));
  inputFrequencies.push_back(std::pair<int,int>(  8,152));
  inputFrequencies.push_back(std::pair<int,int>(  9,171));
  inputFrequencies.push_back(std::pair<int,int>( 10,190));
  
  inputFrequencies.push_back(std::pair<int,int>(  1, 99));
  inputFrequencies.push_back(std::pair<int,int>(  2,198));
  inputFrequencies.push_back(std::pair<int,int>(  3,297));
  inputFrequencies.push_back(std::pair<int,int>(  4,396));
  inputFrequencies.push_back(std::pair<int,int>(  5,495));
  inputFrequencies.push_back(std::pair<int,int>(  6,594));
  inputFrequencies.push_back(std::pair<int,int>(  7,693));
  inputFrequencies.push_back(std::pair<int,int>(  8,792));
  inputFrequencies.push_back(std::pair<int,int>(  9,891));
  inputFrequencies.push_back(std::pair<int,int>( 10,990));
  
  
  
  std::vector<std::map<int,int> > samplingFrequencies;
  // for(unsigned int vecIt = 0; vecIt < inputFrequencies.size(); ++vecIt)
  // {
  //   std::map<int,int> samplingFrequency;
  //   std::string dummy = inputFrequencies.at(vecIt);
  //   for(std::string::size_type it = 0; it < dummy.size(); ++it)
  //   {
  //     if( dummy[it] == '1' ) samplingFrequency[it] = 1;
  //     if( dummy[it] == '0' ) samplingFrequency[it] = 0;
  //     //std::cout << "samplingFrequency[" << it << "] = " << samplingFrequency[it] << std::endl;
  //   }
    
  //   samplingFrequencies.push_back(samplingFrequency);
  // }
  for(unsigned int vecIt = 0; vecIt < inputFrequencies.size(); ++vecIt)
  {
    std::map<int,int> samplingFrequency;
    std::pair<int,int> dummy = inputFrequencies.at(vecIt);
    for(int it = 0; it < dummy.first; ++it)
    {
      samplingFrequency[it] = 1;
      //std::cout << "samplingFrequency[" << it << "] = " << samplingFrequency[it] << std::endl;
    }
    for(int it = dummy.first; it < dummy.first+dummy.second; ++it)
    {
      samplingFrequency[it] = 0;
      //std::cout << "samplingFrequency[" << it << "] = " << samplingFrequency[it] << std::endl;
    }
    
    samplingFrequencies.push_back(samplingFrequency);
  }
  
  
  //----------------------
  // loop over input files
  
  TFile* outFile;
  std::string outFileName = Form("%s.root",argv[2]);
  outFile = new TFile(outFileName.c_str(), "RECREATE");
  
  TChain* Tevt = new TChain("EventTree","EventTree");
  TChain* Trh  = new TChain("RunTree","RunTree");
  FillChain(Tevt,std::string(argv[1]));
  FillChain(Trh,std::string(argv[1]));

  Event* event = new Event();
  Tevt -> SetBranchAddress("Event",&event);
  // Event* event = new Event();
  // TBranch* b_event = Tevt -> GetBranch("Event");
  // b_event -> SetAddress(&event);
  
  RunHeader* runHeader = new RunHeader();
  TBranch* b_runHeader = Trh -> GetBranch("RunHeader");
  b_runHeader -> SetAddress(&runHeader);
  
  
  // global variables
  std::map<G4String,TDirectory*> VolumeDirsetGlobal;
  std::map<G4String,TDirectory*> VolumeDirsetTime;
  
  std::vector<G4String>* volumeList;
  
  
  // define histograms
  std::map<G4String,TH1F*> h_Edep;
  std::map<G4String,TH1F*> h_Eobs;
  std::map<G4String,TH1F*> h_NCeren;
  
  std::map<G4String,std::map<G4String,TH1F*> > h_Eobs_sampling;
  
  std::map<G4String,std::map<G4String,TH1F*> > Eobs_byPos;
  
  G4float RunTotDepEnergy = 0.0;
  G4float RunTotObsEnergy = 0.0;
  G4float RunTotNCeren = 0.0;
  
  G4float Ein;
  G4float Edep;
  G4float Eobs;
  G4float NCeren;
  
  
  // initialize histograms    
  Trh->GetEntry(0);
  runHeader->Print();
  
  Ein = runHeader->GetParticleEnergy();
  
  volumeList = runHeader -> GetVolumes();
  
  TDirectory* globalDir = outFile->mkdir("globalHistos");
  globalDir->cd();
  
  for(unsigned int volIt = 0; volIt < volumeList->size(); ++volIt)
  {
    std::string volName = volumeList->at(volIt);
    
    // create the subdirectories for all the different SD (overkill): 
    if (VolumeDirsetGlobal.find(volName) == VolumeDirsetGlobal.end())
    {
      VolumeDirsetGlobal.insert(std::make_pair(volName,globalDir->mkdir(volName.c_str())));
    }
    
    VolumeDirsetGlobal[volName]->cd();
    
    h_Edep[volName] = new TH1F(Form("h_Edep%s",volName.c_str()), "Total energy deposited", 10000, 0., 1.1*Ein);
    h_Eobs[volName] = new TH1F(Form("h_Eobs%s",volName.c_str()), "Total observed (Birks suppressed) energy", 10000, 0., 1.1*Ein);
    h_NCeren[volName] = new TH1F(Form("h_NCeren%s",volName.c_str()), "Total nr. of cerenkov photons", 10000, 0., 65000*1.1*Ein);
    
    for(unsigned int vecIt = 0; vecIt < inputFrequencies.size(); ++vecIt)
    {
      std::string inputFrequency = Form("%dmmAct_%dmmPas",inputFrequencies.at(vecIt).first,inputFrequencies.at(vecIt).second);
      
      h_Eobs_sampling[volName][inputFrequency] = new TH1F(Form("h_Eobs%s_%s",volName.c_str(),inputFrequency.c_str()), "Observed (Birks suppressed) energy", 10000, 0., 1.1*Ein);
    }
  }
  
  
  //-----------------
  // loop over events
  
  G4int nEvents = Tevt->GetEntries();
  //nEvents = 1;
  G4cout << "Nr. of Events:  " << nEvents << std::endl;
  for(G4int entry = 0; entry < nEvents; ++entry)
  {
    G4cout << ">>> processing event " << entry << " / " << nEvents << G4endl;
    Tevt -> GetEntry(entry);
    
    Edep = event->GetTotDepEnergy();
    Eobs = event->GetTotObsEnergy();
    NCeren = event->GetNCeren();
    RunTotDepEnergy = RunTotDepEnergy + Edep;
    RunTotObsEnergy = RunTotObsEnergy + Eobs;
    RunTotNCeren = RunTotNCeren + NCeren;
    
    std::map<G4String,G4float>* m_Edep   = event->GetEdepMap();   // map<Detector,energy>
    std::map<G4String,G4float>* m_Eobs   = event->GetEobsMap();   // map<Detector,energy>
    std::map<G4String,G4float>* m_NCeren = event->GetNCerenMap(); // map<Detector,energy>
    //std::map<G4String,std::map<G4ThreeVector,G4float> >* m_EdepByPos   = event->GetEdepByPosMap();   // map<Detector,energy>
    std::map<G4String,std::map<G4ThreeVector,G4float> >* m_EobsByPos   = event->GetEobsByPosMap();   // map<Detector,energy>
    //std::map<G4String,std::map<G4ThreeVector,G4float> >* m_NCerenByPos = event->GetNCerenByPosMap(); // map<Detector,energy>
    
    
    // ------------
    // global plots
    
    for(unsigned int volIt = 0; volIt < volumeList->size(); ++volIt)
    {
      G4String volName = volumeList->at(volIt);
      
      h_Edep[volName]   -> Fill((*m_Edep)[volName]);
      h_Eobs[volName]   -> Fill((*m_Eobs)[volName]);
      h_NCeren[volName] -> Fill((*m_NCeren)[volName]);
    }
    
    
    // -------------------
    // pos dependent plots
    
    for(unsigned int volIt = 0; volIt < volumeList->size(); ++volIt)
    {      
      G4String volName = volumeList->at(volIt);
      
      std::vector<G4float> sumEobs_sampling(samplingFrequencies.size());
      
      std::map<G4ThreeVector,G4float> tmpMap = (*m_EobsByPos)[volName];
      for(std::map<G4ThreeVector,G4float>::const_iterator mapIt = tmpMap.begin(); mapIt != tmpMap.end(); ++mapIt)
      {
        for(unsigned int vecIt = 0; vecIt < samplingFrequencies.size(); ++vecIt)
        {
          std::map<int,int> samplingFrequency = samplingFrequencies.at(vecIt);
          int nSamplings = inputFrequencies.at(vecIt).first + inputFrequencies.at(vecIt).second;
          
          int scint = samplingFrequency[int((mapIt->first).z())%nSamplings];
          if( scint == 1 ) sumEobs_sampling.at(vecIt) += mapIt -> second;
          //std::cout << "samplingFrequencyIt: " << vecIt << "   pos: " << mapIt->first << "   energy: " << mapIt->second << "   sumE: " << sumEobs_sampling.at(vecIt) << "   scint? " << scint << std::endl;
        }
      }
      
      
      for(unsigned int vecIt = 0; vecIt < inputFrequencies.size(); ++vecIt)
      {
        std::string inputFrequency = Form("%dmmAct_%dmmPas",inputFrequencies.at(vecIt).first,inputFrequencies.at(vecIt).second);
        
        h_Eobs_sampling[volName][inputFrequency] -> Fill(sumEobs_sampling.at(vecIt)/6.); 
      }
    }
    
  } // end loop over events
  std::cout << "\n>>> loop over events done" << std::endl;
  
  G4cout << "===========================================" << G4endl;
  G4cout << "nr of  B written:  " << int(outFile -> Write())             << G4endl;
  G4cout << "nr of KB written:  " << int(outFile -> Write()/1024.)       << G4endl;
  G4cout << "nr of MB written:  " << int(outFile -> Write()/1024./1024.) << G4endl;
  G4cout << "===========================================" << G4endl;
  
  outFile -> Close();
  
  delete runHeader;
  delete event;
  
  G4cout << G4endl;
}



bool FillChain(TChain* chain, const std::string& inputFileList)
{
  std::ifstream inFile(inputFileList.c_str());
  std::string buffer;
  std::string bufferbeg;
  
  if( !inFile.is_open() )
  {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return false;
  }
  
  while(1)
  {
    inFile >> buffer;
    bufferbeg = buffer.substr(0,1);
    if( !inFile.good() ) break;
    if( bufferbeg == "#" ) continue;
    chain -> Add( buffer.c_str() );
    G4cout << ">>> FillChain - treeName = " << chain->GetName() << " from file " << buffer << std::endl;
  }
  
  return true;
}
