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
    G4cout << "Program requires at least 2 arguments: input file list, base name of output file, maxEntries (-1)" << G4endl;
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
  
  G4int fileIt = 0;
  TFile* outFile;
  G4int maxEntries = -1;
  if( argc > 3 ) maxEntries = G4int(atoi(argv[3]));
  
  while(1)
  {
    std::getline(inputFileList,buffer);
    if( !inputFileList.good() ) break;
    if( buffer.at(0) == '#' ) continue;
    
    std::string outFileName = Form("%s_%d.root",argv[2],fileIt);
    outFile = new TFile(outFileName.c_str(), "RECREATE");
    
    TFile* inFile = TFile::Open(buffer.c_str(),"READ");
    if( !inFile ) continue;
    else std::cout << ">>> file " << buffer << " opened" << std::endl;
    
    TTree* Tevt = (TTree*)( inFile->Get("EventTree") );
    TTree* Trh  = (TTree*)( inFile->Get("RunTree") );
    
    Event* event = new Event();
    TBranch* b_event = Tevt -> GetBranch("Event");
    b_event -> SetAddress(&event);
    
    RunHeader* runHeader = new RunHeader();
    TBranch* b_runHeader = Trh -> GetBranch("RunHeader");
    b_runHeader -> SetAddress(&runHeader);
    
    Trh->GetEntry(0);
    //runHeader->Print();
    
    
    // global variables
    G4float RunTotDepEnergy = 0.0;
    G4float RunTotObsEnergy = 0.0;
    G4float RunTotNCeren = 0.0;
    
    G4float Ein;
    G4float Edep;
    G4float Eobs;
    G4float NCeren;
    
    std::map<G4String,TDirectory*> VolumeDirsetGlobal;
    std::map<G4String,TDirectory*> VolumeDirsetTime;
    
    std::vector<G4String>* volumeList;
    
    std::vector<G4String> timeTypes;
    timeTypes.push_back("globalTime");
    timeTypes.push_back("localTime1");
    //timeTypes.push_back("localTime2");
    
    std::map<G4String,G4float> timeSliceSizes;
    std::map<G4String,G4float> minTimes;
    std::map<G4String,G4float> maxTimes;
    std::map<G4String,G4int> nTimeSlices;
    
    std::vector<G4float> EobsSumVals;
    EobsSumVals.push_back(0.000010); //  10 KeV
    EobsSumVals.push_back(0.000100); // 100 KeV
    EobsSumVals.push_back(0.000300); // 300 KeV
    EobsSumVals.push_back(0.000500); // 500 KeV
    EobsSumVals.push_back(0.000700); // 700 KeV
    EobsSumVals.push_back(0.000750); // 700 KeV
    EobsSumVals.push_back(0.000800); // 700 KeV
    EobsSumVals.push_back(0.000850); // 700 KeV
    EobsSumVals.push_back(0.000900); // 900 KeV
    EobsSumVals.push_back(0.001000); //   1 MeV
    EobsSumVals.push_back(0.003000); //   3 MeV
    EobsSumVals.push_back(0.005000); //   5 MeV
    EobsSumVals.push_back(0.010000); //  10 MeV
    EobsSumVals.push_back(0.100000); // 100 MeV
    
    
    // compute histogram axis
    std::vector<G4String>* Solids = runHeader->GetSolids();
    std::vector<G4float>* SolidsZHalfLength = runHeader->GetSolidsZHalfLength();
    std::reverse(Solids->begin(),Solids->end());
    std::reverse(SolidsZHalfLength->begin(),SolidsZHalfLength->end());
    
    G4int numZLayers = runHeader->GetZLayerNum();
    G4int nBinsZ = numZLayers*Solids->size();
    double* cellZLength = new double[Solids->size()];
    double layerZLength = 0.;
    
    for(unsigned int solidIt = 0; solidIt < Solids->size(); ++solidIt)
    {
      std::cout << "Solid: " << Solids->at(solidIt) << std::endl;
      cellZLength[solidIt] = 2. * SolidsZHalfLength->at(solidIt);
      layerZLength += cellZLength[solidIt];
      std::cout << "Solid: " << Solids->at(solidIt) << "   z: " << cellZLength[solidIt] << std::endl;
    }
    
    G4double* zAxis = new G4double[nBinsZ+1];
    for(int binz = 0; binz < numZLayers; ++binz)
    {
      for(unsigned int solidIt = 0; solidIt < Solids->size(); ++solidIt)
      {
        zAxis[binz*Solids->size()+solidIt] = binz*layerZLength + solidIt*cellZLength[solidIt>0?solidIt-1:0];
      }
    }
    zAxis[nBinsZ] = numZLayers*layerZLength;
    
    
    // define histograms
    std::map<G4String,TH1F*> h_Edep;
    std::map<G4String,TH1F*> h_Eobs;
    std::map<G4String,TH1F*> h_NCeren;
    
    std::map<G4String,std::map<G4String,std::map<G4float,TProfile*> > > p_EobsSumTime; // map<volume,map<timeType,map<EobsSum,TProfile*> > >
    std::map<G4String,std::map<G4String,std::map<G4float,TTree*> > > t_EobsSumTime; // map<volume,map<timeType,map<EobsSum,TProfile*> > >
    std::map<G4float,std::map<G4float,G4float> > t_time;
    
    
    // initialize histograms    
    
    Ein = runHeader->GetParticleEnergy();
    
    volumeList = runHeader -> GetVolumes();
    volumeList->push_back("AllVol");
    
    timeSliceSizes["Low"] = runHeader->GetTimeSliceSizeLow();
    timeSliceSizes["Med"] = runHeader->GetTimeSliceSizeMed();
    timeSliceSizes["Hig"] = runHeader->GetTimeSliceSizeHig();
    minTimes["Low"] = runHeader->GetMinTimeLow();
    minTimes["Med"] = runHeader->GetMinTimeMed();
    minTimes["Hig"] = runHeader->GetMinTimeHig();
    maxTimes["Low"] = runHeader->GetMaxTimeLow();
    maxTimes["Med"] = runHeader->GetMaxTimeMed();
    maxTimes["Hig"] = runHeader->GetMaxTimeHig();
    nTimeSlices["Low"] = int((maxTimes["Low"]-minTimes["Low"])/timeSliceSizes["Low"]);
    nTimeSlices["Med"] = int((maxTimes["Med"]-minTimes["Med"])/timeSliceSizes["Med"]);
    nTimeSlices["Hig"] = int((maxTimes["Hig"]-minTimes["Hig"])/timeSliceSizes["Hig"]);
    
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
      
      h_Edep[volName]   = new TH1F(Form("h_Edep%s",volName.c_str()),   "Total energy deposited",                   1000, 0., std::max(1.1*Ein,0.05));
      h_Eobs[volName]   = new TH1F(Form("h_Eobs%s",volName.c_str()),   "Total observed (Birks suppressed) energy", 1000, 0., std::max(1.1*Ein,0.05));
      h_NCeren[volName] = new TH1F(Form("h_NCeren%s",volName.c_str()), "Total nr. of cerenkov photons",            1000, 0., 65000*1.1*Ein);
    }
    
    TDirectory* timeDir = outFile->mkdir("timeHistos");
    timeDir->cd();
    
    for(unsigned int volIt = 0; volIt < volumeList->size(); ++volIt)
    {
      std::string volName = volumeList->at(volIt);
      
      // create the subdirectories for all the different SD (overkill): 
      if (VolumeDirsetTime.find(volName) == VolumeDirsetTime.end())
      {
        VolumeDirsetTime.insert(std::make_pair(volName,timeDir->mkdir(volName.c_str())));
      }
      
      VolumeDirsetTime[volName]->cd();
      
      for(unsigned int timeTypeIt = 0; timeTypeIt < timeTypes.size(); ++timeTypeIt)
      {
        G4String timeType = timeTypes.at(timeTypeIt);
        G4String timeSliceName = timeType + "Low";
        
        for(unsigned int EobsSumIt = 0; EobsSumIt < EobsSumVals.size(); ++EobsSumIt)
        {
          G4float EobsSum = EobsSumVals.at(EobsSumIt);
          G4String name = Form("EobsSumTime%s_%s_%.6fGeV",volName.c_str(),timeSliceName.c_str(),EobsSum);
          
          p_EobsSumTime[volName][timeSliceName][EobsSum] = new TProfile(Form("p_%s",name.c_str()),"",nBinsZ,zAxis);
          p_EobsSumTime[volName][timeSliceName][EobsSum] -> BuildOptions(0,0,"s");
          
          t_EobsSumTime[volName][timeSliceName][EobsSum] = new TTree(Form("ntu_%s",name.c_str()),Form("ntu_%s",name.c_str()));
          zAxis[nBinsZ] = numZLayers*layerZLength;
          for(int bin = 0; bin < nBinsZ-1; ++bin)
          {
            float binCenter = 0.5* ( zAxis[bin] + zAxis[bin+1] );
            t_EobsSumTime[volName][timeSliceName][EobsSum] -> Branch(Form("t_z%f",binCenter),&t_time[EobsSum][binCenter],Form("t_zBin%f/F",binCenter));
          }
        }
      }
    }
    
    
    //-----------------
    // loop over events
    
    G4int nEvents = Tevt->GetEntries();
    if( maxEntries != -1 ) nEvents = maxEntries;
    //nEvents = 1;
    G4cout << "Nr. of Events:  " << nEvents << std::endl;
    for(G4int entry = 0; entry < nEvents; ++entry)
    {
      std::cout << ">>> processing event " << entry << " / " << nEvents << "\r" << std::flush;
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
      std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4ThreeVector,G4float> > > >* m_Eobs_byPosAndTime = event->GetEobsByPosAndTimeMap(); // map<Detector,map<timeSliceType,map<timeSlice,map<position,energy> > > >
      std::map<G4String,std::map<G4String,std::map<G4ThreeVector,std::vector<G4VHit*> > > >* HCMap = event->GetHCMap();      
      
      // ------------
      // global plots
      
      G4float sumEdep = 0.;
      G4float sumEobs = 0.;
      G4float sumNCeren = 0.;
      for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
      {
        G4String volName = volumeList->at(volIt);
        
        h_Edep[volName]   -> Fill((*m_Edep)[volName]);
        h_Eobs[volName]   -> Fill((*m_Eobs)[volName]);
        h_NCeren[volName] -> Fill((*m_NCeren)[volName]);
        
        sumEdep   += (*m_Edep)[volName];
        sumEobs   += (*m_Eobs)[volName];
        sumNCeren += (*m_NCeren)[volName];
      }
      h_Edep["AllVol"]   -> Fill(sumEdep);
      h_Eobs["AllVol"]   -> Fill(sumEobs);
      h_NCeren["AllVol"] -> Fill(sumNCeren);
      
      
      //---------------------
      // time dependent plots      
      
      for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
      {
        G4String volName = volumeList->at(volIt);
        
        for(unsigned int timeTypeIt = 0; timeTypeIt < timeTypes.size(); ++timeTypeIt)
        {
          std::map<G4ThreeVector,G4float> EobsSum;
          std::map<G4float,std::map<G4float,G4float> > EobsSumTime;
          
          G4String timeType = timeTypes.at(timeTypeIt);
          G4String timeSliceName = timeType + "Low";
          
          std::map<G4int,std::map<G4ThreeVector,G4float> > tempMap = (*m_Eobs_byPosAndTime)[volName][timeSliceName];
          std::map<G4int,std::map<G4ThreeVector,G4float> >::const_iterator tempMapIt;
          
          for(tempMapIt = tempMap.begin(); tempMapIt != tempMap.end(); ++tempMapIt)
          {
            G4int nBins = nTimeSlices["Low"]+2;
            G4float min = minTimes["Low"];
            G4float timeSliceSize = timeSliceSizes["Low"];
            G4float time = min + 0.5*timeSliceSize + tempMapIt->first*timeSliceSize;
            
            std::map<G4ThreeVector,G4float> tempMap2 = tempMapIt->second;
            std::map<G4ThreeVector,G4float>::const_iterator tempMapIt2;
            
            for(tempMapIt2 = tempMap2.begin(); tempMapIt2 != tempMap2.end(); ++tempMapIt2)
            {
              G4ThreeVector pos = tempMapIt2 -> first;
              Eobs = tempMapIt2 -> second;
              
              EobsSum[pos] += Eobs;
              
              for(unsigned int EobsSumIt = 0; EobsSumIt < EobsSumVals.size(); ++EobsSumIt)
              {
                G4float EobsSumTh = EobsSumVals.at(EobsSumIt);
                
                if( EobsSum[pos] > EobsSumTh && EobsSumTime[EobsSumTh][pos.z()] == 0 )
                {
                  EobsSumTime[EobsSumTh][pos.z()] = time;
                  t_time[EobsSumTh][pos.z()] = time;
                }
              }
            }
          }
          
          for(unsigned int EobsSumIt = 0; EobsSumIt < EobsSumVals.size(); ++EobsSumIt)
          {
            G4float EobsSumTh = EobsSumVals.at(EobsSumIt);
            
            std::map<G4float,G4float> tmpMap = EobsSumTime[EobsSumTh];
            std::map<G4float,G4float>::const_iterator tmpMapIt;
            
            for(tmpMapIt = tmpMap.begin(); tmpMapIt != tmpMap.end(); ++tmpMapIt)
              p_EobsSumTime[volName][timeSliceName][EobsSumTh] -> Fill(tmpMapIt->first,tmpMapIt->second);
            
            t_EobsSumTime[volName][timeSliceName][EobsSumTh] -> Fill();
          }
        }
      }
      
      
      // //---------------------
      // // time dependent plots      
      
      // for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
      // {
      //   G4String volName = volumeList->at(volIt);
        
      //   for(unsigned int timeTypeIt = 0; timeTypeIt < timeTypes.size(); ++timeTypeIt)
      //   {
      //     G4String timeType = timeTypes.at(timeTypeIt);
      //     G4String timeSliceName = timeType + "Low";
      //     G4float min = minTimes["Low"];
      //     G4float timeSliceSize = timeSliceSizes["Low"];
          
      //     std::map<G4float,G4float> EobsSum;
      //     std::map<G4float,std::map<G4float,G4float> > EobsSumTime;
          
      //     std::map<G4ThreeVector,std::vector<G4VHit*> > posHitMap = (*HCMap)[volName][timeSliceName];
      //     for(std::map<G4ThreeVector,std::vector<G4VHit*> >::const_iterator mapIt = posHitMap.begin(); mapIt != posHitMap.end(); ++mapIt)
      //     {
      //       double z = (mapIt->first).z();
      //       std::vector<G4VHit*> hitVec = mapIt->second;
      //       for(unsigned int vecIt = 0; vecIt < hitVec.size(); ++vecIt)
      //       {
      //         DRTSCalorimeterHit2* aHit = dynamic_cast<DRTSCalorimeterHit2*>(hitVec.at(vecIt));
      //         G4int timeSlice = aHit -> GetTimeSlice();
      //         G4float time = min + (timeSlice-1)*timeSliceSize + 0.5*timeSliceSize;
              
      //         EobsSum[z] += aHit->GetEobsbirks();
              
      //         for(unsigned int EobsSumIt = 0; EobsSumIt < EobsSumVals.size(); ++EobsSumIt)
      //         {
      //           G4float EobsSumTh = EobsSumVals.at(EobsSumIt);
                
      //           if( EobsSum[z] > EobsSumTh && EobsSumTime[EobsSumTh][z] == 0 )
      //           {
      //             EobsSumTime[EobsSumTh][z] = time;
      //           }
      //         }
      //       }
      //     }
          
      //     for(unsigned int EobsSumIt = 0; EobsSumIt < EobsSumVals.size(); ++EobsSumIt)
      //     {
      //       G4float EobsSumTh = EobsSumVals.at(EobsSumIt);
            
      //       std::map<G4float,G4float> tmpMap = EobsSumTime[EobsSumTh];
      //       std::map<G4float,G4float>::const_iterator tmpMapIt;
            
      //       for(tmpMapIt = tmpMap.begin(); tmpMapIt != tmpMap.end(); ++tmpMapIt)
      //         p_EobsSumTime[volName][timeSliceName][EobsSumTh] -> Fill(tmpMapIt->first,tmpMapIt->second);
      //     }
      //   }
      // }
      
    } // end loop over events
    std::cout << "\n>>> loop over events done" << std::endl;
    
    int bytes = outFile -> Write();
    G4cout << "============================================"  << G4endl;
    G4cout << "nr of  B written:  " << int(bytes)             << G4endl;
    G4cout << "nr of KB written:  " << int(bytes/1024.)       << G4endl;
    G4cout << "nr of MB written:  " << int(bytes/1024./1024.) << G4endl;
    G4cout << "============================================" << G4endl;
    
    outFile -> Close();
    ++fileIt;
    
    delete runHeader;
    delete event;
    inFile -> Close();
  }
  
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
