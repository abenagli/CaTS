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

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
bool FillChain(TChain* chain, const std::string& inputFileList);
TH1F* GetCumulative(TProfile* prof, const int& normalize, const float& Eobs);



int main(int argc, char** argv)
{
  if( argc < 3 )
  {
    G4cout << "Program requires 2 arguments: input file list, name of output file" << G4endl;
    exit(1);
  }
  
  
  
  
  
  // //-----------
  // // histograms
  
  TH1F* h_Edep;
  TH1F* h_Eobs;
  TH1F* h_NCeren;
  
  std::map<G4String,std::map<G4String,TProfile*> > Eobs_byTime;                               // map<volume,map<timeSliceName,TProfile*> >
  std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > >Eobs_byParticleAndTime; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
  std::map<G4String,std::map<G4String,TProfile*> > EobsFrac_byTime;                               // map<volume,map<timeSliceName,TProfile*> >
  std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > >EobsFrac_byParticleAndTime; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
  
  std::map<G4String,std::map<G4String,TProfile*> > Eobs_byTime_cumul;                               // map<volume,map<timeSliceName,TProfile*> >
  std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > >Eobs_byParticleAndTime_cumul; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
  std::map<G4String,std::map<G4String,TProfile*> > EobsFrac_byTime_cumul;                               // map<volume,map<timeSliceName,TProfile*> >
  std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > >EobsFrac_byParticleAndTime_cumul; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
  
  
  std::map<G4String,std::map<G4String,TProfile2D*> > Eobs_time_vs_z; // map<volume,map<timeSliceName,TProfile2D*> >
  std::map<G4String,std::map<G4String,TProfile2D*> > Eobs_time_vs_R; // map<volume,map<timeSliceName,TProfile2D*> >
  std::map<G4String,std::map<G4String,TProfile2D*> > EobsFrac_time_vs_z; // map<volume,map<timeSliceName,TProfile2D*> >
  std::map<G4String,std::map<G4String,TProfile2D*> > EobsFrac_time_vs_R; // map<volume,map<timeSliceName,TProfile2D*> >
  
  // std::map<G4String,TProfile2D*> h2_xy_firstTime;
  // std::map<G4String,TProfile2D*> h2_zy_firstTime;
  
  G4float RunTotDepEnergy = 0.0;
  G4float RunTotObsEnergy = 0.0;
  G4float RunTotNCeren = 0.0;
  
  
  // initialize ROOT
  TSystem ts;
  gSystem->Load("libCintex");
  ROOT::Cintex::Cintex::Enable();
  gSystem->Load("libClassesDict");
  //  ROOT::Cintex::Cintex::SetDebug(2);
  
  
  // define global variables
  Event* event;
  RunHeader* runHeader;
  
  std::map<G4String,TDirectory*> VolumeDirset;
  
  std::vector<G4String>* volumeList;
  std::vector<G4String>* particleList;
  
  std::vector<G4String> timeTypes;
  timeTypes.push_back("globalTime");
  timeTypes.push_back("localTime1");
  //timeTypes.push_back("localTime2");
  
  std::vector<G4String> timeSliceTypes;
  timeSliceTypes.push_back("Low");
  timeSliceTypes.push_back("Med");
  timeSliceTypes.push_back("Hig");
  
  std::map<G4String,G4float> timeSliceSizes;
  std::map<G4String,G4float> minTimes;
  std::map<G4String,G4float> maxTimes;
  std::map<G4String,G4int> nTimeSlices;
  
  bool firstFile = true;
  
  
  // read files
  TFile* outfile = new TFile(argv[2], "RECREATE");
  
  std::ifstream inputFileList(argv[1],std::ios::in);
  std::string buffer;
  if( !inputFileList.is_open() )
  {
    std::cerr << "** ERROR: Can't open '" << argv[1] << "' for input file list" << std::endl;
    return false;
  }
  while(1)
  {
    std::getline(inputFileList,buffer);
    if( !inputFileList.good() ) break;
    if( buffer.at(0) == '#' ) continue;
    
    TFile* inFile = TFile::Open(buffer.c_str(),"READ");
    if( !inFile ) continue;
    
    TTree* Tevt = (TTree*)( inFile->Get("EventTree") );
    TTree* Trh  = (TTree*)( inFile->Get("RunTree") );
    
    event = new Event();
    //TBranch* b_event = Tevt -> GetBranch("Event");
    //b_event -> SetAddress(&event);
    Tevt -> SetBranchAddress("Event",&event);
    
    runHeader = new RunHeader();
    TBranch* b_runHeader = Trh -> GetBranch("RunHeader");
    b_runHeader -> SetAddress(&runHeader);
    
    Tevt -> SetBranchStatus("HCMap",0);
    
    
    
    G4float Ein;
    G4float Edep;
    G4float Eobs;
    G4float NCeren;

    //----------------------
    // initialize histograms    
    
    if( firstFile )
    {
      Trh->GetEntry(0);
      runHeader->Print();
      
      
      Ein = runHeader->GetParticleEnergy();
      
      volumeList = runHeader -> GetVolumes();
      volumeList->push_back("AllVol");
      particleList = runHeader -> GetParticleList();
      
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
      
      
      TDirectory* globalDir = outfile->mkdir("globalHistos");
      globalDir->cd();
      
      h_Edep = new TH1F("h_Edep", "Total energy deposited", 100, 0., 1.1*Ein);
      h_Eobs = new TH1F("h_Eobs", "Total observed (Birks suppressed) energy", 100, 0., 1.1*Ein);
      h_NCeren = new TH1F("h_NCeren", "Total nr. of cerenkov photons", 100, 0., 65000*1.1*Ein);
      
      TDirectory* timeDir = outfile->mkdir("timeHistos");
      timeDir->cd();
      
      for(unsigned int volIt = 0; volIt < volumeList->size(); ++volIt)
      {
        std::string volName = volumeList->at(volIt);
        
        // create the subdirectories for all the different SD (overkill): 
        if (VolumeDirset.find(volName) == VolumeDirset.end())
        {
          VolumeDirset.insert(std::make_pair(volName,timeDir->mkdir(volName.c_str())));
          VolumeDirset[volName]->mkdir("particleHistos");
        }
        
        for(unsigned int timeTypeIt = 0; timeTypeIt < timeTypes.size(); ++timeTypeIt)
        {
          G4String timeType = timeTypes.at(timeTypeIt);
          
          for(unsigned int timeSliceTypeIt = 0; timeSliceTypeIt < timeSliceTypes.size(); ++timeSliceTypeIt)
          {
            G4String timeSliceType = timeSliceTypes.at(timeSliceTypeIt);
            G4String timeSliceName = timeType+timeSliceType;
            G4int nBins = nTimeSlices[timeSliceType]+2;
            G4float min = minTimes[timeSliceType]-timeSliceSizes[timeSliceType];
            G4float max = maxTimes[timeSliceType]+timeSliceSizes[timeSliceType];
            
            VolumeDirset[volName]->cd();
            
            Eobs_byTime[volName][timeSliceName]     = new TProfile(Form("Eobs%s_%s",    volName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
            EobsFrac_byTime[volName][timeSliceName] = new TProfile(Form("EobsFrac%s_%s",volName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
            
            Eobs_byTime_cumul[volName][timeSliceName]     = new TProfile(Form("Eobs%s_cumul_%s",    volName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
            EobsFrac_byTime_cumul[volName][timeSliceName] = new TProfile(Form("EobsFrac%s_cumul_%s",volName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
            
            Eobs_time_vs_z[volName][timeSliceName] = new TProfile2D(Form("Eobs%s_%s_vs_z", volName.c_str(),timeSliceName.c_str()),"Eobs",1000,0.,3000.,nBins,min,max);
            Eobs_time_vs_R[volName][timeSliceName] = new TProfile2D(Form("Eobs%s_%s_vs_R", volName.c_str(),timeSliceName.c_str()),"Eobs",1000,0.,3000.,nBins,min,max);
            EobsFrac_time_vs_z[volName][timeSliceName] = new TProfile2D(Form("EobsFrac%s_%s_vs_z", volName.c_str(),timeSliceName.c_str()),"Eobs",1000,0.,3000.,nBins,min,max);
            EobsFrac_time_vs_R[volName][timeSliceName] = new TProfile2D(Form("EobsFrac%s_%s_vs_R", volName.c_str(),timeSliceName.c_str()),"Eobs",1000,0.,3000.,nBins,min,max);
            
            for(unsigned int partIt = 0; partIt < particleList->size(); ++partIt)
            {
              std::string partName = particleList->at(partIt);
              
              VolumeDirset[volName]->cd("particleHistos");
              
              Eobs_byParticleAndTime[volName][timeSliceName][partName]     = new TProfile(Form("Eobs%s_%s_%s",    volName.c_str(),partName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
              EobsFrac_byParticleAndTime[volName][timeSliceName][partName] = new TProfile(Form("EobsFrac%s_%s_%s",volName.c_str(),partName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
              
              Eobs_byParticleAndTime_cumul[volName][timeSliceName][partName]     = new TProfile(Form("Eobs%s_cumul_%s_%s",    volName.c_str(),partName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
              EobsFrac_byParticleAndTime_cumul[volName][timeSliceName][partName] = new TProfile(Form("EobsFrac%s_cumul_%s_%s",volName.c_str(),partName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
            }
          }
        }
      }
      
      firstFile = false;
    }
    
    
    //-----------------
    // loop over events
    
    G4int nEvents = Tevt->GetEntries();
    //nEvents = 1;
    G4cout << "Nr. of Events:  " << nEvents << " in input file "<< buffer << G4endl;
    for(G4int entry = 0; entry < nEvents; ++entry)
    {
      std::cout << ">>> processing event " << entry << " / " << nEvents << "\r" << std::flush;
      Tevt -> GetEntry(entry);
      
      
      // std::map<G4String, std::map<G4ThreeVector, std::vector<G4VHit*> > >* hcMap = event->GetHCMap();
      
      Edep = event->GetTotDepEnergy();
      Eobs = event->GetTotObsEnergy();
      NCeren = event->GetNCeren();
      RunTotDepEnergy = RunTotDepEnergy + Edep;
      RunTotObsEnergy = RunTotObsEnergy + Eobs;
      RunTotNCeren = RunTotNCeren + NCeren;
      h_Edep->Fill(Edep);
      h_Eobs->Fill(Eobs);
      h_NCeren->Fill(NCeren);
      
      
      std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* m_Eobs_byTime = event->GetEobsByTimeMap(); // map<Detector,map<timeSliceType,map<timeSlice,energy> > >
      std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > >* m_Eobs_byParticleAndTime = event->GetEobsByParticleAndTimeMap(); // map<Detector,map<timeSliceType,map<timeSlice,map<particle,energy> > > >
      
      
      // for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
      // {
      //   std::string volName = volumeList->at(volIt);
      
      //   for(unsigned int timeTypeIt = 0; timeTypeIt < timeTypes.size(); ++timeTypeIt)
      //   {
      //     G4String timeType = timeTypes.at(timeTypeIt);
      
      //     for(unsigned int timeSliceTypeIt = 0; timeSliceTypeIt < timeSliceTypes.size(); ++timeSliceTypeIt)
      //     {
      //       G4String timeSliceType = timeSliceTypes.at(timeSliceTypeIt);
      //       G4String timeSliceName = timeType+timeSliceType;
      //       G4int nBins = nTimeSlices[timeSliceType]+2;
      //       G4float min = minTimes[timeSliceType];
      //       G4float timeSliceSize = timeSliceSizes[timeSliceType];
      
      //       G4float sum = 0.;
      //       std::map<G4String,G4float> sumPart;
      //       std::map<G4String,G4float> sumProc;
      //       for(int bin = 0; bin < nBins; ++bin)
      //       {
      //         G4float time = min + (bin-1)*timeSliceSize + 0.5*timeSliceSize;
      //         sum += (*m_Eobs_byTime)[volName][timeSliceName][bin];
      
      //         Eobs_byTime[volName][timeSliceName]      -> Fill( time,(*m_Eobs_byTime)[volName][timeSliceName][bin] );
      //         EobsFrac_byTime[volName][timeSliceName]  -> Fill( time,(*m_Eobs_byTime)[volName][timeSliceName][bin] / Eobs );
      
      //         Eobs_byTime_cumul[volName][timeSliceName]      -> Fill( time,sum );
      //         EobsFrac_byTime_cumul[volName][timeSliceName]  -> Fill( time,sum / Eobs );
      
      //         for(unsigned int partIt = 0; partIt < particleList->size(); ++partIt)
      //         {
      //           std::string partName = particleList->at(partIt);
      //           sumPart[partName] += (*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName];
      
      //           Eobs_byParticleAndTime[volName][timeSliceName][partName]      -> Fill( time,(*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName] );
      //           EobsFrac_byParticleAndTime[volName][timeSliceName][partName]  -> Fill( time,(*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName] / Eobs );
      
      //           Eobs_byParticleAndTime_cumul[volName][timeSliceName][partName]      -> Fill( time,sumPart[partName] );
      //           EobsFrac_byParticleAndTime_cumul[volName][timeSliceName][partName]  -> Fill( time,sumPart[partName] / Eobs );
      //         }
      //       }
      //     }
      //   }
      // }
      
      for(unsigned int timeTypeIt = 0; timeTypeIt < timeTypes.size(); ++timeTypeIt)
      {
        G4String timeType = timeTypes.at(timeTypeIt);
        
        for(unsigned int timeSliceTypeIt = 0; timeSliceTypeIt < timeSliceTypes.size(); ++timeSliceTypeIt)
        {
          G4String timeSliceType = timeSliceTypes.at(timeSliceTypeIt);
          G4String timeSliceName = timeType+timeSliceType;
          G4int nBins = nTimeSlices[timeSliceType]+2;
          G4float min = minTimes[timeSliceType];
          G4float timeSliceSize = timeSliceSizes[timeSliceType];
          
          G4float sum = 0.;
          std::map<G4String,G4float> sumPart;
          std::map<G4String,G4float> sumProc;
          for(int bin = 0; bin < nBins; ++bin)
          {
            G4float time = min + (bin-1)*timeSliceSize + 0.5*timeSliceSize;
            G4float volSum = 0.;
            for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
            {
              G4String volName = volumeList->at(volIt);
              sum += (*m_Eobs_byTime)[volName][timeSliceName][bin];
              volSum += (*m_Eobs_byTime)[volName][timeSliceName][bin];
              
              Eobs_byTime[volName][timeSliceName]      -> Fill( time,(*m_Eobs_byTime)[volName][timeSliceName][bin] );
              EobsFrac_byTime[volName][timeSliceName]  -> Fill( time,(*m_Eobs_byTime)[volName][timeSliceName][bin] / Eobs );
              
              Eobs_byTime_cumul[volName][timeSliceName]      -> Fill( time,sum );
              EobsFrac_byTime_cumul[volName][timeSliceName]  -> Fill( time,sum / Eobs );
            }
            
            Eobs_byTime["AllVol"][timeSliceName]     -> Fill( time,volSum );
            EobsFrac_byTime["AllVol"][timeSliceName] -> Fill( time,volSum / Eobs );
            
            Eobs_byTime_cumul["AllVol"][timeSliceName]     -> Fill( time,sum );
            EobsFrac_byTime_cumul["AllVol"][timeSliceName] -> Fill( time,sum / Eobs );
            
            for(unsigned int partIt = 0; partIt < particleList->size(); ++partIt)
            {
              std::string partName = particleList->at(partIt);
              volSum = 0.;
              for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
              {
                G4String volName = volumeList->at(volIt);
                sumPart[partName] += (*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName];
                volSum += (*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName];
                
                Eobs_byParticleAndTime[volName][timeSliceName][partName]      -> Fill( time,(*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName] );
                EobsFrac_byParticleAndTime[volName][timeSliceName][partName]  -> Fill( time,(*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName] / Eobs );
                
                Eobs_byParticleAndTime_cumul[volName][timeSliceName][partName]      -> Fill( time,sumPart[partName] );
                EobsFrac_byParticleAndTime_cumul[volName][timeSliceName][partName]  -> Fill( time,sumPart[partName] / Eobs );
              }
              
              Eobs_byParticleAndTime["AllVol"][timeSliceName][partName]     -> Fill( time,volSum );
              EobsFrac_byParticleAndTime["AllVol"][timeSliceName][partName] -> Fill( time,volSum / Eobs );
              
              Eobs_byParticleAndTime_cumul["AllVol"][timeSliceName][partName]     -> Fill( time,sumPart[partName] );
              EobsFrac_byParticleAndTime_cumul["AllVol"][timeSliceName][partName] -> Fill( time,sumPart[partName] / Eobs );
            }
            
          }
        }
      }
    } // end loop over events
    std::cout << "\n>>> loop over events done" << std::endl;
    
    delete runHeader;
    std::cout << "qui1" << std::endl;
    delete event;
    std::cout << "qui2" << std::endl;
    inFile -> Close();
    std::cout << "qui3" << std::endl;
  }
  std::cout << "AAAA" << std::endl;
  G4cout << G4endl;
  
  
  
  // for(unsigned int volIt = 0; volIt < volumeList->size(); ++volIt)
  // {
  //   std::string volName = volumeList->at(volIt);
    
  //   for(unsigned int timeTypeIt = 0; timeTypeIt < timeTypes.size(); ++timeTypeIt)
  //   {
  //     G4String timeType = timeTypes.at(timeTypeIt);
      
  //     for(unsigned int timeSliceTypeIt = 0; timeSliceTypeIt < timeSliceTypes.size(); ++timeSliceTypeIt)
  //     {
  //       G4String timeSliceType = timeSliceTypes.at(timeSliceTypeIt);
  //       G4String timeSliceName = timeType+timeSliceType;
        
  //       VolumeDirset[volName] -> cd();
        
  //       TProfile* prof = (TProfile*)( gDirectory->Get(Form("Eobs%s_%s",volName.c_str(),timeSliceName.c_str())) );
  //       GetCumulative(prof,1,event->GetTotObsEnergy());
        
  //       for(unsigned int partIt = 0; partIt < particleList->size(); ++partIt)
  //       {
  //         G4String partName = particleList->at(partIt);
          
  //         VolumeDirset[volName] -> cd("particleHistos");
      
  //         prof = (TProfile*)( gDirectory->Get(Form("Eobs%s_%s_%s",volName.c_str(),partName.c_str(),timeSliceName.c_str())) );
  //         GetCumulative(prof,1,event->GetTotObsEnergy());
  //       }
  //     }
  //   }
  // }
  
  
  G4cout << "===========================================" << G4endl;
  G4cout << "nr of MB written:  " << outfile->Write()/1024/1024 << G4endl;
  G4cout << "===========================================" << G4endl;
  
}



std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
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



TH1F* GetCumulative(TProfile* prof, const int& normalize, const float& Eobs)
{
  int nBinsX  = prof -> GetNbinsX();
  double xMin = prof -> GetBinLowEdge(1);
  double xMax = prof -> GetBinLowEdge(nBinsX) + prof -> GetBinWidth(nBinsX);
  
  TH1F* histo_cumul;
  if     ( normalize == 0 ) histo_cumul = new TH1F(Form("cumul%s",prof->GetName()),     "",nBinsX,xMin,xMax);
  else if( normalize == 1 ) histo_cumul = new TH1F(Form("normCumul%s",prof->GetName()), "",nBinsX,xMin,xMax);
  else if( normalize == 2 ) histo_cumul = new TH1F(Form("normECumul%s",prof->GetName()),"",nBinsX,xMin,xMax);
  
  if( prof->Integral() == 0 ) return histo_cumul;

  prof -> ComputeIntegral();
  double* integral = prof -> GetIntegral();
  
  double norm = -1.;
  if( normalize == 0 ) norm = prof -> Integral();
  if( normalize == 1 ) norm = 1.;
  if( normalize == 2 ) norm = prof->Integral()/Eobs;
  
  if( normalize == 2 )
  {
    std::cout << prof->GetName() << "   Eobs: " << Eobs << "   int: " << prof->Integral() << std::endl;
  }
  
  for(int bin = 1; bin <= nBinsX; ++bin)
  {
    histo_cumul -> SetBinContent(bin,integral[bin]*norm);
  }
  
  return histo_cumul;
}
