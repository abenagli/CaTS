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



int main(int argc, char** argv)
{
  if( argc < 3 )
  {
    G4cout << "Program requires 2 arguments: input file list, name of output file" << G4endl;
    exit(1);
  }
  
  
  
  
  
  // //-----------
  // // histograms
  
  std::map<G4String,std::map<G4String,TProfile*> > Eobs_byTime;                               // map<volume,map<timeSliceName,TProfile*> >
  std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > >Eobs_byParticleAndTime; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
  
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
  
  
  TFile* outfile = new TFile(argv[2], "RECREATE");
  
  
  TChain* Tevt = new TChain("EventTree","EventTree");
  TChain* Trh = new TChain("RunTree","RunTree");
  FillChain(Tevt,argv[1]);
  FillChain(Trh,argv[1]);
  
  Event* event = new Event();
  TBranch* b_event = Tevt -> GetBranch("Event");
  b_event -> SetAddress(&event);
  
  RunHeader* runHeader = new RunHeader();
  TBranch* b_runHeader = Trh -> GetBranch("RunHeader");
  b_runHeader -> SetAddress(&runHeader);
  
  Tevt -> SetBranchStatus("HCMap",0);
  
  
  
  Trh->GetEntry(0);
  runHeader->Print();
  G4float Ein = runHeader->GetParticleEnergy();
  G4float Edep;
  G4float Eobs;
  G4float NCeren;
  std::vector<G4String> volumeList = runHeader -> GetVolumes();
  volumeList.push_back("AllVol");
  std::vector<G4String> particleList = runHeader -> GetParticleList();
  
  std::vector<G4String> timeTypes;
  timeTypes.push_back("globalTime");
  timeTypes.push_back("localTime1");
  timeTypes.push_back("localTime2");
  
  std::vector<G4String> timeSliceTypes;
  timeSliceTypes.push_back("Low");
  timeSliceTypes.push_back("Med");
  timeSliceTypes.push_back("Hig");
  
  std::map<G4String,G4float> timeSliceSizes;
  std::map<G4String,G4float> minTimes;
  std::map<G4String,G4float> maxTimes;
  std::map<G4String,G4int> nTimeSlices;
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
  
  
  
  //------------------
  // define histograms
  
  std::map<G4String, TDirectory*> VolumeDirset;
  
  TDirectory* globalDir = outfile->mkdir("globalHistos");
  globalDir->cd();
  
  TH1F* h_Edep = new TH1F("h_Edep", "Total energy deposited", 100, 0., 1.1*Ein);
  TH1F* h_Eobs = new TH1F("h_Eobs", "Total observed (Birks suppressed) energy", 100, 0., 1.1*Ein);
  TH1F* h_NCeren = new TH1F("h_NCeren", "Total nr. of cerenkov photons", 100, 0., 65000*1.1*Ein);
  
  
  
  TDirectory* timeDir = outfile->mkdir("timeHistos");
  timeDir->cd();
  
  for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
  {
    std::string volName = volumeList.at(volIt);
    
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

        Eobs_byTime[volName][timeSliceName] = new TProfile(Form("Eobs%s_%s",volName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max);
        
        for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
        {
          std::string partName = particleList.at(partIt);
          
          VolumeDirset[volName]->cd("particleHistos");
        
          Eobs_byParticleAndTime[volName][timeSliceName][partName] = new TProfile(Form("Eobs%s_%s_%s",volName.c_str(),partName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max);
        }
      }
    }
  }
  
  
  
  //-----------------
  // loop over events
  
  G4int nEvents = Tevt->GetEntries();
  Tevt -> GetListOfBranches();
  G4cout << " Nr. of Events:  " << nEvents << " in input file "<< G4endl;
  for(G4int entry = 0; entry < nEvents; ++entry)
  {
    std::cout << ">>> processing event " << entry << " / " << nEvents << "\r" << std::flush;
    
    Long64_t localEntry = Tevt->LoadTree(entry);
    b_event->GetEntry(localEntry);
    
    
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
    
    
    std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* m_Eobs_byTime = event->GetEdepByTimeMap(); // map<Detector,map<timeSliceType,map<timeSlice,energy> > >
    std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > >* m_Eobs_byParticleAndTime = event->GetEdepByParticleAndTimeMap(); // map<Detector,map<timeSliceType,map<timeSlice,map<particle,energy> > > >
    
    
    for(unsigned int volIt = 0; volIt < volumeList.size()-1; ++volIt)
    {
      std::string volName = volumeList.at(volIt);
      
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
          
          for(int bin = 0; bin <= nBins+1; ++bin)
          {
            G4float time = min + (bin-1)*timeSliceSize + 0.5*timeSliceSize;
            Eobs_byTime[volName][timeSliceName] -> Fill( time,(*m_Eobs_byTime)[volName][timeSliceName][bin] );
            
            for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
            {
              std::string partName = particleList.at(partIt);
              Eobs_byParticleAndTime[volName][timeSliceName][partName] -> Fill( time,(*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName] );
            }
          }
        }
      }
    }
    // std::map<G4String, std::map<G4String, G4float> >* E_byPart = event->GetE_byParticle();
    // std::map<G4String, std::map<G4String, G4float> >::iterator Eiter;
    // for (Eiter = E_byPart->begin(); Eiter != E_byPart->end(); ++Eiter) {
    //   std::map<G4String, G4float>::iterator partiter;
    //   for (partiter = E_byPart->at((*Eiter).first).begin(); partiter != E_byPart->at((*Eiter).first).end(); partiter++) {
    //     E_byParticle[(*Eiter).first][(*partiter).first] += (*partiter).second;
    //   }
    // }
    
    // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >* E_byPartAndGlobTimeMed = event->GetE_byParticleAndGlobalTimeMed();
    // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >::iterator ETimeMediter;
    // for (ETimeMediter = E_byPartAndGlobTimeMed->begin(); ETimeMediter != E_byPartAndGlobTimeMed->end(); ETimeMediter++) {
    //   std::map<G4int, std::map<G4String, G4float> >::iterator timeiter;        
    //   for (timeiter = ((*ETimeMediter).second).begin(); timeiter != ((*ETimeMediter).second).end(); timeiter++) {
    //     std::map<G4String, G4float>::iterator partiter;
    //     for (partiter = ((*timeiter).second).begin(); partiter != ((*timeiter).second).end(); ++partiter) {
    //       E_byParticleAndGlobalTimeMed[(*ETimeMediter).first][(*timeiter).first][(*partiter).first] += (*partiter).second;
    //     }
    //   }
    // }
    
    
    // //
    // // histogram the relative energy contribution by particle:
    // //        
    // //G4cout << ">>> histogram the relative energy contribution by particle:" << G4endl;
    
    // for (Eiter = E_byPart->begin(); Eiter != E_byPart->end(); Eiter++) {
    //   G4String volName = (*Eiter).first;
    //   VolumeDirset[volName]->cd();
    //   if (EvtEnergybySDandpart.find(volName) == EvtEnergybySDandpart.end()) {
    //     std::map<G4String, TH1F*> tmpmap;
    //     EvtEnergybySDandpart.insert(std::make_pair(volName, tmpmap));
    //     std::map<G4String, G4float>::iterator partiter;
    //     for (partiter = E_byPart->at(volName).begin(); partiter != E_byPart->at(volName).end(); partiter++) {
    //       const char* hisname = (*partiter).first.c_str();
    //       TH1F *histo = new TH1F(hisname, "energy by particle", 200, 0, 100);
    //       EvtEnergybySDandpart[volName].insert(std::make_pair((*partiter).first, histo));
    //       EvtEnergybySDandpart[volName].at((*partiter).first)->Fill(((*partiter).second * 0.1) / Ein);
    //     }
    //   } else {
    //     std::map<G4String, G4float>::iterator partiter;
    //     for (partiter = E_byPart->at(volName).begin(); partiter != E_byPart->at(volName).end(); partiter++) {
    //       G4String PName = (*partiter).first;
    //       if (EvtEnergybySDandpart[volName].find(PName) == EvtEnergybySDandpart[volName].end()) {
    //         const char* hisname = PName.c_str();
    //         TH1F *histo = new TH1F(hisname, "energy by particle", 200, 0, 100);
    //         EvtEnergybySDandpart[volName].insert(std::make_pair(PName, histo));
    //         EvtEnergybySDandpart[volName].at(PName)->Fill(((*partiter).second * 0.1) / Ein);
    //       } else {
    //         EvtEnergybySDandpart[volName].at(PName)->Fill(((*partiter).second * 0.1) / Ein);
    //       }
    //     }
    //   }
    // }
    

    // //------------------------------------------------------------------------------------------------------
    // std::map<G4String, std::map<G4String, G4float> >* Eobs_byPart = event->GetEobs_byParticle();
    // std::map<G4String, std::map<G4String, G4float> >::iterator Eobsiter;
    // for (Eobsiter = Eobs_byPart->begin(); Eobsiter != Eobs_byPart->end(); Eobsiter++) {
    //   std::map<G4String, G4float>::iterator partobsiter;
    //   for (partobsiter = Eobs_byPart->at((*Eobsiter).first).begin(); partobsiter != Eobs_byPart->at((*Eobsiter).first).end(); partobsiter++) {
    //     Eobs_byParticle[(*Eobsiter).first][(*partobsiter).first] += (*partobsiter).second;
    //   }
    // }
  
    
    // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >* Eobs_byPartAndGlobTimeLow = event->GetEobs_byParticleAndGlobalTimeLow();
    // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >* Eobs_byPartAndGlobTimeMed = event->GetEobs_byParticleAndGlobalTimeMed();
    // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >* Eobs_byPartAndGlobTimeHig = event->GetEobs_byParticleAndGlobalTimeHig();
    
    // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >* Eobs_byPartAndLocTime1Low = event->GetEobs_byParticleAndLocalTime1Low();
    // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >* Eobs_byPartAndLocTime1Med = event->GetEobs_byParticleAndLocalTime1Med();
    // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >* Eobs_byPartAndLocTime1Hig = event->GetEobs_byParticleAndLocalTime2Hig();
    
    // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >* Eobs_byPartAndLocTime2Low = event->GetEobs_byParticleAndLocalTime2Low();
    // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >* Eobs_byPartAndLocTime2Med = event->GetEobs_byParticleAndLocalTime2Med();
    // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >* Eobs_byPartAndLocTime2Hig = event->GetEobs_byParticleAndLocalTime2Hig();
    
    
    // for(int sliceLowIt = 0; sliceLowIt <= nSlicesLow+1; ++sliceLowIt)
    // {
    //   G4float time = sliceLowIt * timeSliceLow - timeSliceLow;
    //   G4float sum_globTime = 0.;
    //   G4float sum_locTime1 = 0.;
    //   G4float sum_locTime2 = 0.;
    //   for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
    //   {
    //     std::string partName = particleList.at(partIt);
    //     for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
    //     {
    //       std::string volName = volumeList.at(volIt);
    //       sum_globTime += (*Eobs_byPartAndGlobTimeLow)[volName][sliceLowIt][partName];
    //       sum_locTime1 += (*Eobs_byPartAndLocTime1Low)[volName][sliceLowIt][partName];
    //       sum_locTime2 += (*Eobs_byPartAndLocTime2Low)[volName][sliceLowIt][partName];
    //     }
    //   }
    //   EvtEnergybyglobaltimelow -> Fill(time,sum_globTime);
    //   EvtEnergybylocaltime1low -> Fill(time,sum_locTime1);
    //   EvtEnergybylocaltime2low -> Fill(time,sum_locTime2);
    // }
    // for(int sliceMedIt = 0; sliceMedIt <= nSlicesMed+1; ++sliceMedIt)
    // {
    //   G4float time = sliceMedIt * timeSliceMed - timeSliceMed;
    //   G4float sum_globTime = 0.;
    //   G4float sum_locTime1 = 0.;
    //   G4float sum_locTime2 = 0.;
    //   for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
    //   {
    //     std::string partName = particleList.at(partIt);
    //     for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
    //     {
    //       std::string volName = volumeList.at(volIt);
    //       sum_globTime += (*Eobs_byPartAndGlobTimeMed)[volName][sliceMedIt][partName];
    //       sum_locTime1 += (*Eobs_byPartAndLocTime1Med)[volName][sliceMedIt][partName];
    //       sum_locTime2 += (*Eobs_byPartAndLocTime2Med)[volName][sliceMedIt][partName];
    //     }
    //   }
    //   EvtEnergybyglobaltimemed -> Fill(time,sum_globTime);
    //   EvtEnergybylocaltime1med -> Fill(time,sum_locTime1);
    //   EvtEnergybylocaltime2med -> Fill(time,sum_locTime2);
    // }  
    // for(int sliceHigIt = 0; sliceHigIt <= nSlicesHig+1; ++sliceHigIt)
    // {
    //   G4float time = sliceHigIt * timeSliceHig - timeSliceHig;
    //   G4float sum_globTime = 0.;
    //   G4float sum_locTime1 = 0.;
    //   G4float sum_locTime2 = 0.;
    //   for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
    //   {
    //     std::string partName = particleList.at(partIt);
    //     for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
    //     {
    //       std::string volName = volumeList.at(volIt);
    //       sum_globTime += (*Eobs_byPartAndGlobTimeHig)[volName][sliceHigIt][partName];
    //       sum_locTime1 += (*Eobs_byPartAndLocTime1Hig)[volName][sliceHigIt][partName];
    //       sum_locTime2 += (*Eobs_byPartAndLocTime2Hig)[volName][sliceHigIt][partName];
    //     }
    //   }
    //   EvtEnergybyglobaltimehig -> Fill(time,sum_globTime);
    //   EvtEnergybylocaltime1hig -> Fill(time,sum_locTime1);
    //   EvtEnergybylocaltime2hig -> Fill(time,sum_locTime2);
    // }
  
    
    // for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
    // {
    //   std::string partName = particleList.at(partIt);
    //   //G4cout << "partName: " << partName << G4endl;
      
    //   for(int sliceLowIt = 0; sliceLowIt <= nSlicesLow+1; ++sliceLowIt)
    //   {
    //     G4float time = sliceLowIt * timeSliceLow - timeSliceLow;
    //     G4float sum_globTime = 0.;
    //     G4float sum_locTime1 = 0.;
    //     G4float sum_locTime2 = 0.;
    //     for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
    //     {
    //       std::string volName = volumeList.at(volIt);
    //       sum_globTime += (*Eobs_byPartAndGlobTimeLow)[volName][sliceLowIt][partName];
    //       sum_locTime1 += (*Eobs_byPartAndLocTime1Low)[volName][sliceLowIt][partName];
    //       sum_locTime2 += (*Eobs_byPartAndLocTime2Low)[volName][sliceLowIt][partName];
    //     }
    //     EvtEnergybypartandglobaltimelow[partName] -> Fill(time,sum_globTime);
    //     EvtEnergybypartandlocaltime1low[partName] -> Fill(time,sum_locTime1);
    //     EvtEnergybypartandlocaltime2low[partName] -> Fill(time,sum_locTime2);
    //   }
    //   for(int sliceMedIt = 0; sliceMedIt <= nSlicesMed+1; ++sliceMedIt)
    //   {
    //     G4float time = sliceMedIt * timeSliceMed - timeSliceMed;
    //     G4float sum_globTime = 0.;
    //     G4float sum_locTime1 = 0.;
    //     G4float sum_locTime2 = 0.;
    //     for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
    //     {
    //       std::string volName = volumeList.at(volIt);
    //       sum_globTime += (*Eobs_byPartAndGlobTimeMed)[volName][sliceMedIt][partName];
    //       sum_locTime1 += (*Eobs_byPartAndLocTime1Med)[volName][sliceMedIt][partName];
    //       sum_locTime2 += (*Eobs_byPartAndLocTime2Med)[volName][sliceMedIt][partName];
    //     }
    //     EvtEnergybypartandglobaltimemed[partName] -> Fill(time,sum_globTime);
    //     EvtEnergybypartandlocaltime1med[partName] -> Fill(time,sum_locTime1);
    //     EvtEnergybypartandlocaltime2med[partName] -> Fill(time,sum_locTime2);
    //   }
      
    //   for(int sliceHigIt = 0; sliceHigIt <= nSlicesHig+1; ++sliceHigIt)
    //   {
    //     G4float time = sliceHigIt * timeSliceHig - timeSliceHig;
    //     G4float sum_globTime = 0.;
    //     G4float sum_locTime1 = 0.;
    //     G4float sum_locTime2 = 0.;
    //     for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
    //     {
    //       std::string volName = volumeList.at(volIt);
    //       sum_globTime += (*Eobs_byPartAndGlobTimeHig)[volName][sliceHigIt][partName];
    //       sum_locTime1 += (*Eobs_byPartAndLocTime1Hig)[volName][sliceHigIt][partName];
    //       sum_locTime2 += (*Eobs_byPartAndLocTime2Hig)[volName][sliceHigIt][partName];
    //     }
    //     EvtEnergybypartandglobaltimehig[partName] -> Fill(time,sum_globTime);
    //     EvtEnergybypartandlocaltime1hig[partName] -> Fill(time,sum_locTime1);
    //     EvtEnergybypartandlocaltime2hig[partName] -> Fill(time,sum_locTime2);
    //   }
    // }
    
    
    // for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
    // {
    //   std::string volName = volumeList.at(volIt);
    //   //G4cout << "volName: " << volName << G4endl;
      
    //   for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
    //   {
    //     std::string partName = particleList.at(partIt);
    //     //G4cout << "partName: " << partName << G4endl;
        
    //     for(int sliceLowIt = 0; sliceLowIt <= nSlicesLow+1; ++sliceLowIt)
    //     {
    //       G4float time = sliceLowIt * timeSliceLow - timeSliceLow;
    //       EvtEnergybySDandpartandglobaltimelow[volName][partName] -> Fill(time,(*Eobs_byPartAndGlobTimeLow)[volName][sliceLowIt][partName]);
    //       EvtEnergybySDandpartandlocaltime1low[volName][partName] -> Fill(time,(*Eobs_byPartAndLocTime1Low)[volName][sliceLowIt][partName]);
    //       EvtEnergybySDandpartandlocaltime2low[volName][partName] -> Fill(time,(*Eobs_byPartAndLocTime2Low)[volName][sliceLowIt][partName]);
    //     }
        
    //     for(int sliceMedIt = 0; sliceMedIt <= nSlicesMed+1; ++sliceMedIt)
    //     {
    //       G4float time = sliceMedIt * timeSliceMed - timeSliceMed;
    //       EvtEnergybySDandpartandglobaltimemed[volName][partName] -> Fill(time,(*Eobs_byPartAndGlobTimeMed)[volName][sliceMedIt][partName]);
    //       EvtEnergybySDandpartandlocaltime1med[volName][partName] -> Fill(time,(*Eobs_byPartAndLocTime1Med)[volName][sliceMedIt][partName]);
    //       EvtEnergybySDandpartandlocaltime2med[volName][partName] -> Fill(time,(*Eobs_byPartAndLocTime2Med)[volName][sliceMedIt][partName]);
    //     }
        
    //     for(int sliceHigIt = 0; sliceHigIt <= nSlicesHig+1; ++sliceHigIt)
    //     {
    //       G4float time = sliceHigIt * timeSliceHig - timeSliceHig;
    //       EvtEnergybySDandpartandglobaltimehig[volName][partName] -> Fill(time,(*Eobs_byPartAndGlobTimeHig)[volName][sliceHigIt][partName]);
    //       EvtEnergybySDandpartandlocaltime1hig[volName][partName] -> Fill(time,(*Eobs_byPartAndLocTime1Hig)[volName][sliceHigIt][partName]);
    //       EvtEnergybySDandpartandlocaltime2hig[volName][partName] -> Fill(time,(*Eobs_byPartAndLocTime2Hig)[volName][sliceHigIt][partName]);
    //     }
    //   }
    // }
    
  
    // // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >::iterator EobsTimeLowiter;
    // // for (EobsTimeLowiter = Eobs_byPartAndGlobTimeLow->begin(); EobsTimeLowiter != Eobs_byPartAndGlobTimeLow->end(); EobsTimeLowiter++) {
    // //   std::map<G4int, std::map<G4String, G4float> >::iterator timeiter;
    // //   for (timeiter = ((*EobsTimeLowiter).second).begin(); timeiter != ((*EobsTimeLowiter).second).end(); timeiter++) {
    // //     std::map<G4String, G4float>::iterator partobsiter;
    // //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
    // //       Eobs_byParticleAndGlobalTimeLow[(*EobsTimeLowiter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
    // //       EvtEnergybySDandpartandglobaltimelow[(*EobsTimeLowiter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
    // //     }
    // //   }
    // // }
    // // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >::iterator EobsTimeMediter;
    // // for (EobsTimeMediter = Eobs_byPartAndGlobTimeMed->begin(); EobsTimeMediter != Eobs_byPartAndGlobTimeMed->end(); EobsTimeMediter++) {
    // //   std::map<G4int, std::map<G4String, G4float> >::iterator timeiter;
    // //   for (timeiter = ((*EobsTimeMediter).second).begin(); timeiter != ((*EobsTimeMediter).second).end(); timeiter++) {
    // //     std::map<G4String, G4float>::iterator partobsiter;
    // //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
    // //       Eobs_byParticleAndGlobalTimeMed[(*EobsTimeMediter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
    // //       EvtEnergybySDandpartandglobaltimemed[(*EobsTimeMediter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
    // //     }
    // //   }
    // // }
    // // std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >::iterator EobsTimeHigiter;
    // // for (EobsTimeHigiter = Eobs_byPartAndGlobTimeHig->begin(); EobsTimeHigiter != Eobs_byPartAndGlobTimeHig->end(); EobsTimeHigiter++) {
    // //   std::map<G4int, std::map<G4String, G4float> >::iterator timeiter;
    // //   for (timeiter = ((*EobsTimeHigiter).second).begin(); timeiter != ((*EobsTimeHigiter).second).end(); timeiter++) {
    // //     std::map<G4String, G4float>::iterator partobsiter;
    // //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
    // //       Eobs_byParticleAndGlobalTimeHig[(*EobsTimeHigiter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
    // //       EvtEnergybySDandpartandglobaltimehig[(*EobsTimeHigiter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
    // //     }
    // //   }
    // // }
    
    
    // // for (EobsTimeLowiter = Eobs_byPartAndLocTime1Low->begin(); EobsTimeLowiter != Eobs_byPartAndLocTime1Low->end(); EobsTimeLowiter++) {
    // //   std::map<G4int, std::map<G4String, G4float> >::iterator timeiter;
    // //   for (timeiter = ((*EobsTimeLowiter).second).begin(); timeiter != ((*EobsTimeLowiter).second).end(); timeiter++) {
    // //     std::map<G4String, G4float>::iterator partobsiter;
    // //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
    // //       Eobs_byParticleAndLocalTime1Low[(*EobsTimeLowiter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
    // //       EvtEnergybySDandpartandlocaltime1low[(*EobsTimeLowiter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
    // //     }
    // //   }
    // // }
    // // for (EobsTimeMediter = Eobs_byPartAndLocTime1Med->begin(); EobsTimeMediter != Eobs_byPartAndLocTime1Med->end(); EobsTimeMediter++) {
    // //   std::map<G4int, std::map<G4String, G4float> >::iterator timeiter;
    // //   for (timeiter = ((*EobsTimeMediter).second).begin(); timeiter != ((*EobsTimeMediter).second).end(); timeiter++) {
    // //     std::map<G4String, G4float>::iterator partobsiter;
    // //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
    // //       Eobs_byParticleAndLocalTime1Med[(*EobsTimeMediter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
    // //       EvtEnergybySDandpartandlocaltime1med[(*EobsTimeMediter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
    // //     }
    // //   }
    // // }
    // // for (EobsTimeHigiter = Eobs_byPartAndLocTime1Hig->begin(); EobsTimeHigiter != Eobs_byPartAndLocTime1Hig->end(); EobsTimeHigiter++) {
    // //   std::map<G4int, std::map<G4String, G4float> >::iterator timeiter;
    // //   for (timeiter = ((*EobsTimeHigiter).second).begin(); timeiter != ((*EobsTimeHigiter).second).end(); timeiter++) {
    // //     std::map<G4String, G4float>::iterator partobsiter;
    // //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
    // //       Eobs_byParticleAndLocalTime1Hig[(*EobsTimeHigiter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
    // //       EvtEnergybySDandpartandlocaltime1hig[(*EobsTimeHigiter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
    // //     }
    // //   }
    // // }
    
    
    // // for (EobsTimeLowiter = Eobs_byPartAndLocTime2Low->begin(); EobsTimeLowiter != Eobs_byPartAndLocTime2Low->end(); EobsTimeLowiter++) {
    // //   std::map<G4int, std::map<G4String, G4float> >::iterator timeiter;
    // //   for (timeiter = ((*EobsTimeLowiter).second).begin(); timeiter != ((*EobsTimeLowiter).second).end(); timeiter++) {
    // //     std::map<G4String, G4float>::iterator partobsiter;
    // //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
    // //       Eobs_byParticleAndLocalTime2Low[(*EobsTimeLowiter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
    // //       EvtEnergybySDandpartandlocaltime2low[(*EobsTimeLowiter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
    // //     }
    // //   }
    // // }
    // // for (EobsTimeMediter = Eobs_byPartAndLocTime2Med->begin(); EobsTimeMediter != Eobs_byPartAndLocTime2Med->end(); EobsTimeMediter++) {
    // //   std::map<G4int, std::map<G4String, G4float> >::iterator timeiter;
    // //   for (timeiter = ((*EobsTimeMediter).second).begin(); timeiter != ((*EobsTimeMediter).second).end(); timeiter++) {
    // //     std::map<G4String, G4float>::iterator partobsiter;
    // //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
    // //       Eobs_byParticleAndLocalTime2Med[(*EobsTimeMediter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
    // //       EvtEnergybySDandpartandlocaltime2med[(*EobsTimeMediter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
    // //     }
    // //   }
    // // }
    // // for (EobsTimeHigiter = Eobs_byPartAndLocTime2Hig->begin(); EobsTimeHigiter != Eobs_byPartAndLocTime2Hig->end(); EobsTimeHigiter++) {
    // //   std::map<G4int, std::map<G4String, G4float> >::iterator timeiter;
    // //   for (timeiter = ((*EobsTimeHigiter).second).begin(); timeiter != ((*EobsTimeHigiter).second).end(); timeiter++) {
    // //     std::map<G4String, G4float>::iterator partobsiter;
    // //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
    // //       Eobs_byParticleAndLocalTime2Hig[(*EobsTimeHigiter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
    // //       EvtEnergybySDandpartandlocaltime2hig[(*EobsTimeHigiter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
    // //     }
    // //   }
    // // }
    
    
    // //
    // // histogram the relative observed energy contribution by particle:
    // //  
    // //G4cout << ">>> histogram the relative observed energy contribution by particle:" << G4endl;
    
    // for (Eobsiter = Eobs_byPart->begin(); Eobsiter != Eobs_byPart->end(); Eobsiter++) {
    //   G4String volName = (*Eobsiter).first;
    //   VolumeDirset[volName]->cd();
    //   if (EvtObsEnergybySDandpart.find(volName) == EvtObsEnergybySDandpart.end()) {
    //     std::map<G4String, TH1F*> tmpmap;
    //     EvtObsEnergybySDandpart.insert(std::make_pair(volName, tmpmap));
    //     std::map<G4String, G4float>::iterator partiter;
    //     for (partiter = Eobs_byPart->at(volName).begin(); partiter != Eobs_byPart->at(volName).end(); partiter++) {
    //       G4String hname = "Obs" + (*partiter).first;
    //       const char* hisname = hname.c_str();
    //       TH1F *histo = new TH1F(hisname, "observed energy contribution by particle", 200, 0, 100);
    //       EvtObsEnergybySDandpart[volName].insert(std::make_pair((*partiter).first, histo));
    //       EvtObsEnergybySDandpart[volName].at((*partiter).first)->Fill(((*partiter).second * 0.1) / Ein);
    //     }
    //   } else {
    //     std::map<G4String, G4float>::iterator partiter;
    //     for (partiter = Eobs_byPart->at(volName).begin(); partiter != Eobs_byPart->at(volName).end(); partiter++) {
    //       G4String PName = (*partiter).first;
    //       if (EvtObsEnergybySDandpart[volName].find(PName) == EvtObsEnergybySDandpart[volName].end()) {
    //         G4String hname = "Obs" + PName;
    //         const char* hisname = hname.c_str();
    //         TH1F *histo = new TH1F(hisname, "observed Energy contribution by particle", 200, 0, 100);
    //         EvtObsEnergybySDandpart[volName].insert(std::make_pair(PName, histo));
    //         EvtObsEnergybySDandpart[volName].at(PName)->Fill(((*partiter).second * 0.1) / Ein);
    //       } else {
    //         EvtObsEnergybySDandpart[volName].at(PName)->Fill(((*partiter).second * 0.1) / Ein);
    //       }
    //     }
    //   }
    // }
    
    // std::map<G4String, std::map<G4String, G4float> >* NCeren_byPart = event->GetNCeren_byParticle();
    // std::map<G4String, std::map<G4String, G4float> >::iterator Niter;
    // for (Niter = NCeren_byPart->begin(); Niter != NCeren_byPart->end(); Niter++) {
    //   if (NCeren_byParticle.find((*Niter).first) == NCeren_byParticle.end()) {
    //     NCeren_byParticle[(*Niter).first] = (*Niter).second;
    //   } else {
    //     std::map<G4String, G4float>::iterator partiter;
    //     for (partiter = NCeren_byPart->at((*Niter).first).begin(); partiter != NCeren_byPart->at((*Niter).first).end(); partiter++) {
    //       NCeren_byParticle[(*Niter).first][(*partiter).first] = NCeren_byParticle[(*Niter).first][(*partiter).first] + (*partiter).second;
    //     }
    //   }
    // }
    
    // //
    // // histogram the relative NCeren contribution by particle:
    // //  
    // //G4cout << ">>> histogram the relative NCeren contribution by particle:" << G4endl;
      
    // for (Niter = NCeren_byPart->begin(); Niter != NCeren_byPart->end(); Niter++) {
    //   G4String volName = (*Niter).first;
    //   VolumeDirset[volName]->cd();
    //   if (EvtNcerenbySDandpart.find(volName) == EvtNcerenbySDandpart.end()) {
    //     std::map<G4String, TH1F*> tmpmap;
    //     EvtNcerenbySDandpart.insert(std::make_pair(volName, tmpmap));
    //     std::map<G4String, G4float>::iterator partiter;
    //     for (partiter = NCeren_byPart->at(volName).begin(); partiter != NCeren_byPart->at(volName).end(); partiter++) {
    //       G4String hname = "NC" + (*partiter).first;
    //       const char* hisname = hname.c_str();
    //       TH1F *histo = new TH1F(hisname, "nceren contribution by particle", 200, 0, 100);
    //       EvtNcerenbySDandpart[volName].insert(std::make_pair((*partiter).first, histo));
    //       EvtNcerenbySDandpart[volName].at((*partiter).first)->Fill(100.*(*partiter).second / event->GetNCeren());
    //     }
    //   } else {
    //     std::map<G4String, G4float>::iterator partiter;
    //     for (partiter = NCeren_byPart->at(volName).begin(); partiter != NCeren_byPart->at(volName).end(); partiter++) {
    //       G4String PName = (*partiter).first;
    //       if (EvtNcerenbySDandpart[volName].find(PName) == EvtNcerenbySDandpart[volName].end()) {
    //         G4String hname = "NC" + PName;
    //         const char* hisname = hname.c_str();
    //         TH1F *histo = new TH1F(hisname, "nceren contribution by particle", 200, 0, 100);
    //         EvtNcerenbySDandpart[volName].insert(std::make_pair(PName, histo));
    //         EvtNcerenbySDandpart[volName].at(PName)->Fill(100.*(*partiter).second / event->GetNCeren());
    //       } else {
    //         EvtNcerenbySDandpart[volName].at(PName)->Fill(100.*(*partiter).second / event->GetNCeren());
    //       }
    //     }
    //   }
    // }
    
    // std::map<G4String, std::map<G4String, G4int > >*procmult = event->GetProcessMult();
    // std::map<G4String, std::map<G4String, G4int> >::iterator prociter;
    // for (prociter = procmult->begin(); prociter != procmult->end(); prociter++) {
    //   if (processMult.find((*prociter).first) == processMult.end()) {
    //     processMult[(*prociter).first] = (*prociter).second;
    //   } else {
    //     std::map<G4String, G4int>::iterator partiter;
    //     for (partiter = procmult->at((*prociter).first).begin(); partiter != procmult->at((*prociter).first).end(); partiter++) {
    //       processMult[(*prociter).first][(*partiter).first] = processMult[(*prociter).first][(*partiter).first] + (*partiter).second;
    //     }
    //   }
    // }
    
    // std::map<G4String, std::map<G4String, std::map<G4String, products > > >*procmap = event->GetProcessMap();
    // std::map<G4String, std::map<G4String, std::map<G4String, products> > >::iterator SDiter; // iterator over sensitive detectors
    // std::map<G4String, std::map<G4String, products> >::iterator procmapiter; // Iterator over processes
    // std::map<G4String, products>::iterator partmpiter; // Iterator over particles
    
    // for (SDiter = procmap->begin(); SDiter != procmap->end(); SDiter++) { // check if sensitive det already included 
    //   if (processMap.find((*SDiter).first) == processMap.end()) {
    //     processMap[(*SDiter).first] = (*SDiter).second;
    //   } else {
    //     for (procmapiter = (*SDiter).second.begin(); procmapiter != (*SDiter).second.end(); procmapiter++) {
    //       if (processMap[(*SDiter).first].find((*procmapiter).first) == processMap[(*SDiter).first].end()) // check if process already included  
    //       {
    //         processMap[(*SDiter).first][(*procmapiter).first] = (*procmapiter).second;
    //       } else {
    //         for (partmpiter = (*procmapiter).second.begin(); partmpiter != (*procmapiter).second.end(); partmpiter++) {
    //           if (processMap[(*SDiter).first][(*procmapiter).first].find((*partmpiter).first)
    //               == processMap[(*SDiter).first][(*procmapiter).first].end()) {
    //             processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first] = (*partmpiter).second;
    //           } else {
    //             processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first].NParticles =
    //               processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first].NParticles
    //               + (*partmpiter).second.NParticles;
    //             processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first].kinE =
    //               processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first].kinE
    //               + (*partmpiter).second.kinE;
    //             processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first].totE =
    //               processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first].totE
    //               + (*partmpiter).second.totE;
    //           }
    //         }
    //       }
    //     }
    //   }
    // }
    
    
    // 
    // //-------------------------------------------
    // //
    // // Now we deal with the Hit Collections.
    // //
    // //-------------------------------------------
    // //G4cout << ">>> Now we deal with the Hit Collections." << G4endl;
    
    // std::map<G4String, std::vector<G4VHit*> >::iterator hciter;
    // for (hciter = hcMap->begin(); hciter != hcMap->end(); hciter++)
    // {
    //   std::vector<G4VHit*> hits = (*hciter).second;
    //   G4int NbHits = hits.size();
    //   //G4cout << "NbHits: " << NbHits << G4endl;
    //   std::vector<std::string> y = split((*hciter).first, '_');
    //   std::string Classname = y[1];
    //   std::string Colltype = y[2];
    //   std::string volName = y[0] + "_" + y[1];
    //   EvtEnergy[volName] = 0.0;
    //   EvtObsEnergy[volName] = 0.0;
    //   EvtNceren[volName] = 0.0;
      
    //   if (RunEnergy.find(volName) == RunEnergy.end()) {
    //     RunEnergy[volName] = 0.0;
    //   }
    //   if (RunObsEnergy.find(volName) == RunObsEnergy.end()) {
    //     RunObsEnergy[volName] = 0.0;
    //   }
    //   if (RunNceren.find(volName) == RunNceren.end()) {
    //     RunNceren[volName] = 0.0;
    //   }
    //   if (Classname == "DRCalorimeter") {
    //     for (G4int ii = 0; ii < NbHits; ii++) {
    //       DRCalorimeterHit* DRHit = dynamic_cast<DRCalorimeterHit*> (hits[ii]);
    //       tot_en = tot_en + DRHit->GetEdep();
    //       tot_nceren = tot_nceren + DRHit->GetNCeren();
    //       EvtEnergy[volName] = EvtEnergy[volName] + DRHit->GetEdep();
    //       EvtNceren[volName] = EvtNceren[volName] + DRHit->GetNCeren();
    //     }
    //   } else if (Classname == "DRTSCalorimeter" && Colltype == "HC") {
    //     for (G4int ii = 0; ii < NbHits; ii++) {
    //       DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
    //       tot_en = tot_en + DRHit->GetEdep();
    //       obs_en = obs_en + DRHit->GetEobsbirks();
    //       tot_nceren = tot_nceren + DRHit->GetNCeren();
    //       EvtEnergy[volName] = EvtEnergy[volName] + DRHit->GetEdep();
    //       EvtObsEnergy[volName] = EvtObsEnergy[volName] + DRHit->GetEobsbirks();
    //       EvtNceren[volName] = EvtNceren[volName] + DRHit->GetNCeren();
    //       h2_xy_firstTime_all -> Fill(DRHit->GetPos().x(),DRHit->GetPos().y(),DRHit->GetGlobalTime());
    //       h2_zy_firstTime_all -> Fill(DRHit->GetPos().z(),DRHit->GetPos().y(),DRHit->GetGlobalTime());
    //       h2_xy_firstTime[volName] -> Fill(DRHit->GetPos().x(),DRHit->GetPos().y(),DRHit->GetGlobalTime());
    //       h2_zy_firstTime[volName] -> Fill(DRHit->GetPos().z(),DRHit->GetPos().y(),DRHit->GetGlobalTime());
    //       //std::vector<DRTimeSliceHit>* tshc = DRHit->GetDRGlobalTimeSliceHitCol();
    //       // for (unsigned int iii = 0; iii < tshc->size(); iii++) {
    //       //   tedephisto->Fill((tshc->at(iii)).GetSlice(), (tshc->at(iii)).GetEdep());
    //       //   tedepEMhisto->Fill((tshc->at(iii)).GetSlice(), (tshc->at(iii)).GetEdepEM());
    //       //   tedepnEMhisto->Fill((tshc->at(iii)).GetSlice(), (tshc->at(iii)).GetEdepnonEM());
    //       //   tNChisto->Fill((tshc->at(i)).GetSlice(), (tshc->at(iii)).GetNCeren());
    //       // }
    //     }
    //   } else if (Classname == "Calorimeter") {
    //     for (G4int ii = 0; ii < NbHits; ii++) {
    //       CalorimeterHit* Hit = dynamic_cast<CalorimeterHit*> (hits[ii]);
    //       tot_en = tot_en + Hit->GetEdep();
    //     }
    //   } else if (Classname == "Tracker") {
    //     for (G4int ii = 0; ii < NbHits; ii++) {
    //       TrackerHit* THit = dynamic_cast<TrackerHit*> (hits[ii]);
    //       tot_en = tot_en + THit->GetEdep();
    //     }
    //   } else if (Classname == "PhotonDetector") {
    //     tot_nceren = tot_nceren + NbHits;
    //     G4cout << "Photon Detector hits: " << NbHits << G4endl;
    //     for (G4int ii = 0; ii < NbHits; ii++) {
    //       PhotonHit* PHit = dynamic_cast<PhotonHit*> (hits[ii]);
    //       PHit->Print(); // HJW just for now 
    //     }
    //   } else {
    //     //G4cout << "SD type: " << Classname << " unknown" << G4endl;
    //   }
    // } // end loop over Hit collections
    // 
    
    // std::map<G4String, G4float>::iterator e_itere;
    // for (e_itere = EvtEnergy.begin(); e_itere != EvtEnergy.end(); e_itere++) {
    //   G4String VolumeName = (*e_itere).first;
    //   if (EHistosbySD.find(VolumeName) == EHistosbySD.end()) {
    //     G4String hname = "E_" + (*e_itere).first;
    //     const char* hisname = hname.c_str();
    //     TH1F* tmphisto = new TH1F(hisname, "Total energy deposited in", 100, 0.0, Ein);
    //     EHistosbySD[VolumeName].insert(std::make_pair(VolumeName, tmphisto));
    //     EHistosbySD[VolumeName].at(VolumeName)->Fill(EvtEnergy[VolumeName]*0.001);
    //   } else {
    //     EHistosbySD[VolumeName].at(VolumeName)->Fill(EvtEnergy[VolumeName]*0.001);
    //   }
    // }
    // for (e_itere = EvtObsEnergy.begin(); e_itere != EvtObsEnergy.end(); e_itere++) {
    //   G4String VolumeName = (*e_itere).first;
    //   if (EObsHistosbySD.find(VolumeName) == EObsHistosbySD.end()) {
    //     G4String hname = "EObs_" + (*e_itere).first;
    //     const char* hisname = hname.c_str();
    //     TH1F* tmphisto = new TH1F(hisname, "Total observed energy deposited in", 100, 0.0, Ein);
    //     EObsHistosbySD[VolumeName].insert(std::make_pair(VolumeName, tmphisto));
    //     EObsHistosbySD[VolumeName].at(VolumeName)->Fill(EvtObsEnergy[VolumeName]*0.001);
    //   } else {
    //     EObsHistosbySD[VolumeName].at(VolumeName)->Fill(EvtObsEnergy[VolumeName]*0.001);
    //   }
    // }
    
    // for (e_itere = EvtNceren.begin(); e_itere != EvtNceren.end(); e_itere++) {
    //   G4String VolumeName = (*e_itere).first;
    //   if (NCHistosbySD.find(VolumeName) == NCHistosbySD.end()) {
    //     G4String hname = "NC_" + (*e_itere).first;
    //     const char* hisname = hname.c_str();
    //     TH1F* tmphisto = new TH1F(hisname, "Total Nr of Cerenkov photons in", 100, 0.0, 65000 * (Ein + 0.1 * Ein));
    //     NCHistosbySD[VolumeName].insert(std::make_pair(VolumeName, tmphisto));
    //     NCHistosbySD[VolumeName].at(VolumeName)->Fill(EvtObsEnergy[VolumeName]*0.001);
    //   } else {
    //     NCHistosbySD[VolumeName].at(VolumeName)->Fill(EvtObsEnergy[VolumeName]*0.001);
    //   }
    // }
    // std::map<G4String, G4float>::iterator e_iter;
    // for (e_iter = EvtEnergy.begin(); e_iter != EvtEnergy.end(); e_iter++) {
    //   RunEnergy[(*e_iter).first] = RunEnergy[(*e_iter).first]+(*e_iter).second;
    // }
    // std::map<G4String, G4float>::iterator o_iter;
    // for (o_iter = EvtObsEnergy.begin(); o_iter != EvtObsEnergy.end(); o_iter++) {
    //   RunObsEnergy[(*o_iter).first] = RunObsEnergy[(*o_iter).first]+(*o_iter).second;
    // }
    // std::map<G4String, G4float>::iterator c_iter;
    // for (c_iter = EvtObsEnergy.begin(); o_iter != EvtObsEnergy.end(); o_iter++) {
    //   RunObsEnergy[(*o_iter).first] = RunObsEnergy[(*o_iter).first]+(*o_iter).second;
    // }
  } // end loop over events
  
  
  /*
  G4cout << "==============================================================================" << G4endl;
  G4cout << "=========    Average deposited Energy/Event:   " << std::setw(12) << RunTotEnergy / G4float(nEvents) << "  [GeV]" << G4endl;
  G4cout << "=========    Average observed Energy/Event:    " << std::setw(12) << RunTotObsEnergy / G4float(nEvents) << "  [GeV]" << G4endl;
  G4cout << "==============================================================================" << G4endl;
  
  std::map<G4String, std::map<G4String, G4float> >::iterator Eiter;
  for (Eiter = E_byParticle.begin(); Eiter != E_byParticle.end(); Eiter++) {
    G4float sum = 0.0;
    G4cout << "-----------------------------------------------" << G4endl;
    G4cout << "Sensitive Volume:  " << (*Eiter).first << G4endl;
    G4cout << "Total Energy: (per Evt)     " << RunEnergy[(*Eiter).first] / G4float(nEvents) << G4endl;
    G4cout << "Observed Energy: (per Evt)  " << RunObsEnergy[(*Eiter).first] / G4float(nEvents) << G4endl;
    G4cout << "-----------------------------------------------" << G4endl;
    std::map<G4String, G4float> ::iterator p_iter;
    for (p_iter = E_byParticle[(*Eiter).first].begin(); p_iter != E_byParticle[(*Eiter).first].end(); p_iter++) {
      sum = sum + 100. * ((*p_iter).second / RunEnergy[(*Eiter).first]);
      G4cout << std::setw(20) << (*p_iter).first << std::setw(20) << 100. * ((*p_iter).second / RunEnergy[(*Eiter).first]) << " [%]" << G4endl;
    }
    G4cout << "==============================================================================" << G4endl;
    G4cout << std::setw(20) << "sum:" << std::setw(20) << sum << " [%]" << G4endl;
    G4cout << G4endl;
    G4cout << G4endl;
  }
  
  
  G4cout << "==============================================================================" << G4endl;
  G4cout << "=========    Eobs contribution by particle " << G4endl;
  G4cout << "==============================================================================" << G4endl;
  
  std::map<G4String, std::map<G4String, G4float> >::iterator Eobsiter;
  for (Eobsiter = Eobs_byParticle.begin(); Eobsiter != Eobs_byParticle.end(); Eobsiter++) {
    G4cout << "==============================================================================" << G4endl;
    G4cout << "Sensitive Volume:  " << (*Eobsiter).first << G4endl;
    G4cout << "==============================================================================" << G4endl;
    G4float sumObs = 0.0;
    std::map<G4String, G4float> ::iterator pobs_iter;
    for (pobs_iter = Eobs_byParticle[(*Eobsiter).first].begin(); pobs_iter != Eobs_byParticle[(*Eobsiter).first].end(); pobs_iter++) {
      sumObs = sumObs + 100. * ((*pobs_iter).second / RunObsEnergy[(*Eobsiter).first]);
      G4cout << std::setw(20) << (*pobs_iter).first << std::setw(20) << 100. * ((*pobs_iter).second / RunObsEnergy[(*Eobsiter).first]) << " [%]" << G4endl;
    }
    G4cout << "==============================================================================" << G4endl;
    G4cout << std::setw(20) << "sumObs:" << std::setw(20) << sumObs << " [%]" << G4endl;
    G4cout << G4endl;
    G4cout << G4endl;
  }
  
  if (RunTotNCeren > 0.0) {
    G4cout << "==============================================================================" << G4endl;
    G4cout << "========   NCeren contribution by particle          ==========================" << G4endl;
    G4cout << "==============================================================================" << G4endl;
    G4cout << "Average number of Cerenkov Photons/Event:" << RunTotNCeren / G4float(nEvents) << G4endl;
    G4float sumN = 0.0;
    for (Eiter = NCeren_byParticle.begin(); Eiter != NCeren_byParticle.end(); Eiter++) {
      G4cout << "Sensitive Volume:  " << (*Eiter).first << G4endl;
      std::map<G4String, G4float> ::iterator p_iter;
      G4cout << "-----------------------------------------------" << G4endl;
      for (p_iter = NCeren_byParticle[(*Eiter).first].begin(); p_iter != NCeren_byParticle[(*Eiter).first].end(); p_iter++) {
        sumN = sumN + 100. * ((*p_iter).second / RunTotNCeren);
        G4cout << std::setw(20) << (*p_iter).first << std::setw(20) << 100. * ((*p_iter).second / RunTotNCeren) << " [%]" << G4endl;
      }
    }
    G4cout << "==============================================================================" << G4endl;
    G4cout << std::setw(20) << "sumN:" << std::setw(20) << sumN << " [%]" << G4endl;
    G4cout << "==============================================================================" << G4endl;
    G4cout << G4endl;
    G4cout << G4endl;
    G4cout << G4endl;
  }
  
  
  
  G4cout << "==============================================================================" << G4endl;
  G4cout << "=========    Eobs contribution by particle and time (med) ====================" << G4endl;
  G4cout << "==============================================================================" << G4endl;
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4float> > >::iterator EobsTimeMediter;
  for (EobsTimeMediter = Eobs_byParticleAndGlobalTimeMed.begin(); EobsTimeMediter != Eobs_byParticleAndGlobalTimeMed.end(); ++EobsTimeMediter++) {
    G4cout << "==============================================================================" << G4endl;
    G4cout << "Sensitive Volume:  " << (*EobsTimeMediter).first << G4endl;
    G4cout << "==============================================================================" << G4endl;
    std::map<G4int, std::map<G4String, G4float> >::iterator timeiter;
    for (timeiter = ((*EobsTimeMediter).second).begin(); timeiter != ((*EobsTimeMediter).second).end(); ++timeiter)
    {
      G4cout << ">>> global time slice (med): " << (*timeiter).first << G4endl;
      G4float sumObs = 0.0;
      std::map<G4String, G4float> ::iterator pobs_iter;
      for (pobs_iter = ((*timeiter).second).begin(); pobs_iter != ((*timeiter).second).end(); pobs_iter++) {
        sumObs = sumObs + 100. * ((*pobs_iter).second / RunObsEnergy[(*EobsTimeMediter).first]);
        G4cout << std::setw(20) << (*pobs_iter).first << std::setw(20) << 100. * ((*pobs_iter).second / RunObsEnergy[(*EobsTimeMediter).first]) << " [%]" << G4endl;
      }
      G4cout << "==============================================================================" << G4endl;
      G4cout << std::setw(20) << "sumObs:" << std::setw(20) << sumObs << " [%]" << G4endl;
      G4cout << G4endl;
      G4cout << G4endl;
    }
  }
  
  // if (RunTotNCeren > 0.0) {
  //     G4cout << "==============================================================================" << G4endl;
  //     G4cout << "========   NCeren contribution by particle and time ==========================" << G4endl;
  //     G4cout << "==============================================================================" << G4endl;
  //     G4cout << "Average number of Cerenkov Photons/Event:" << RunTotNCeren / G4float(nEvents) << G4endl;
  //     G4float sumN = 0.0;
  //     for (Eiter = NCeren_byParticle.begin(); Eiter != NCeren_byParticle.end(); Eiter++) {
  //         G4cout << "Sensitive Volume:  " << (*Eiter).first << G4endl;
  //         std::map<G4String, G4float> ::iterator p_iter;
  //         G4cout << "-----------------------------------------------" << G4endl;
  //         for (p_iter = NCeren_byParticle[(*Eiter).first].begin(); p_iter != NCeren_byParticle[(*Eiter).first].end(); p_iter++) {
  //             sumN = sumN + 100. * ((*p_iter).second / RunTotNCeren);
  //             G4cout << std::setw(20) << (*p_iter).first << std::setw(20) << 100. * ((*p_iter).second / RunTotNCeren) << " [%]" << G4endl;
  //         }
  //     }
  //     G4cout << "==============================================================================" << G4endl;
  //     G4cout << std::setw(20) << "sumN:" << std::setw(20) << sumN << " [%]" << G4endl;
  //     G4cout << "==============================================================================" << G4endl;
  //     G4cout << G4endl;
  //     G4cout << G4endl;
  //     G4cout << G4endl;
  // }
  
    
  
  G4cout << "==============================================================================" << G4endl;
  G4cout << "========   Multiplicity of processes/event          ==========================" << G4endl;
  G4cout << "==============================================================================" << G4endl;
  std::map<G4String, G4int >::iterator processiter;
  std::map<G4String, std::map<G4String, G4int > >::iterator processMultiter;
  for (processMultiter = processMult.begin(); processMultiter != processMult.end(); processMultiter++) {
    G4cout << "-----------------------------------------------" << G4endl;
    G4cout << "Sensitive Volume:  " << (*processMultiter).first << G4endl;
    G4cout << "-----------------------------------------------" << G4endl;
    std::map<G4String, G4int> npartmap = (*processMultiter).second;
    std::map<G4String, G4int>::iterator npartiter;
    for (npartiter = npartmap.begin(); npartiter != npartmap.end(); npartiter++) {
      G4cout << std::setw(20) << (*npartiter).first << std::setw(20) << (*npartiter).second / G4float(nEvents) << G4endl;
    }
  }
  G4cout << "==============================================================================" << G4endl;
  G4cout << G4endl;
  G4cout << G4endl;
  G4cout << G4endl;
  G4cout << "==============================================================================" << G4endl;
  G4cout << "========   Particles/event produced by specific process         ==============" << G4endl;
  G4cout << "==============================================================================" << G4endl;
  std::map<G4String, std::map<G4String, std::map<G4String, products> > >::iterator SDiter1; // iterator over sensitive detectors
  std::map<G4String, std::map<G4String, products> >::iterator prociter1; // Iterator over processes
  std::map<G4String, products>::iterator pariter1; // Iterator over particles
  for (SDiter1 = processMap.begin(); SDiter1 != processMap.end(); SDiter1++) {
    G4cout << "SD:  " << (*SDiter1).first << G4endl;
    for (prociter1 = (*SDiter1).second.begin(); prociter1 != (*SDiter1).second.end(); prociter1++) {
      G4cout << "------------------------------------------------------------------------------" << G4endl;
      G4cout << "Process:  " << (*prociter1).first << G4endl;
      for (pariter1 = (*prociter1).second.begin(); pariter1 != (*prociter1).second.end(); pariter1++) {
        G4cout << "Particle:  " << std::setw(15) << (*pariter1).first
               << " #: " << std::setw(10) << (*pariter1).second.NParticles / G4float(nEvents)
               << " kinE: " << std::setw(10) << (*pariter1).second.kinE / G4float(nEvents)
               << " totE: " << std::setw(10) << (*pariter1).second.totE / G4float(nEvents)
               << G4endl;
      }
    }
  }
  */
  
  G4cout << "==============================================================================" << G4endl;
  G4cout << G4endl;
  G4cout << G4endl;
  G4cout << G4endl;
  
  G4cout << "===========================================" << G4endl;
  G4cout << " nr of bytes written:  " << outfile->Write() << G4endl;
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
