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

void GetProcessList(std::vector<G4String>* vec);
bool FillChain(TChain* chain, const std::string& inputFileList);
TH1F* GetCumulative(TProfile* prof, const int& normalize, const float& Eobs);



int main(int argc, char** argv)
{
  if( argc < 3 )
  {
    G4cout << "Program requires at least 2 arguments: input file list, base name of output file, maxEntries (-1), HCMap (false)" << G4endl;
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
  G4bool doHCMap = false;
  if( argc > 3 ) maxEntries = G4int(atoi(argv[3]));
  if( argc > 4 ) doHCMap = G4bool(atoi(argv[4]));
  
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
    
    if( !doHCMap ) Tevt -> SetBranchStatus("HCMap",0);
    
    
    // global variables
    std::map<G4String,TDirectory*> VolumeDirsetGlobal;
    std::map<G4String,TDirectory*> VolumeDirsetTime;
    
    std::vector<G4String>* volumeList;
    std::vector<G4String>* particleList  = new std::vector<G4String>;
    std::vector<G4String>* processList = new std::vector<G4String>;
    
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
    
    
    // define histograms
    std::map<G4String,TH1F*> h_Edep;
    std::map<G4String,TH1F*> h_Eobs;
    std::map<G4String,TH1F*> h_NCeren;
    
    std::map<G4String,std::map<G4String,TProfile*> > nonZeroEobs_byTime;                                // map<volume,map<timeSliceName,TProfile*> >
    std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > > nonZeroEobs_byParticleAndTime; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
    
    std::map<G4String,std::map<G4String,TProfile*> > Eobs_byTime;                                // map<volume,map<timeSliceName,TProfile*> >
    std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > > Eobs_byParticleAndTime; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
    std::map<G4String,std::map<G4String,TProfile*> > EobsFrac_byTime;                                // map<volume,map<timeSliceName,TProfile*> >
    std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > > EobsFrac_byParticleAndTime; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
    
    std::map<G4String,std::map<G4String,TProfile*> > Eobs_byTime_cumul;                                // map<volume,map<timeSliceName,TProfile*> >
    std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > > Eobs_byParticleAndTime_cumul; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
    std::map<G4String,std::map<G4String,TProfile*> > EobsFrac_byTime_cumul;                                // map<volume,map<timeSliceName,TProfile*> >
    std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > > EobsFrac_byParticleAndTime_cumul; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
    
    std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > > mult_byProcessAndTime; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
    std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > > Etot_byProcessAndTime; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
    std::map<G4String,std::map<G4String,std::map<G4String,TProfile*> > > Ekin_byProcessAndTime; // map<volume,map<timeSliceName,map<particle,TProfile*> > >
    
    std::map<G4String,std::map<G4String,TProfile2D*> > Eobs_time_vs_z; // map<volume,map<timeSliceName,TProfile2D*> >
    std::map<G4String,std::map<G4String,TProfile2D*> > Eobs_time_vs_R; // map<volume,map<timeSliceName,TProfile2D*> >
    std::map<G4String,std::map<G4String,TProfile2D*> > EobsFrac_time_vs_z; // map<volume,map<timeSliceName,TProfile2D*> >
    std::map<G4String,std::map<G4String,TProfile2D*> > EobsFrac_time_vs_R; // map<volume,map<timeSliceName,TProfile2D*> >
    
    // std::map<G4String,TProfile2D*> h2_xy_firstTime;
    // std::map<G4String,TProfile2D*> h2_zy_firstTime;
    
    G4float RunTotDepEnergy = 0.0;
    G4float RunTotObsEnergy = 0.0;
    G4float RunTotNCeren = 0.0;
    
    G4float Ein;
    G4float Edep;
    G4float Eobs;
    G4float NCeren;
    
    
    // initialize histograms    
    Trh->GetEntry(0);
    //runHeader->Print();
    
    Ein = runHeader->GetParticleEnergy();
    
    volumeList = runHeader -> GetVolumes();
    volumeList->push_back("AllVol");
    particleList = runHeader -> GetParticleList();
    GetProcessList(processList);
    
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
        VolumeDirsetTime[volName]->mkdir("particleHistos");
        VolumeDirsetTime[volName]->mkdir("processHistos");
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
          
          VolumeDirsetTime[volName]->cd();
          
          nonZeroEobs_byTime[volName][timeSliceName] = new TProfile(Form("nonZeroEobs%s_%s",    volName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
          
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
            
            VolumeDirsetTime[volName]->cd("particleHistos");
            
            nonZeroEobs_byParticleAndTime[volName][timeSliceName][partName] = new TProfile(Form("nonZeroEobs%s_%s_%s", volName.c_str(),partName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
            
            Eobs_byParticleAndTime[volName][timeSliceName][partName]     = new TProfile(Form("Eobs%s_%s_%s",    volName.c_str(),partName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
            EobsFrac_byParticleAndTime[volName][timeSliceName][partName] = new TProfile(Form("EobsFrac%s_%s_%s",volName.c_str(),partName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
            
            Eobs_byParticleAndTime_cumul[volName][timeSliceName][partName]     = new TProfile(Form("Eobs%s_cumul_%s_%s",    volName.c_str(),partName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
            EobsFrac_byParticleAndTime_cumul[volName][timeSliceName][partName] = new TProfile(Form("EobsFrac%s_cumul_%s_%s",volName.c_str(),partName.c_str(),timeSliceName.c_str()),"Eobs",nBins,min,max,"s");
          }
          
          VolumeDirsetTime[volName]->cd();
          
          for(unsigned int procIt = 0; procIt < processList->size(); ++procIt)
          {
            std::string procName = processList->at(procIt);
            
            VolumeDirsetTime[volName]->cd("processHistos");
            
            mult_byProcessAndTime[volName][timeSliceName][procName] = new TProfile(Form("mult%s_%s_%s",volName.c_str(),procName.c_str(),timeSliceName.c_str()),"mult",nBins,min,max,"s");
            Etot_byProcessAndTime[volName][timeSliceName][procName] = new TProfile(Form("Etot%s_%s_%s",volName.c_str(),procName.c_str(),timeSliceName.c_str()),"Etot",nBins,min,max,"s");
            Ekin_byProcessAndTime[volName][timeSliceName][procName] = new TProfile(Form("Ekin%s_%s_%s",volName.c_str(),procName.c_str(),timeSliceName.c_str()),"Ekin",nBins,min,max,
"s");
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
      //std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* m_Edep_byTime   = event->GetEdepByTimeMap();   // map<Detector,map<timeSliceType,map<timeSlice,energy> > >
      std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* m_Eobs_byTime   = event->GetEobsByTimeMap();   // map<Detector,map<timeSliceType,map<timeSlice,energy> > >
      //std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* m_NCeren_byTime = event->GetNCerenByTimeMap(); // map<Detector,map<timeSliceType,map<timeSlice,energy> > >
      std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > >* m_Eobs_byParticleAndTime = event->GetEobsByParticleAndTimeMap(); // map<Detector,map<timeSliceType,map<timeSlice,map<particle,energy> > > >
      std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4int> > > >* m_mult_byProcessAndTime = event->GetProcessAndTimeMult(); // map<Detector,map<timeSliceType,map<timeSlice,map<processes,mult> > > >
      std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,std::map<G4String,products> > > > >* m_products_byProcessAndTime = event->GetProcessAndTimeMap(); // map<Detector,map<timeSliceType,map<timeSlice,map<process,map< particle,products> > > > >
      std::map<G4String,std::map<G4String,std::map<G4ThreeVector,std::vector<G4VHit*> > > >* HCMap;
      if( doHCMap ) HCMap = event->GetHCMap();
      
      
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
          
          
          //-----------------------------------
          // global, particle and process plots
          std::map<G4String,G4float> sum;
          G4float sumAllVol = 0.;
          std::map<G4String,std::map<G4String,G4float> > sumPart;
          std::map<G4String,G4float> sumPartAllVol;
          for(int bin = 0; bin < nBins; ++bin)
          {
            G4float time = min + (bin-1)*timeSliceSize + 0.5*timeSliceSize;
            G4float volSum = 0.;
            G4float volSum2 = 0.;
            G4float volMultSum = 0.;
            for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
            {
              G4String volName = volumeList->at(volIt);
              sum[volName] += (*m_Eobs_byTime)[volName][timeSliceName][bin];
              sumAllVol += (*m_Eobs_byTime)[volName][timeSliceName][bin];
              volSum += (*m_Eobs_byTime)[volName][timeSliceName][bin];
              
              if( (*m_Eobs_byTime)[volName][timeSliceName][bin] > 0. )
                nonZeroEobs_byTime[volName][timeSliceName] -> Fill( time,(*m_Eobs_byTime)[volName][timeSliceName][bin] );
              
              Eobs_byTime[volName][timeSliceName]      -> Fill( time,(*m_Eobs_byTime)[volName][timeSliceName][bin] );
              EobsFrac_byTime[volName][timeSliceName]  -> Fill( time,(*m_Eobs_byTime)[volName][timeSliceName][bin] / Eobs );
              
              Eobs_byTime[volName][timeSliceName]      -> Fill( time,(*m_Eobs_byTime)[volName][timeSliceName][bin] );
              EobsFrac_byTime[volName][timeSliceName]  -> Fill( time,(*m_Eobs_byTime)[volName][timeSliceName][bin] / Eobs );
              
              Eobs_byTime_cumul[volName][timeSliceName]      -> Fill( time,sum[volName] );
              EobsFrac_byTime_cumul[volName][timeSliceName]  -> Fill( time,sum[volName] / Eobs );
            }
            
            Eobs_byTime["AllVol"][timeSliceName]     -> Fill( time,volSum );
            EobsFrac_byTime["AllVol"][timeSliceName] -> Fill( time,volSum / Eobs );
            
            Eobs_byTime_cumul["AllVol"][timeSliceName]     -> Fill( time,sumAllVol );
            EobsFrac_byTime_cumul["AllVol"][timeSliceName] -> Fill( time,sumAllVol / Eobs );
            
            for(unsigned int partIt = 0; partIt < particleList->size(); ++partIt)
            {
              std::string partName = particleList->at(partIt);
              volSum = 0.;
              for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
              {
                G4String volName = volumeList->at(volIt);
                sumPart[volName][partName] += (*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName];
                sumPartAllVol[partName] += (*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName];
                volSum += (*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName];
                
                if( (*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName] > 0. )
                  nonZeroEobs_byParticleAndTime[volName][timeSliceName][partName] -> Fill( time,(*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName] );
                
                Eobs_byParticleAndTime[volName][timeSliceName][partName]      -> Fill( time,(*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName] );
                EobsFrac_byParticleAndTime[volName][timeSliceName][partName]  -> Fill( time,(*m_Eobs_byParticleAndTime)[volName][timeSliceName][bin][partName] / Eobs );
                
                Eobs_byParticleAndTime_cumul[volName][timeSliceName][partName]      -> Fill( time,sumPart[volName][partName] );
                EobsFrac_byParticleAndTime_cumul[volName][timeSliceName][partName]  -> Fill( time,sumPart[volName][partName] / Eobs );
              }
              
              Eobs_byParticleAndTime["AllVol"][timeSliceName][partName]     -> Fill( time,volSum );
              EobsFrac_byParticleAndTime["AllVol"][timeSliceName][partName] -> Fill( time,volSum / Eobs );
              
              Eobs_byParticleAndTime_cumul["AllVol"][timeSliceName][partName]     -> Fill( time,sumPartAllVol[partName] );
              EobsFrac_byParticleAndTime_cumul["AllVol"][timeSliceName][partName] -> Fill( time,sumPartAllVol[partName] / Eobs );
            }
            
            for(unsigned int procIt = 0; procIt < processList->size(); ++procIt)
            {
              std::string procName = processList->at(procIt);
              volSum = 0.;
              volSum2 = 0.;
              volMultSum = 0.;
              for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
              {
                G4String volName = volumeList->at(volIt);
                
                std::map<G4String,products> tmpMap = (*m_products_byProcessAndTime)[volName][timeSliceName][bin][procName];
                G4float tmpSum = 0.;
                G4float tmpSum2 = 0.;
                for(std::map<G4String,products>::const_iterator tmpMapIt = tmpMap.begin(); tmpMapIt != tmpMap.end(); ++tmpMapIt)
                {
                  tmpSum  += tmpMapIt->second.totE;
                  tmpSum2 += tmpMapIt->second.kinE;
                }
                
                volMultSum += (*m_mult_byProcessAndTime)[volName][timeSliceName][bin][procName];
                volSum += tmpSum;
                volSum2 += tmpSum2;
                
                mult_byProcessAndTime[volName][timeSliceName][procName] -> Fill( time,(*m_mult_byProcessAndTime)[volName][timeSliceName][bin][procName] );
                Etot_byProcessAndTime[volName][timeSliceName][procName] -> Fill( time,tmpSum );
                Ekin_byProcessAndTime[volName][timeSliceName][procName] -> Fill( time,tmpSum2 );
              }
              
              mult_byProcessAndTime["AllVol"][timeSliceName][procName] -> Fill( time,volMultSum );
              Etot_byProcessAndTime["AllVol"][timeSliceName][procName] -> Fill( time,volSum );
              Ekin_byProcessAndTime["AllVol"][timeSliceName][procName] -> Fill( time,volSum2 );
            }
          }
          
          
          //------------
          // HCMap plots
          if( doHCMap )
          {
            for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
            {
              G4String volName = volumeList->at(volIt);
              
              std::map<G4ThreeVector,std::vector<G4VHit*> > posHitMap = (*HCMap)[volName][timeSliceName];
              for(std::map<G4ThreeVector,std::vector<G4VHit*> >::const_iterator mapIt = posHitMap.begin(); mapIt != posHitMap.end(); ++mapIt)
              {
                double x = (mapIt->first).x();
                double y = (mapIt->first).y();
                double z = (mapIt->first).z();
                std::vector<G4VHit*> hitVec = mapIt->second;
                for(unsigned int vecIt = 0; vecIt < hitVec.size(); ++vecIt)
                {
                  DRTSCalorimeterHit2* aHit = dynamic_cast<DRTSCalorimeterHit2*>(hitVec.at(vecIt));
                  G4int timeSlice = aHit -> GetTimeSlice();
                  G4float time = min + (timeSlice-1)*timeSliceSize + 0.5*timeSliceSize;
                  
                  Eobs_time_vs_z[volName][timeSliceName] -> Fill( z,time,aHit->GetEdep() );
                  Eobs_time_vs_R[volName][timeSliceName] -> Fill( sqrt(x*x+y*y+z*z),time,aHit->GetEobsbirks() );
                  EobsFrac_time_vs_z[volName][timeSliceName] -> Fill( z,time,aHit->GetEdep()/Eobs );
                  EobsFrac_time_vs_R[volName][timeSliceName] -> Fill( sqrt(x*x+y*y+z*z),time,aHit->GetEobsbirks()/Eobs );
                }
              }
            }
          } // HCMap plots
          
          
        }
      }
      
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



void GetProcessList(std::vector<G4String>* vec)
{
  vec -> push_back("alphaInelastic");
  vec -> push_back("annihil");
  vec -> push_back("anti-lambdaInelastic");
  vec -> push_back("anti_neutronInelastic");
  vec -> push_back("anti_omega-Inelastic");
  vec -> push_back("anti_protonInelastic");
  vec -> push_back("anti_sigma-Inelastic");
  vec -> push_back("anti_sigma+Inelastic");
  vec -> push_back("anti_xi0Inelastic");
  vec -> push_back("anti_xi-Inelastic");
  vec -> push_back("compt");
  vec -> push_back("conv");
  vec -> push_back("CoulombScat");
  vec -> push_back("Decay");
  vec -> push_back("dInelastic");
  vec -> push_back("eBrem");
  vec -> push_back("eIoni");
  vec -> push_back("electronNuclear");
  vec -> push_back("hadElastic");
  vec -> push_back("hBertiniCaptureAtRest");
  vec -> push_back("hBrems");
  vec -> push_back("He3Inelastic");
  vec -> push_back("hFritiofCaptureAtRest");
  vec -> push_back("hIoni");
  vec -> push_back("hPairProd");
  vec -> push_back("ionInelastic");
  vec -> push_back("ionIoni");
  vec -> push_back("kaon0LInelastic");
  vec -> push_back("kaon0SInelastic");
  vec -> push_back("kaon-Inelastic");
  vec -> push_back("kaon+Inelastic");
  vec -> push_back("lambdaInelastic");
  vec -> push_back("msc");
  vec -> push_back("muBrems");
  vec -> push_back("muIoni");
  vec -> push_back("muMinusCaptureAtRest");
  vec -> push_back("muonNuclear");
  vec -> push_back("muPairProd");
  vec -> push_back("nCapture");
  vec -> push_back("neutronInelastic");
  vec -> push_back("nKiller");
  vec -> push_back("omega-Inelastic");
  vec -> push_back("phot");
  vec -> push_back("photonNuclear");
  vec -> push_back("pi-Inelastic");
  vec -> push_back("pi+Inelastic");
  vec -> push_back("positronNuclear");
  vec -> push_back("protonInelastic");
  vec -> push_back("sigma-Inelastic");
  vec -> push_back("sigma+Inelastic");
  vec -> push_back("tInelastic");
  vec -> push_back("Transportation");
  vec -> push_back("xi0Inelastic");
  vec -> push_back("xi-Inelastic");
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
