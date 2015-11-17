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

double NCerenCalib = 5.70084e+04;



int main(int argc, char** argv)
{
  if( argc < 3 )
  {
    G4cout << "Program requires at least 4 arguments: input file list, base name of output file, particle, computeCorr (1), applyCorr (0), maxEntries (-1)" << G4endl;
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
  // initialize variables
  
  std::vector<int> energies;
  energies.push_back(3);
  energies.push_back(10);
  energies.push_back(30);
  energies.push_back(50);
  energies.push_back(100);
  energies.push_back(300);
  
  
  
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
  G4String particle = "";
  G4int computeCorr = 1;
  G4int applyCorr = 0;
  G4int maxEntries = -1;
  
  if( argc > 3 ) particle    = G4String(argv[3]);
  if( argc > 4 ) computeCorr = G4int(atoi(argv[4]));
  if( argc > 5 ) applyCorr   = G4int(atoi(argv[5]));
  if( argc > 6 ) maxEntries  = G4int(atoi(argv[6]));
  
  
  while(1)
  {
    std::getline(inputFileList,buffer);
    if( !inputFileList.good() ) break;
    if( buffer.at(0) == '#' ) continue;
    
    std::string outFileName = "";
    if( applyCorr == 0 ) outFileName = Form("%s_%d.root",argv[2],fileIt);
    else                 outFileName = Form("%s_corr_%d.root",argv[2],fileIt);
    outFile = new TFile(outFileName.c_str(), "RECREATE");
    
    TFile* inFile = TFile::Open(buffer.c_str(),"READ");
    if( !inFile ) continue;
    else std::cout << ">>> file " << buffer << " opened" << std::endl;
    
    TTree* Tevt = (TTree*)( inFile->Get("EventTree") );
    TTree* Trh  = (TTree*)( inFile->Get("RunTree") );
    
    if( !Tevt->GetListOfBranches()->FindObject("Event") )
    {
      std::cout << ">>> Event not found <<<" << std::endl;
      continue;
    }
    if( !Trh->GetListOfBranches()->FindObject("RunHeader") )
    {
      std::cout << ">>> RunHeader not found <<<" << std::endl;
      continue;
    }
    
    Event* event = new Event();
    TBranch* b_event = Tevt -> GetBranch("Event");
    b_event -> SetAddress(&event);
    
    RunHeader* runHeader = new RunHeader();
    TBranch* b_runHeader = Trh -> GetBranch("RunHeader");
    b_runHeader -> SetAddress(&runHeader);
    
    Trh->GetEntry(0);
    runHeader->Print();
    
    
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
    
    std::vector<G4float> intTimeVals;
    intTimeVals.push_back(0.100); //  100 ps
    intTimeVals.push_back(0.200); //  200 ps
    intTimeVals.push_back(0.300); //  300 ps
    intTimeVals.push_back(0.500); //  500 ps
    intTimeVals.push_back(0.750); //  750 ps
    intTimeVals.push_back(1.000); // 1.00 ns
    intTimeVals.push_back(1.250); // 1.25 ns
    intTimeVals.push_back(1.500); // 1.50 ns
    intTimeVals.push_back(1.750); // 1.75 ns
    intTimeVals.push_back(2.000); // 2.00 ns
    intTimeVals.push_back(2.250); // 2.25 ns
    intTimeVals.push_back(2.500); // 2.50 ns
    intTimeVals.push_back(2.750); // 2.75 ns
    intTimeVals.push_back(3.000); // 3.00 ns
    intTimeVals.push_back(3.500); // 3.50 ns
    intTimeVals.push_back(4.000); // 4.00 ns
    intTimeVals.push_back(4.500); // 4.50 ns
    
    
    
    // define histograms
    std::map<G4String,TH1F*> h_Edep;
    std::map<G4String,TH1F*> h_Eobs;
    std::map<G4String,TH1F*> h_NCeren;
    
    std::map<G4String,TH1F*> h_Eobs_corr1;
    std::map<G4String,TH1F*> h_Eobs_corr2_cer;
    std::map<G4String,TH1F*> h_Eobs_corr3_cer;
    std::map<G4String,std::map<G4String,std::map<G4float,TH1F*> > > h_Eobs_corr2_time;
    std::map<G4String,std::map<G4String,std::map<G4float,TH1F*> > > h_Eobs_corr3_time;
    
    std::map<G4String,TH2F*> h2_EobsMaxTime_vs_NCeren; // map<volume,TH2F*>
    std::map<G4String,TProfile*> p_EobsMaxTime_vs_NCeren; // map<volume,TProfile*>
    std::map<G4String,TH2F*> h2_EobsMaxTime_vs_NCerenCalib; // map<volume,TH2F*>
    std::map<G4String,TProfile*> p_EobsMaxTime_vs_NCerenCalib; // map<volume,TProfile*>
    
    std::map<G4String,std::map<G4String,std::map<G4float,TH2F*> > > h2_EobsMaxTime_vs_EobsIntTime; // map<volume,map<timeType,map<timeVal,TH2F*> > >
    std::map<G4String,std::map<G4String,std::map<G4float,TProfile*> > > p_EobsMaxTime_vs_EobsIntTime; // map<volume,map<timeType,map<timeVal,TProfile*> > >
    
    std::map<G4String,std::map<int,TF1*> > f_EobsMaxTime_vs_NCerenCalib;
    std::map<G4String,std::map<G4String,std::map<int,std::map<G4float,TF1*> > > > f_EobsMaxTime_vs_EobsIntTime;
    
    
    // initialize histograms    
    
    Ein = runHeader->GetParticleEnergy();
    std::cout << "Ein: " << Ein << std::endl;
    
    volumeList = runHeader -> GetVolumes();
    volumeList->push_back("AllVol");
    
    // timeSliceSizes["Low"] = runHeader->GetTimeSliceSizeLow();
    // timeSliceSizes["Med"] = runHeader->GetTimeSliceSizeMed();
    // timeSliceSizes["Hig"] = runHeader->GetTimeSliceSizeHig();
    // minTimes["Low"] = runHeader->GetMinTimeLow();
    // minTimes["Med"] = runHeader->GetMinTimeMed();
    // minTimes["Hig"] = runHeader->GetMinTimeHig();
    // maxTimes["Low"] = runHeader->GetMaxTimeLow();
    // maxTimes["Med"] = runHeader->GetMaxTimeMed();
    // maxTimes["Hig"] = runHeader->GetMaxTimeHig();
    timeSliceSizes["Low"] = 0.010;
    timeSliceSizes["Med"] = 0.2;
    timeSliceSizes["Hig"] = 10.;
    minTimes["Low"] = 0.;
    minTimes["Med"] = 5.;
    minTimes["Hig"] = 0.;
    maxTimes["Low"] = 100.;
    maxTimes["Med"] = 0.;
    maxTimes["Hig"] = 5000.;
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
      
      h_Edep[volName]   = new TH1F(Form("h_Edep%s",volName.c_str()),   "Total energy deposited",                   2000, 0., std::max(2.*Ein,0.05));
      h_Eobs[volName]   = new TH1F(Form("h_Eobs%s",volName.c_str()),   "Total observed (Birks suppressed) energy", 2000, 0., std::max(2.*Ein,0.05));
      h_NCeren[volName] = new TH1F(Form("h_NCeren%s",volName.c_str()), "Total nr. of cerenkov photons",            2000, 0., 65000*1.1*Ein);
      
      h_Eobs_corr1[volName] = new TH1F(Form("h_Eobs%s_corr1",volName.c_str()), "Total observed (Birks suppressed) energy", 2000, 0., std::max(2.*Ein,0.05));
      h_Eobs_corr2_cer[volName] = new TH1F(Form("h_Eobs%s_corr2_cer",volName.c_str()), "Total observed (Birks suppressed) energy", 2000, 0., std::max(2.*Ein,0.05));
      h_Eobs_corr3_cer[volName] = new TH1F(Form("h_Eobs%s_corr3_cer",volName.c_str()), "Total observed (Birks suppressed) energy", 2000, 0., std::max(2.*Ein,0.05));
      
      h2_EobsMaxTime_vs_NCeren[volName] = new TH2F(Form("h2_EobsMaxTime_vs_NCeren%s",volName.c_str()),"",2000, 0., 65000*1.1*Ein,1000,0.,1.);
      p_EobsMaxTime_vs_NCeren[volName] = new TProfile(Form("p_EobsMaxTime_vs_NCeren%s",volName.c_str()),"",2000, 0., 65000*1.1*Ein);
      h2_EobsMaxTime_vs_NCerenCalib[volName] = new TH2F(Form("h2_EobsMaxTime_vs_NCerenCalib%s",volName.c_str()),"",4000, 0., 2.,1000,0.,1.);
      p_EobsMaxTime_vs_NCerenCalib[volName] = new TProfile(Form("p_EobsMaxTime_vs_NCerenCalib%s",volName.c_str()),"",4000, 0., 2.);
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
        
        for(unsigned int intTimeIt = 0; intTimeIt < intTimeVals.size(); ++intTimeIt)
        {
          G4float intTime = intTimeVals.at(intTimeIt);
          G4String name = Form("intTime%s_%s_%.3fns",volName.c_str(),timeSliceName.c_str(),intTime);
          
          h2_EobsMaxTime_vs_EobsIntTime[volName][timeSliceName][intTime] = new TH2F(Form("h2_%s",name.c_str()),"",1000,0.,1.,1000,0.,1.);
          p_EobsMaxTime_vs_EobsIntTime[volName][timeSliceName][intTime] = new TProfile(Form("p_%s",name.c_str()),"",1000,0.,1.);
          
          h_Eobs_corr2_time[volName][timeSliceName][intTime] = new TH1F(Form("h_Eobs%s_corr2_time",name.c_str()), "Total observed (Birks suppressed) energy", 2000, 0., std::max(2.*Ein,0.05));
          h_Eobs_corr3_time[volName][timeSliceName][intTime] = new TH1F(Form("h_Eobs%s_corr3_time",name.c_str()), "Total observed (Birks suppressed) energy", 2000, 0., std::max(2.*Ein,0.05));
        }
      }
    }
    
    
    //-----------------
    // loop over events
    
    if( computeCorr )
    {
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
        std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* m_Edep_byTime = event->GetEdepByTimeMap();
        
        
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
          
          h2_EobsMaxTime_vs_NCeren[volName] -> Fill( (*m_NCeren)[volName],(*m_Eobs)[volName]/Ein );
          p_EobsMaxTime_vs_NCeren[volName]  -> Fill( (*m_NCeren)[volName],(*m_Eobs)[volName]/Ein );
          
          h2_EobsMaxTime_vs_NCerenCalib[volName] -> Fill( (*m_NCeren)[volName]/(NCerenCalib*(*m_Eobs)[volName]),(*m_Eobs)[volName]/Ein );
          p_EobsMaxTime_vs_NCerenCalib[volName]  -> Fill( (*m_NCeren)[volName]/(NCerenCalib*(*m_Eobs)[volName]),(*m_Eobs)[volName]/Ein );
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
            G4String timeType = timeTypes.at(timeTypeIt);
            G4String timeSliceName = timeType + "Low";
            
            
            // sum Eobs
            float EobsMaxTime = 0.;
            std::map<float,float> EobsIntTime;
            for(int bin = 1; bin <= nTimeSlices["Low"]+1; ++bin)
            {
              EobsMaxTime += (*m_Edep_byTime)[volName][timeSliceName][bin];
              
              for(unsigned int intTimeIt = 0; intTimeIt < intTimeVals.size(); ++intTimeIt)
              {
                G4float intTime = intTimeVals.at(intTimeIt);
                
                if( minTimes["Low"]+(bin-1)*timeSliceSizes["Low"] < intTime )
                  EobsIntTime[intTime] += (*m_Edep_byTime)[volName][timeSliceName][bin];
              }
            }
            
            
            // fill histo
            for(unsigned int intTimeIt = 0; intTimeIt < intTimeVals.size(); ++intTimeIt)
            {
              G4float intTime = intTimeVals.at(intTimeIt);
              
              h2_EobsMaxTime_vs_EobsIntTime[volName][timeSliceName][intTime] -> Fill( EobsIntTime[intTime]/EobsMaxTime,EobsMaxTime/Ein );
              p_EobsMaxTime_vs_EobsIntTime[volName][timeSliceName][intTime] -> Fill( EobsIntTime[intTime]/EobsMaxTime,EobsMaxTime/Ein );
            }
          }
        }
        
      } // end loop over events
      std::cout << "\n>>> loop over events done" << std::endl;
    }
    
    
    if( applyCorr )
    {
      //-------------------
      // compute correction
      
      std::map<G4String,float> corr1;
      TFile* inFile_corr = TFile::Open("/afs/cern.ch/work/a/abenagli/CaTS/CaTS_source/dualReadout__tiledsamplingcal_1000x1000x1mmFe_rindex__QGSP_BERT_HP_plots.root","READ");
      
      for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
      {
        G4String volName = volumeList->at(volIt);
        
        h_Eobs[volName] = (TH1F*)( inFile_corr->Get(Form("h_Eobs%s_%s_%dGeV",volName.c_str(),particle.c_str(),G4int(Ein))) );
        corr1[volName] = Ein / h_Eobs[volName] -> GetMean();
        
        for(unsigned int timeTypeIt = 0; timeTypeIt < timeTypes.size(); ++timeTypeIt)
        {
          G4String timeType = timeTypes.at(timeTypeIt);
          G4String timeSliceName = timeType + "Low";
          
          for(unsigned int energyIt = 0; energyIt < energies.size(); ++energyIt)
          {
            int energy = energies.at(energyIt);
            
            for(unsigned int intTimeIt = 0; intTimeIt < intTimeVals.size(); ++intTimeIt)
            {
              G4float intTime = intTimeVals.at(intTimeIt);
              
              TF1* fitFunc = (TF1*)( inFile_corr->Get(Form("fitFunc_cer_tiledsamplingcal_1000x1000x%dmm%s__%s__%dGeV",
                                                           1,"Fe","QGSP_BERT_HP",energy)) );
              f_EobsMaxTime_vs_NCerenCalib[volName][energy] = fitFunc;
              
              fitFunc = (TF1*)( inFile_corr->Get(Form("fitFunc_%s_tiledsamplingcal_1000x1000x%dmm%s__%s_%3dps__%dGeV",
                                                      timeSliceName.c_str(),1,"Fe","QGSP_BERT_HP",int(1000*intTime),energy)) );
              f_EobsMaxTime_vs_EobsIntTime[volName][timeSliceName][energy][intTime] = fitFunc;
            }
          }
        }
      } // end compute correction
      std::cout << "\n>>> compute correction done" << std::endl;
      
      
      
      //-------------------
      // loop over events 2
      
      G4int nEvents = Tevt->GetEntries();      
      G4cout << "Nr. of Events:  " << nEvents << std::endl;
      for(G4int entry = 0; entry < nEvents; ++entry)
      {
        std::cout << ">>> processing event " << entry << " / " << nEvents << "\r" << std::flush;
        Tevt -> GetEntry(entry);
        
        
        Eobs = event->GetTotObsEnergy();
        
        std::map<G4String,G4float>* m_Edep   = event->GetEdepMap();   // map<Detector,energy>
        std::map<G4String,G4float>* m_Eobs   = event->GetEobsMap();   // map<Detector,energy>
        std::map<G4String,G4float>* m_NCeren = event->GetNCerenMap(); // map<Detector,energy>
        std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* m_Edep_byTime = event->GetEdepByTimeMap();
        
        
        // ------------
        // global plots
        
        for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
        {
          G4String volName = volumeList->at(volIt);
          float x = (*m_NCeren)[volName]/(NCerenCalib*(*m_Eobs)[volName]);
          
          h_Edep[volName] -> Fill((*m_Edep)[volName]);
          h_Eobs[volName] -> Fill((*m_Eobs)[volName]);
          h_Eobs_corr1[volName] -> Fill( (*m_Eobs)[volName]*corr1[volName] );
          
          if( f_EobsMaxTime_vs_NCerenCalib[volName][Ein] != NULL )
            h_Eobs_corr2_cer[volName] -> Fill( (*m_Eobs)[volName] / f_EobsMaxTime_vs_NCerenCalib[volName][Ein]->Eval(x) );
          
          if( f_EobsMaxTime_vs_NCerenCalib[volName][50] != NULL )
            h_Eobs_corr3_cer[volName] -> Fill( (*m_Eobs)[volName] / f_EobsMaxTime_vs_NCerenCalib[volName][50]->Eval(x) );
        }
        
        
        //---------------------
        // time dependent plots      
        
        for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
        {
          G4String volName = volumeList->at(volIt);
          
          for(unsigned int timeTypeIt = 0; timeTypeIt < timeTypes.size(); ++timeTypeIt)
          {
            G4String timeType = timeTypes.at(timeTypeIt);
            G4String timeSliceName = timeType + "Low";
            
            
            // sum Eobs
            float EobsMaxTime = 0.;
            std::map<float,float> EobsIntTime;
            for(int bin = 1; bin <= nTimeSlices["Low"]+1; ++bin)
            {
              EobsMaxTime += (*m_Edep_byTime)[volName][timeSliceName][bin];
              
              for(unsigned int intTimeIt = 0; intTimeIt < intTimeVals.size(); ++intTimeIt)
              {
                G4float intTime = intTimeVals.at(intTimeIt);
                
                if( minTimes["Low"]+(bin-1)*timeSliceSizes["Low"] < intTime )
                  EobsIntTime[intTime] += (*m_Edep_byTime)[volName][timeSliceName][bin];
              }
            }
            
            
            // fill histo
            for(unsigned int intTimeIt = 0; intTimeIt < intTimeVals.size(); ++intTimeIt)
            {
              G4float intTime = intTimeVals.at(intTimeIt);
              float x = EobsIntTime[intTime]/EobsMaxTime;
              
              if( f_EobsMaxTime_vs_EobsIntTime[volName][timeSliceName][Ein][intTime] != NULL )
                h_Eobs_corr2_time[volName][timeSliceName][intTime] -> Fill( EobsMaxTime / f_EobsMaxTime_vs_EobsIntTime[volName][timeSliceName][Ein][intTime]->Eval(x) );
              
              if( f_EobsMaxTime_vs_EobsIntTime[volName][timeSliceName][50][intTime] != NULL )
                h_Eobs_corr3_time[volName][timeSliceName][intTime] -> Fill( EobsMaxTime / f_EobsMaxTime_vs_EobsIntTime[volName][timeSliceName][50][intTime]->Eval(x) );
            }
          }
        }
        
      } // end loop over events
      std::cout << "\n>>> loop over events done" << std::endl;
    }
    
    
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
