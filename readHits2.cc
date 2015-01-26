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
  
  std::map<G4String,TProfile2D*> Eobs_time_vs_z; // map<timeSliceName,TProfile2D*>
  std::map<G4String,TProfile2D*> Eobs_time_vs_R; // map<timeSliceName,TProfile2D*>
  std::map<G4String,TProfile2D*> EobsFrac_time_vs_z; // map<timeSliceName,TProfile2D*>
  std::map<G4String,TProfile2D*> EobsFrac_time_vs_R; // map<timeSliceName,TProfile2D*>
  
  G4float RunTotDepEnergy = 0.0;
  G4float RunTotObsEnergy = 0.0;
  G4float RunTotNCeren = 0.0;
  
  
  // initialize ROOT
  TSystem ts;
  gSystem->Load("libCintex");
  ROOT::Cintex::Cintex::Enable();
  gSystem->Load("libClassesDict");
  //ROOT::Cintex::Cintex::SetDebug(2);
  
  
  TFile* outfile = new TFile(argv[2], "RECREATE");
  
  
  TChain* Tevt = new TChain("EventTree","EventTree");
  TChain* Trh = new TChain("RunTree","RunTree");
  FillChain(Tevt,argv[1]);
  FillChain(Trh,argv[1]);
  
  Event* event = new Event();
  //TBranch* b_event = Tevt -> GetBranch("Event");
  //b_event -> SetAddress(&event);
  Tevt -> SetBranchAddress("Event",&event);
  
  RunHeader* runHeader = new RunHeader();
  TBranch* b_runHeader = Trh -> GetBranch("RunHeader");
  b_runHeader -> SetAddress(&runHeader);
  
  
  
  Trh->GetEntry(0);
  runHeader->Print();
  G4float Ein = runHeader->GetParticleEnergy();
  G4float Edep;
  G4float Eobs;
  G4float NCeren;
  std::vector<G4String>* volumeList = runHeader -> GetVolumes();
  std::vector<G4String>* particleList = runHeader -> GetParticleList();
  std::vector<G4String>* processList  = runHeader -> GetProcessList();
  
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
  
  
  // compute histogram axis
  std::vector<G4String>* Solids = runHeader->GetSolids();
  int numXCells  = runHeader->GetXCellNum();
  int numYCells  = runHeader->GetYCellNum();
  int numZLayers = runHeader->GetZLayerNum();
  int nBinsX = numXCells;
  int nBinsY = numYCells;
  int nBinsZ = numZLayers*Solids->size();
  std::cout << "nBinsX: " << nBinsX << "   nBinsY: " << nBinsY << "   nBinsZ: " << nBinsZ << std::endl;
  double* cellXLength = new double[Solids->size()];
  double* cellYLength = new double[Solids->size()];
  double* cellZLength = new double[Solids->size()];
  double layerZLength = 0.;  
  
  for(unsigned int solidIt = 0; solidIt < Solids->size(); ++solidIt)
  {
    std::cout << "Solid: " << Solids->at(solidIt) << std::endl;
    cellXLength[solidIt] = 2. * (runHeader->GetSolidsXHalfLength())->at(solidIt);
    cellYLength[solidIt] = 2. * (runHeader->GetSolidsYHalfLength())->at(solidIt);
    cellZLength[solidIt] = 2. * (runHeader->GetSolidsZHalfLength())->at(solidIt);
    layerZLength += cellZLength[solidIt];

    if( solidIt > 0 && ( (cellXLength[solidIt] != cellXLength[0]) ||
                         (cellYLength[solidIt] != cellYLength[0]) ) )
    {
      std::cout << ">>> Detector geometry not supported. Terminating... <<<" << std::endl;
      exit(-1);
    }
  }
  std::cout << "cell x length: " << cellXLength[0] << "   cell y length: " << cellYLength[0] << "   layer z length: " << layerZLength << std::endl;
  
  double* xAxis = new double[numXCells+1];
  for(int binx = 0; binx < numXCells+1; ++binx)
  {
    xAxis[binx] = -0.5*numXCells*cellXLength[0] + binx*cellXLength[0];
  }
  double* yAxis = new double[numYCells+1];
  for(int biny = 0; biny < numYCells+1; ++biny)
  {
    yAxis[biny] = -0.5*numYCells*cellYLength[0] + biny*cellYLength[0];
  }
  double* zAxis = new double[numZLayers*Solids->size()+1];
  for(int binz = 0; binz < numZLayers; ++binz)
  {
    for(unsigned int solidIt = 0; solidIt < Solids->size(); ++solidIt)
    {
      zAxis[binz*Solids->size()+solidIt] = binz*layerZLength + solidIt*cellZLength[solidIt>0?solidIt-1:0];
    }
  }
  zAxis[numZLayers*Solids->size()] = numZLayers*layerZLength;
  
  
  //------------------
  // define histograms
  
  TDirectory* globalDir = outfile->mkdir("globalHistos");
  globalDir->cd();
  
  TH1F* h_Edep = new TH1F("h_Edep", "Total energy deposited", 100, 0., 1.1*Ein);
  TH1F* h_Eobs = new TH1F("h_Eobs", "Total observed (Birks suppressed) energy", 100, 0., 1.1*Ein);
  TH1F* h_NCeren = new TH1F("h_NCeren", "Total nr. of cerenkov photons", 100, 0., 65000*1.1*Ein);
  
  
  
  TDirectory* timeDir = outfile->mkdir("timeHistos");
  timeDir->cd();
  
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
      double* Eaxis = new double[nBins+1];
      for(int bin = 0; bin <= nBins; ++bin) Eaxis[bin] = min+(max-min)/nBins*bin;
      
      Eobs_time_vs_z[timeSliceName] = new TProfile2D(Form("Eobs_%s_vs_z",timeSliceName.c_str()),"Eobs",nBinsZ,zAxis,nBins,Eaxis);
      Eobs_time_vs_R[timeSliceName] = new TProfile2D(Form("Eobs_%s_vs_R",timeSliceName.c_str()),"Eobs",nBinsZ,zAxis,nBins,Eaxis);
      EobsFrac_time_vs_z[timeSliceName] = new TProfile2D(Form("EobsFrac_%s_vs_z",timeSliceName.c_str()),"Eobs",nBinsZ,zAxis,nBins,Eaxis);
      EobsFrac_time_vs_R[timeSliceName] = new TProfile2D(Form("EobsFrac_%s_vs_R",timeSliceName.c_str()),"Eobs",nBinsZ,zAxis,nBins,Eaxis);
    }
  }
  
  
  
  //-----------------
  // loop over events
  
  G4int nEvents = Tevt->GetEntries();
  //nEvents = 10;
  G4cout << " Nr. of Events:  " << nEvents << " in input file "<< G4endl;
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
    h_Edep->Fill(Edep);
    h_Eobs->Fill(Eobs);
    h_NCeren->Fill(NCeren);
    
    
    std::map<G4String,std::map<G4String,std::map<G4ThreeVector,std::vector<G4VHit*> > > >* HCMap = event->GetHCMap();
    
    for(unsigned int volIt = 0; volIt < volumeList->size()-1; ++volIt)
    {
      std::string volName = volumeList->at(volIt);
      
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
              
              Eobs_time_vs_z[timeSliceName] -> Fill( z,time,aHit->GetEdep() );
              Eobs_time_vs_R[timeSliceName] -> Fill( sqrt(x*x+y*y+z*z),time,aHit->GetEobsbirks() );
              EobsFrac_time_vs_z[timeSliceName] -> Fill( z,time,aHit->GetEdep()/Eobs );
              EobsFrac_time_vs_R[timeSliceName] -> Fill( sqrt(x*x+y*y+z*z),time,aHit->GetEobsbirks()/Eobs );
            }
          }
        }
      }
    }
  }
  
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



TH1F* GetCumulative(TProfile* prof, const int& normalize, const float& Eobs)
{
  int nBinsX  = prof -> GetNbinsX();
  double xMin = prof -> GetBinLowEdge(1);
  double xMax = prof -> GetBinLowEdge(nBinsX) + prof -> GetBinWidth(nBinsX);
  
  TH1F* histo_cumul;
  if     ( normalize == 0 ) histo_cumul = new TH1F(Form("cumul%s",prof->GetName()),     "",nBinsX,xMin,xMax);
  else if( normalize == 1 ) histo_cumul = new TH1F(Form("normCumul%s",prof->GetName()), "",nBinsX,xMin,xMax);
  else if( normalize == 2 ) histo_cumul = new TH1F(Form("normECumul%s",prof->GetName()),"",nBinsX,xMin,xMax);
  
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
