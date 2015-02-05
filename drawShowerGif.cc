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

#include "Classes.hh"
#include "Event.hh"
#include "RunHeader.hh"
#include "DRCalorimeterHit.hh"
#include "setTDRStyle.hh"

#include "Cintex/Cintex.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TLatex.h"

#include <cstdio>
#include <vector>



int main(int argc, char** argv)
{
  //----------------
  // initialize ROOT
  
  setTDRStyle();
  TSystem ts;
  gSystem->Load("libCintex");
  ROOT::Cintex::Cintex::Enable();
  gSystem->Load("libClassesDict");
  
  if( argc < 5 )
  {
    G4cout << "Program requires 4 arguments: name of input file, eventId, timeType (globalTime | localTime1), timeSliceType (Low | Med | Hig)" << G4endl;
    exit(1);
  }
  
  G4String inFileName = std::string(argv[1]);
  G4int eventId = atoi(argv[2]);
  std::string timeType = std::string(argv[3]);
  std::string timeSliceType = std::string(argv[4]);
  std::string timeSliceName = timeType + timeSliceType;
  
  
  
  //----------------
  // initialize tree
  
  TFile* inFile = TFile::Open(inFileName.c_str(),"READ");
  
  Event* event = new Event();
  TTree* Tevt = (TTree*)(inFile->Get("EventTree"));
  Tevt -> SetBranchAddress("Event",&event);
  RunHeader* runHeader = new RunHeader();
  TTree* Trh = (TTree*)(inFile->Get("RunTree"));
  Trh -> SetBranchAddress("RunHeader",&runHeader);
  Trh -> GetEntry(0);
  runHeader->Print();
  
  // initialize variables
  G4float Ein = runHeader-> GetParticleEnergy();
  G4float minZaxis = 0.1/1000.;
  G4float maxZaxis = Ein/10.;
  
  std::vector<G4String>* volumeList = runHeader -> GetVolumes();
  std::map<G4String,G4bool> volumeEnabled;
  for(unsigned int volIt = 0; volIt < volumeList->size(); ++volIt)
  {
    volumeEnabled[volumeList->at(volIt)] = true;
  }
  
  std::vector<G4String>* particleList = runHeader->GetParticleList();
  std::map<G4String,G4bool> particleEnabled;
  for(unsigned int particleIt = 0; particleIt < particleList->size(); ++particleIt)
  {
    particleEnabled[particleList->at(particleIt)] = true;
  }
  
  std::vector<G4String>* processList = runHeader->GetProcessList();
  std::map<G4String,G4bool> processEnabled;
  for(unsigned int processIt = 0; processIt < processList->size(); ++processIt)
  {
    processEnabled[processList->at(processIt)] = true;
  }
  
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
  nTimeSlices["Low"] = (maxTimes["Low"]-minTimes["Low"])/(timeSliceSizes["Low"]);
  nTimeSlices["Med"] = (maxTimes["Med"]-minTimes["Med"])/(timeSliceSizes["Med"]);
  nTimeSlices["Hig"] = (maxTimes["Hig"]-minTimes["Hig"])/(timeSliceSizes["Hig"]);
  
  // compute histogram axis
  std::vector<G4String>* Solids = runHeader->GetSolids();
  std::vector<G4float>* SolidsXHalfLength = runHeader->GetSolidsXHalfLength();
  std::vector<G4float>* SolidsYHalfLength = runHeader->GetSolidsYHalfLength();
  std::vector<G4float>* SolidsZHalfLength = runHeader->GetSolidsZHalfLength();
  std::reverse(Solids->begin(),Solids->end());
  std::reverse(SolidsXHalfLength->begin(),SolidsXHalfLength->end());
  std::reverse(SolidsYHalfLength->begin(),SolidsYHalfLength->end());
  std::reverse(SolidsZHalfLength->begin(),SolidsZHalfLength->end());
  
  int numXCells  = runHeader->GetXCellNum();
  int numYCells  = runHeader->GetYCellNum();
  int numZLayers = runHeader->GetZLayerNum();
  G4int nBinsX = numXCells;
  G4int nBinsY = numYCells;
  G4int nBinsZ = numZLayers*Solids->size();
  std::cout << "nBinsX: " << nBinsX << "   nBinsY: " << nBinsY << "   nBinsZ: " << nBinsZ << std::endl;
  double* cellXLength = new double[Solids->size()];
  double* cellYLength = new double[Solids->size()];
  double* cellZLength = new double[Solids->size()];
  double layerZLength = 0.;
  
  for(unsigned int solidIt = 0; solidIt < Solids->size(); ++solidIt)
  {
    std::cout << "Solid: " << Solids->at(solidIt) << std::endl;
    cellXLength[solidIt] = 2. * SolidsXHalfLength->at(solidIt);
    cellYLength[solidIt] = 2. * SolidsYHalfLength->at(solidIt);
    cellZLength[solidIt] = 2. * SolidsZHalfLength->at(solidIt);
    layerZLength += cellZLength[solidIt];
    std::cout << "Solid: " << Solids->at(solidIt) << "   z: " << cellZLength[solidIt] << std::endl;
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
  
  
  // Get a particular event
  Tevt -> GetEntry(eventId);
  std::cout << ">>> reading event " << eventId << " / " << Tevt->GetEntriesFast() << std::endl;

  std::map<G4String,std::map<G4String,std::map<G4ThreeVector,std::vector<G4VHit*> > > >* HCMap = event->GetHCMap();
  
  
  
  //----------------
  // Draw histograms
  
  TCanvas* c1 = new TCanvas("c1","c1",1200,600);
  c1->SetFillColor(0);
  c1->SetGrid();
  c1 -> Divide(2,1);
  
  c1 -> cd(1);
  gPad -> SetLogz();
  gPad->SetRightMargin(0.20);
  
  c1 -> cd(2);
  gPad -> SetLogz();
  gPad->SetRightMargin(0.20);
  
  
  // Fill histograms
  TH2F* h2_yx;
  TH2F* h2_yz;
  
  h2_yx = new TH2F(Form("h2_yx"),"y-x view",nBinsY,yAxis,nBinsX,xAxis);
  h2_yx -> SetTitle(";x (mm);y (mm);E_{obs} (GeV)");
  h2_yx -> GetZaxis() -> SetTitleOffset(1.50);
  h2_yx -> SetMinimum(minZaxis);
  h2_yx -> SetMaximum(maxZaxis);
  
  h2_yz = new TH2F(Form("h2_yz"),"y-z view",nBinsZ,zAxis,nBinsY,yAxis);
  h2_yz -> SetTitle(";z (mm);y (mm);E_{obs} (GeV)");
  h2_yz -> GetZaxis() -> SetTitleOffset(1.50); 
  h2_yz -> SetMinimum(minZaxis);
  h2_yz -> SetMaximum(maxZaxis);
  
  
  G4double EobsTot = 0.;
  TLatex* latex1 = NULL;
  TLatex* latex2 = NULL;
  for(int tsIt = 0; tsIt < nTimeSlices[timeSliceType]; ++tsIt)
  {
    //std::cout << ">>> processing time slice " << tsIt << " / " << nTimeSlices[timeSliceType] << "\r" << std::flush;
    std::cout << ">>> processing time slice " << tsIt << " / " << nTimeSlices[timeSliceType] << std::endl;
    
    for(unsigned int volIt = 0; volIt < volumeList->size(); ++volIt)
    {
      G4String volumeName = volumeList->at(volIt);
      if( volumeEnabled[volumeName] == false ) continue;
      
      std::map<G4ThreeVector,std::vector<G4VHit*> > posHitMap = (*HCMap)[volumeName][timeSliceName];
      for(std::map<G4ThreeVector,std::vector<G4VHit*> >::const_iterator mapIt = posHitMap.begin(); mapIt != posHitMap.end(); ++mapIt)
      {
        G4double x = (mapIt->first).x();
        G4double y = (mapIt->first).y();
        G4double z = (mapIt->first).z();
        std::vector<G4VHit*> hitVec = mapIt->second;
        G4double Eobs = 0.;
        for(unsigned int vecIt = 0; vecIt < hitVec.size(); ++vecIt)
        {
          DRTSCalorimeterHit2* aHit = dynamic_cast<DRTSCalorimeterHit2*>(hitVec.at(vecIt));
          G4int timeSlice = aHit -> GetTimeSlice();
          G4String particleName = aHit -> GetParticleName();
          G4String processName = aHit -> GetProcessName();
          if( particleEnabled[particleName] == false ) continue;
          if( processEnabled[processName] == false ) continue;
          if( timeSlice == tsIt )
          {
            Eobs += aHit -> GetEobsbirks(); 
            EobsTot += aHit -> GetEobsbirks();
            //std::cout << "(x,y,z) = (" << x << "," << y << "," << z << ")" << "   ts: " << timeSlice << "   Eobs: " << Eobs << "   EobsTot: " << EobsTot << std::endl;
          }
        }
        
        if( Eobs > 0 )
        {
          h2_yx -> Fill(x,y,Eobs);
          h2_yz -> Fill(z,y,Eobs);
        }
      }
    }
    
    if( latex1 != NULL ) delete latex1;
    latex1 = new TLatex(0.13,0.96,Form("%s: %.2f ns",timeType.c_str(),minTimes[timeSliceType]+tsIt*timeSliceSizes[timeSliceType]));
    latex1 -> SetNDC();
    latex1 -> SetTextFont(42);
    latex1 -> SetTextSize(0.04);
    
    if( latex2 != NULL ) delete latex2;
    latex2 = new TLatex(0.55,0.96,Form("E_{obs}: %.1f GeV",EobsTot));
    latex2 -> SetNDC();
    latex2 -> SetTextFont(42);
    latex2 -> SetTextSize(0.04);
    
    c1 -> cd(1);
    h2_yx -> Draw("COLZ");
    latex1 -> Draw("same");
    latex2 -> Draw("same");
    c1 -> cd(2);
    h2_yz -> Draw("COLZ");
    latex1 -> Draw("same");
    latex2 -> Draw("same");
    
    if( tsIt == 0 )
    {
      remove(Form("c_drawShowerGif_event%d_%s.gif",eventId,timeSliceName.c_str()));
    }
    
    //c1 -> Print(Form("c_drawShowerGif_event%d_%s%03d.png",eventId,timeSliceName.c_str(),tsIt));
    c1 -> Print(Form("c_drawShowerGif_event%d_%s.gif+10",eventId,timeSliceName.c_str()));
  }
  
  
  
  
  // TCanvas::Update() draws the frame, after which it can be changed
  c1->Update();
  c1->GetFrame()->SetFillColor(kWhite);
  c1->GetFrame()->SetBorderSize(12);
}
