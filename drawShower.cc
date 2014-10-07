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
#include "TSystem.h"
#include "TKey.h"
#include "TH1.h"
#include "TH2.h"
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
#include "DRTSCalorimeterHit2.hh"
#include "PhotonHit.hh"
#include "TrackerHit.hh"
//
//STL:
#include <vector>
#include<iomanip>
using namespace std;



int main(int argc, char** argv)
{
  //
  // initialize ROOT
  TSystem ts;
  gSystem->Load("libCintex");
  ROOT::Cintex::Cintex::Enable();
  gSystem->Load("libClassesDict");
  
  if (argc < 3) {
    G4cout << "Program requires 2 arguments: name of input file, name of output file" << G4endl;
    exit(1);
  }
  
  TFile* outfile = new TFile(argv[2], "RECREATE");
  TFile fo(argv[1]);
  fo.GetListOfKeys()->Print();
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
  double Ein = runh->GetParticleEnergy() / 1000.;
  double Edep;
  double Eobs;
  double RunTotEnergy = 0.0;
  double RunTotObsEnergy = 0.0;
  double RunTotNCeren = 0.0;  
  
  std::vector<G4String> volumeList = runh -> GetVolumes();
  std::vector<G4String> particleList = runh -> GetParticleList();
  std::vector<G4String> particleListCeren = runh -> GetParticleListCeren();
  std::vector<G4String> particleListShort = runh -> GetParticleList();
  
  G4double timeSliceLow = runh->GetTimeSliceLow();
  G4double minTimeLow   = runh->GetMinTimeLow();
  G4double maxTimeLow   = runh->GetMaxTimeLow();
  G4int nSlicesLow = int((maxTimeLow-minTimeLow)/timeSliceLow);
  G4double timeSliceMed = runh->GetTimeSliceMed();
  G4double minTimeMed   = runh->GetMinTimeMed();
  G4double maxTimeMed   = runh->GetMaxTimeMed();
  G4int nSlicesMed = int((maxTimeMed-minTimeMed)/timeSliceMed);
  G4double timeSliceHig = runh->GetTimeSliceHig();
  G4double minTimeHig   = runh->GetMinTimeHig();
  G4double maxTimeHig   = runh->GetMaxTimeHig();
  G4int nSlicesHig = int((maxTimeHig-minTimeHig)/timeSliceHig);  
  
  
  
  //-----------------
  // loop over events
  
  Int_t nevent = fevtbranch->GetEntries();
  G4cout << " Nr. of Events:  " << nevent << " in input file "<< G4endl;
  for (Int_t i = 0; i < nevent; i++)
  {
    G4cout << ">>> processing event " << i << " / " << nevent << G4endl;
    fevtbranch->GetEntry(i);
    
    
    outfile -> cd();
    
    std::map<G4int,TH2F*> h2_xy_low;
    std::map<G4int,TH2F*> h2_xy_med;
    std::map<G4int,TH2F*> h2_xy_hig;
    std::map<G4int,TH2F*> h2_zy_low;
    std::map<G4int,TH2F*> h2_zy_med;
    std::map<G4int,TH2F*> h2_zy_hig;
    
    std::map<G4int, std::map<G4String, TH2F*> > h2part_xy_low;
    std::map<G4int, std::map<G4String, TH2F*> > h2part_xy_med;
    std::map<G4int, std::map<G4String, TH2F*> > h2part_xy_hig;
    std::map<G4int, std::map<G4String, TH2F*> > h2part_zy_low;
    std::map<G4int, std::map<G4String, TH2F*> > h2part_zy_med;
    std::map<G4int, std::map<G4String, TH2F*> > h2part_zy_hig;
    
    for(G4int slice = 1; slice <= nSlicesLow; ++slice)
    {
      h2_xy_low[slice] = new TH2F(Form("h2_xy_globalSliceLow%d_evt%d",slice,i),"",40,-500.,500.,40.,-500.,500.);
      h2_zy_low[slice] = new TH2F(Form("h2_zy_globalSliceLow%d_evt%d",slice,i),"",200,0.,3000.,40.,-500.,500.);
      
      for(unsigned int partIt = 0; partIt < particleListShort.size(); ++partIt)
      {
        G4String partName = particleListShort.at(partIt);
        
        (h2part_xy_low[slice])[partName] = new TH2F(Form("h2%s_xy_globalSliceLow%d_evt%d",partName.c_str(),slice,i),"",40,-500.,500.,40.,-500.,500.);
        (h2part_zy_low[slice])[partName] = new TH2F(Form("h2%s_zy_globalSliceLow%d_evt%d",partName.c_str(),slice,i),"",200,0.,3000.,40.,-500.,500.);
      }
    }
    for(G4int slice = 1; slice <= nSlicesMed; ++slice)
    {
      h2_xy_med[slice] = new TH2F(Form("h2_xy_globalSliceMed%d_evt%d",slice,i),"",40,-500.,500.,40.,-500.,500.);
      h2_zy_med[slice] = new TH2F(Form("h2_zy_globalSliceMed%d_evt%d",slice,i),"",200,0.,3000.,40.,-500.,500.);
      
      for(unsigned int partIt = 0; partIt < particleListShort.size(); ++partIt)
      {
        G4String partName = particleListShort.at(partIt);
        
        (h2part_xy_med[slice])[partName] = new TH2F(Form("h2%s_xy_globalSliceMed%d_evt%d",partName.c_str(),slice,i),"",40,-500.,500.,40.,-500.,500.);
        (h2part_zy_med[slice])[partName] = new TH2F(Form("h2%s_zy_globalSliceMed%d_evt%d",partName.c_str(),slice,i),"",200,0.,3000.,40.,-500.,500.);
      }
    }
    for(G4int slice = 1; slice <= nSlicesHig; ++slice)
    {
      h2_xy_hig[slice] = new TH2F(Form("h2_xy_globalSliceHig%d_evt%d",slice,i),"",40,-500.,500.,40.,-500.,500.);
      h2_zy_hig[slice] = new TH2F(Form("h2_zy_globalSliceHig%d_evt%d",slice,i),"",200,0.,3000.,40.,-500.,500.);
      
      for(unsigned int partIt = 0; partIt < particleListShort.size(); ++partIt)
      {
        G4String partName = particleListShort.at(partIt);
        
        (h2part_xy_hig[slice])[partName] = new TH2F(Form("h2%s_xy_globalSliceHig%d_evt%d",partName.c_str(),slice,i),"",40,-500.,500.,40.,-500.,500.);
        (h2part_zy_hig[slice])[partName] = new TH2F(Form("h2%s_zy_globalSliceHig%d_evt%d",partName.c_str(),slice,i),"",200,0.,3000.,40.,-500.,500.);
      }
    }    
    
    double tot_en = 0;
    double obs_en = 0;
    int tot_nceren = 0;
    std::map<G4String, std::vector<G4VHit*> >* hcMap  = event->GetHCMap();
    std::map<G4String, std::vector<G4VHit*> >::const_iterator mapIt;
    for(mapIt = hcMap->begin(); mapIt != hcMap->end(); ++mapIt)
    {
      G4cout << mapIt -> first << " " << (mapIt->second).size() << G4endl;
    }
    for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
    {
      G4String volName = volumeList.at(volIt) + "_HC";
      std::vector<G4VHit*> hcVec  = (*hcMap)[volName];
      
      volName = volumeList.at(volIt) + "_HC2";
      std::vector<G4VHit*> hcVec2 = (*hcMap)[volName];
      
      for(unsigned int vecIt = 0; vecIt < hcVec.size(); ++vecIt)
      {
        DRTSCalorimeterHit* aHit = dynamic_cast<DRTSCalorimeterHit*>(hcVec.at(vecIt));
        // aHit -> Print();
      }
      
      for(unsigned int vecIt = 0; vecIt < hcVec2.size(); ++vecIt)
      {
        DRTSCalorimeterHit2* aHit = dynamic_cast<DRTSCalorimeterHit2*>(hcVec2.at(vecIt));
        G4ThreeVector pos = aHit -> GetPos();
        G4int globalSliceLow = aHit -> GetGlobalSliceLow();
        G4int globalSliceMed = aHit -> GetGlobalSliceMed();
        G4int globalSliceHig = aHit -> GetGlobalSliceHig();
        G4String partName = aHit -> GetParticleName();
        Edep = aHit -> GetEdep();

        for( int iSlice = nSlicesLow; iSlice >= globalSliceLow; --iSlice)
        {
          if( globalSliceLow >= 1 && globalSliceLow <= nSlicesLow )
          {
            h2_xy_low[iSlice] -> Fill(pos.x(),pos.y(),Edep);
            h2_zy_low[iSlice] -> Fill(pos.z(),pos.y(),Edep);
            (h2part_xy_low[iSlice])[partName] -> Fill(pos.x(),pos.y(),Edep);
            (h2part_zy_low[iSlice])[partName] -> Fill(pos.z(),pos.y(),Edep);
          }
        }
        
        for( int iSlice = nSlicesMed; iSlice >= globalSliceMed; --iSlice)
        {
          if( globalSliceMed >= 1 && globalSliceMed <= nSlicesMed )
          {
            h2_xy_med[iSlice] -> Fill(pos.x(),pos.y(),Edep);
            h2_zy_med[iSlice] -> Fill(pos.z(),pos.y(),Edep);
            (h2part_xy_med[iSlice])[partName] -> Fill(pos.x(),pos.y(),Edep);
            (h2part_zy_med[iSlice])[partName] -> Fill(pos.z(),pos.y(),Edep);
          }
        }
        
        for( int iSlice = nSlicesHig; iSlice >= globalSliceHig; --iSlice)
        {
          if( globalSliceHig >= 1 && globalSliceHig <= nSlicesHig )
          {
            h2_xy_hig[iSlice] -> Fill(pos.x(),pos.y(),Edep);
            h2_zy_hig[iSlice] -> Fill(pos.z(),pos.y(),Edep);
            (h2part_xy_hig[iSlice])[partName] -> Fill(pos.x(),pos.y(),Edep);
            (h2part_zy_hig[iSlice])[partName] -> Fill(pos.z(),pos.y(),Edep);
          }
        }
        
        //aHit -> Print();
      }
    }
    
    outfile -> cd();
    outfile -> mkdir(Form("event%d",i));
    outfile -> cd(Form("event%d",i));
    for(G4int slice = 1; slice <= nSlicesLow; ++slice)
    {
      h2_xy_low[slice] -> Write();
      h2_zy_low[slice] -> Write();
      
      delete h2_xy_low[slice];
      delete h2_zy_low[slice];
      
      for(unsigned int partIt = 0; partIt < particleListShort.size(); ++partIt)
      {
        G4String partName = particleListShort.at(partIt);
        
        (h2part_xy_low[slice])[partName] -> Write();
        (h2part_zy_low[slice])[partName] -> Write();
        
        delete (h2part_xy_low[slice])[partName];
        delete (h2part_zy_low[slice])[partName];
      }
    }
    for(G4int slice = 1; slice <= nSlicesMed; ++slice)
    {
      h2_xy_med[slice] -> Write();
      h2_zy_med[slice] -> Write();
      
      delete h2_xy_med[slice];
      delete h2_zy_med[slice];
      
      for(unsigned int partIt = 0; partIt < particleListShort.size(); ++partIt)
      {
        G4String partName = particleListShort.at(partIt);
        
        (h2part_xy_med[slice])[partName] -> Write();
        (h2part_zy_med[slice])[partName] -> Write();
        
        delete (h2part_xy_med[slice])[partName];
        delete (h2part_zy_med[slice])[partName];
      }
    }
    for(G4int slice = 1; slice <= nSlicesHig; ++slice)
    {
      h2_xy_hig[slice] -> Write();
      h2_zy_hig[slice] -> Write();
      
      delete h2_xy_hig[slice];
      delete h2_zy_hig[slice];
      
      for(unsigned int partIt = 0; partIt < particleListShort.size(); ++partIt)
      {
        G4String partName = particleListShort.at(partIt);
        
        (h2part_xy_hig[slice])[partName] -> Write();
        (h2part_zy_hig[slice])[partName] -> Write();
        
        delete (h2part_xy_hig[slice])[partName];
        delete (h2part_zy_hig[slice])[partName];
      }
    }
    outfile -> cd();
    
    Edep = event->GetTotEnergy() / 1000.;
    Eobs = event->GetTotObsEnergy() / 1000.;
    RunTotEnergy = RunTotEnergy + Edep;
    RunTotObsEnergy = RunTotObsEnergy + Eobs;
    RunTotNCeren = RunTotNCeren + event->GetNCeren();
  } // end loop over events
  
  
  outfile -> cd();
  TH1F* h_timeSliceLow = new TH1F("h_timeSliceLow","",nSlicesLow,minTimeLow,maxTimeLow);
  TH1F* h_timeSliceMed = new TH1F("h_timeSliceMed","",nSlicesMed,minTimeMed,maxTimeMed);
  TH1F* h_timeSliceHig = new TH1F("h_timeSliceHig","",nSlicesHig,minTimeHig,maxTimeHig);
  h_timeSliceLow->Write();
  h_timeSliceMed->Write();
  h_timeSliceHig->Write();
  outfile -> Close();
}
