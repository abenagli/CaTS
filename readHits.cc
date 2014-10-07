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
using namespace std;



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



int main(int argc, char** argv)
{
    std::map<G4String, TDirectory *> VolumeDirset;
    std::map<G4String, G4double> EvtEnergy; // Evt. energy by SD
    std::map<G4String, G4double> EvtObsEnergy; // Evt. observed energy by SD 
    std::map<G4String, G4double> EvtNceren; // Evt. cerenkov photons by SD
    std::map<G4String, G4double> RunEnergy; // Run energy by SD
    std::map<G4String, G4double> RunObsEnergy; // Run observed energy by SD 
    std::map<G4String, G4double> RunNceren; // Run cerenkov photons by SD
    std::map<G4String, std::map<G4String, G4double> > E_byParticle; // Run Sum: Energy deposited by particle type
    std::map<G4String, std::map<G4String, G4double> > Eobs_byParticle; // Run Sum: Observed energy by particle type
    std::map<G4String, std::map<G4String, G4double> > NCeren_byParticle; // Run Sum:Cerenkov photons by particle type
    std::map<G4String, std::map<G4String, std::map<G4String, products> > > processMap; // keep track of particles produced by a process
    std::map<G4String, std::map<G4String, G4int > > processMult; // Run Sum: Keeps track how often a specific process occurs
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndGlobalTimeLow; // Run Sum: Energy deposited by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndGlobalTimeMed; // Run Sum: Energy deposited by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndGlobalTimeHig; // Run Sum: Energy deposited by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndGlobalTimeLow; // Run Sum: Observed energy by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndGlobalTimeMed; // Run Sum: Observed energy by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndGlobalTimeHig; // Run Sum: Observed energy by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndGlobalTimeLow; // Run Sum:Cerenkov photons by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndGlobalTimeMed; // Run Sum:Cerenkov photons by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndGlobalTimeHig; // Run Sum:Cerenkov photons by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndLocalTime1Low; // Run Sum: Energy deposited by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndLocalTime1Med; // Run Sum: Energy deposited by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndLocalTime1Hig; // Run Sum: Energy deposited by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndLocalTime1Low; // Run Sum: Observed energy by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndLocalTime1Med; // Run Sum: Observed energy by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndLocalTime1Hig; // Run Sum: Observed energy by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndLocalTime1Low; // Run Sum:Cerenkov photons by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndLocalTime1Med; // Run Sum:Cerenkov photons by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndLocalTime1Hig; // Run Sum:Cerenkov photons by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndLocalTime2Low; // Run Sum: Energy deposited by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndLocalTime2Med; // Run Sum: Energy deposited by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndLocalTime2Hig; // Run Sum: Energy deposited by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndLocalTime2Low; // Run Sum: Observed energy by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndLocalTime2Med; // Run Sum: Observed energy by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndLocalTime2Hig; // Run Sum: Observed energy by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndLocalTime2Low; // Run Sum:Cerenkov photons by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndLocalTime2Med; // Run Sum:Cerenkov photons by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndLocalTime2Hig; // Run Sum:Cerenkov photons by particle type
    std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > processAndTimeLowMap; // keep track of particles produced by a process
    std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > processAndTimeMedMap; // keep track of particles produced by a process
    std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > processAndTimeHigMap; // keep track of particles produced by a process
    std::map<G4String, std::map<G4int, std::map<G4String, G4int > > > processAndTimeLowMult; // Run Sum: Keeps track how often a specific process occurs
    std::map<G4String, std::map<G4int, std::map<G4String, G4int > > > processAndTimeMedMult; // Run Sum: Keeps track how often a specific process occurs
    std::map<G4String, std::map<G4int, std::map<G4String, G4int > > > processAndTimeHigMult; // Run Sum: Keeps track how often a specific process occurs
    //
    std::map<G4String, std::map<G4String, TH1F*> > EHistosbySD; // Energy deposited in various SD's
    std::map<G4String, std::map<G4String, TH1F*> > EObsHistosbySD; // observed Energy deposited in various SD's
    std::map<G4String, std::map<G4String, TH1F*> > NCHistosbySD; // Number of Cerenkov photons in various SD's
    //
    std::map<G4String, std::map<G4String, TH1F*> > EvtEnergybySDandpart; // <SDname, std::map<particlename, histo of Energy contrib> >
    std::map<G4String, std::map<G4String, TH1F*> > EvtObsEnergybySDandpart; // <SDname, std::map<particlename, histo of observed Energy contrib> >
    std::map<G4String, std::map<G4String, TH1F*> > EvtNcerenbySDandpart; // <SDname, std::map<particlename, histo of Number of Cerenkov photons> >
    
    TProfile* EvtEnergybyglobaltimelow;
    TProfile* EvtEnergybyglobaltimemed;
    TProfile* EvtEnergybyglobaltimehig;
    TProfile* EvtEnergybylocaltime1low;
    TProfile* EvtEnergybylocaltime1med;
    TProfile* EvtEnergybylocaltime1hig;
    TProfile* EvtEnergybylocaltime2low;
    TProfile* EvtEnergybylocaltime2med;
    TProfile* EvtEnergybylocaltime2hig;
      
    std::map<G4String, TProfile*> EvtEnergybypartandglobaltimelow;    // std::map<particlename, histo of Energy contrib>
    std::map<G4String, TProfile*> EvtObsEnergybypartandglobaltimelow; // std::map<particlename, histo of observed Energy contrib>
    std::map<G4String, TProfile*> EvtNcerenbypartandglobaltimelow;    // std::map<particlename, histo of Number of Cerenkov photons>
    std::map<G4String, TProfile*> EvtEnergybypartandglobaltimemed;    // std::map<particlename, histo of Energy contrib>
    std::map<G4String, TProfile*> EvtObsEnergybypartandglobaltimemed; // std::map<particlename, histo of observed Energy contrib>
    std::map<G4String, TProfile*> EvtNcerenbypartandglobaltimemed;    // std::map<particlename, histo of Number of Cerenkov photons>
    std::map<G4String, TProfile*> EvtEnergybypartandglobaltimehig;    // std::map<particlename, histo of Energy contrib>
    std::map<G4String, TProfile*> EvtObsEnergybypartandglobaltimehig; // std::map<particlename, histo of observed Energy contrib>
    std::map<G4String, TProfile*> EvtNcerenbypartandglobaltimehig;    // std::map<particlename, histo of Number of Cerenkov photons>
    std::map<G4String, TProfile*> EvtEnergybypartandlocaltime1low;    // std::map<particlename, histo of Energy contrib>
    std::map<G4String, TProfile*> EvtObsEnergybypartandlocaltime1low; // std::map<particlename, histo of observed Energy contrib>
    std::map<G4String, TProfile*> EvtNcerenbypartandlocaltime1low;    // std::map<particlename, histo of Number of Cerenkov photons>
    std::map<G4String, TProfile*> EvtEnergybypartandlocaltime1med;    // std::map<particlename, histo of Energy contrib>
    std::map<G4String, TProfile*> EvtObsEnergybypartandlocaltime1med; // std::map<particlename, histo of observed Energy contrib>
    std::map<G4String, TProfile*> EvtNcerenbypartandlocaltime1med;    // std::map<particlename, histo of Number of Cerenkov photons>
    std::map<G4String, TProfile*> EvtEnergybypartandlocaltime1hig;    // std::map<particlename, histo of Energy contrib>
    std::map<G4String, TProfile*> EvtObsEnergybypartandlocaltime1hig; // std::map<particlename, histo of observed Energy contrib>
    std::map<G4String, TProfile*> EvtNcerenbypartandlocaltime1hig;    // std::map<particlename, histo of Number of Cerenkov photons>
    std::map<G4String, TProfile*> EvtEnergybypartandlocaltime2low;    // std::map<particlename, histo of Energy contrib>
    std::map<G4String, TProfile*> EvtObsEnergybypartandlocaltime2low; // std::map<particlename, histo of observed Energy contrib>
    std::map<G4String, TProfile*> EvtNcerenbypartandlocaltime2low;    // std::map<particlename, histo of Number of Cerenkov photons>
    std::map<G4String, TProfile*> EvtEnergybypartandlocaltime2med;    // std::map<particlename, histo of Energy contrib>
    std::map<G4String, TProfile*> EvtObsEnergybypartandlocaltime2med; // std::map<particlename, histo of observed Energy contrib>
    std::map<G4String, TProfile*> EvtNcerenbypartandlocaltime2med;    // std::map<particlename, histo of Number of Cerenkov photons>
    std::map<G4String, TProfile*> EvtEnergybypartandlocaltime2hig;    // std::map<particlename, histo of Energy contrib>
    std::map<G4String, TProfile*> EvtObsEnergybypartandlocaltime2hig; // std::map<particlename, histo of observed Energy contrib>
    std::map<G4String, TProfile*> EvtNcerenbypartandlocaltime2hig;    // std::map<particlename, histo of Number of Cerenkov photons>
    
    std::map<G4String, std::map<G4String, TProfile*> > EvtEnergybySDandpartandglobaltimelow;    // <SDname, std::map<particlename, histo of Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtObsEnergybySDandpartandglobaltimelow; // <SDname, std::map<particlename, histo of observed Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtNcerenbySDandpartandglobaltimelow;    // <SDname, std::map<particlename, histo of Number of Cerenkov photons> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtEnergybySDandpartandglobaltimemed;    // <SDname, std::map<particlename, histo of Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtObsEnergybySDandpartandglobaltimemed; // <SDname, std::map<particlename, histo of observed Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtNcerenbySDandpartandglobaltimemed;    // <SDname, std::map<particlename, histo of Number of Cerenkov photons> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtEnergybySDandpartandglobaltimehig;    // <SDname, std::map<particlename, histo of Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtObsEnergybySDandpartandglobaltimehig; // <SDname, std::map<particlename, histo of observed Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtNcerenbySDandpartandglobaltimehig;    // <SDname, std::map<particlename, histo of Number of Cerenkov photons> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtEnergybySDandpartandlocaltime1low;    // <SDname, std::map<particlename, histo of Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtObsEnergybySDandpartandlocaltime1low; // <SDname, std::map<particlename, histo of observed Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtNcerenbySDandpartandlocaltime1low;    // <SDname, std::map<particlename, histo of Number of Cerenkov photons> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtEnergybySDandpartandlocaltime1med;    // <SDname, std::map<particlename, histo of Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtObsEnergybySDandpartandlocaltime1med; // <SDname, std::map<particlename, histo of observed Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtNcerenbySDandpartandlocaltime1med;    // <SDname, std::map<particlename, histo of Number of Cerenkov photons> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtEnergybySDandpartandlocaltime1hig;    // <SDname, std::map<particlename, histo of Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtObsEnergybySDandpartandlocaltime1hig; // <SDname, std::map<particlename, histo of observed Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtNcerenbySDandpartandlocaltime1hig;    // <SDname, std::map<particlename, histo of Number of Cerenkov photons> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtEnergybySDandpartandlocaltime2low;    // <SDname, std::map<particlename, histo of Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtObsEnergybySDandpartandlocaltime2low; // <SDname, std::map<particlename, histo of observed Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtNcerenbySDandpartandlocaltime2low;    // <SDname, std::map<particlename, histo of Number of Cerenkov photons> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtEnergybySDandpartandlocaltime2med;    // <SDname, std::map<particlename, histo of Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtObsEnergybySDandpartandlocaltime2med; // <SDname, std::map<particlename, histo of observed Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtNcerenbySDandpartandlocaltime2med;    // <SDname, std::map<particlename, histo of Number of Cerenkov photons> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtEnergybySDandpartandlocaltime2hig;    // <SDname, std::map<particlename, histo of Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtObsEnergybySDandpartandlocaltime2hig; // <SDname, std::map<particlename, histo of observed Energy contrib> >
    std::map<G4String, std::map<G4String, TProfile*> > EvtNcerenbySDandpartandlocaltime2hig;    // <SDname, std::map<particlename, histo of Number of Cerenkov photons> >
    
    std::map<G4String,TProfile2D*> h2_xy_firstTime;
    std::map<G4String,TProfile2D*> h2_zy_firstTime;
    
    //
    double RunTotEnergy = 0.0;
    double RunTotObsEnergy = 0.0;
    E_byParticle.clear();
    Eobs_byParticle.clear();
    double RunTotNCeren = 0.0;
    //
    VolumeDirset.clear();
    //
    EvtEnergy.clear();
    EvtObsEnergy.clear();
    EvtNceren.clear();
    RunEnergy.clear();
    RunObsEnergy.clear();
    RunNceren.clear();
    //
    EHistosbySD.clear();
    EObsHistosbySD.clear();
    NCHistosbySD.clear();
    EvtEnergybySDandpart.clear();
    EvtObsEnergybySDandpart.clear();
    EvtNcerenbySDandpart.clear();
    // EvtEnergybySDandpartandglobaltimelow.clear();
    // EvtObsEnergybySDandpartandglobaltimelow.clear();
    // EvtNcerenbySDandpartandglobaltimelow.clear();
    // EvtEnergybySDandpartandglobaltimemed.clear();
    // EvtObsEnergybySDandpartandglobaltimemed.clear();
    // EvtNcerenbySDandpartandglobaltimemed.clear();
    // EvtEnergybySDandpartandglobaltimehig.clear();
    // EvtObsEnergybySDandpartandglobaltimehig.clear();
    // EvtNcerenbySDandpartandglobaltimehig.clear();
    // EvtEnergybySDandpartandlocaltime1low.clear();
    // EvtObsEnergybySDandpartandlocaltime1low.clear();
    // EvtNcerenbySDandpartandlocaltime1low.clear();
    // EvtEnergybySDandpartandlocaltime1med.clear();
    // EvtObsEnergybySDandpartandlocaltime1med.clear();
    // EvtNcerenbySDandpartandlocaltime1med.clear();
    // EvtEnergybySDandpartandlocaltime1hig.clear();
    // EvtObsEnergybySDandpartandlocaltime1hig.clear();
    // EvtNcerenbySDandpartandlocaltime1hig.clear();
    // EvtEnergybySDandpartandlocaltime2low.clear();
    // EvtObsEnergybySDandpartandlocaltime2low.clear();
    // EvtNcerenbySDandpartandlocaltime2low.clear();
    // EvtEnergybySDandpartandlocaltime2med.clear();
    // EvtObsEnergybySDandpartandlocaltime2med.clear();
    // EvtNcerenbySDandpartandlocaltime2med.clear();
    // EvtEnergybySDandpartandlocaltime2hig.clear();
    // EvtObsEnergybySDandpartandlocaltime2hig.clear();
    // EvtNcerenbySDandpartandlocaltime2hig.clear();
    //
    // initialize ROOT
    TSystem ts;
    gSystem->Load("libCintex");
    ROOT::Cintex::Cintex::Enable();
    gSystem->Load("libClassesDict");
    //  ROOT::Cintex::Cintex::SetDebug(2);

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
    std::vector<G4String> volumeList = runh -> GetVolumes();
    std::vector<G4String> particleList = runh -> GetParticleList();
    std::vector<G4String> particleListCeren = runh -> GetParticleListCeren();
    
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
    
    
    //------------------
    // define histograms
    
    TDirectory *histotop = outfile->mkdir("histos");
    histotop->cd();
    TH1F* Eobshisto = new TH1F("Eobshisto", "Total observed (Birks suppressed) energy", 100, 0.0, Ein + 0.2 * Ein);
    TH1F* Ehisto = new TH1F("Ehisto", "Total energy deposited", 100, 0.0, Ein + Ein);
    TH1F* Ehisto2 = new TH1F("Ehisto2", "Total energy deposited", 100, 0.0, Ein + 0.1 * Ein);
    TH1F* Ehisto3 = new TH1F("Ehisto3", "Total energy observed", 100, 0.0, Ein + 2 * Ein);
    TH1F* Nhisto = new TH1F("Nhisto", "Total nr of cerenkov photons", 100, 0.0, 65000 * (Ein + 0.1 * Ein));
    // TH1F* tedephisto = new TH1F("Edeptime", "time of energy deposit", 100, 0.0, 10000.);
    // TH1F* tedepEMhisto = new TH1F("EdepEMtime", "time of EM energy deposit", 100, 0.0, 10000.);
    // TH1F* tedepnEMhisto = new TH1F("EdepnEMtime", "time of non EM energy deposit", 100, 0.0, 10000.);
    // TH1F* tNChisto = new TH1F("NCtime", "time of Cerenkov photon creation", 100, 0.0, 10000.); 
    TProfile2D* h2_xy_firstTime_all = new TProfile2D(Form("h2_xy_firstTime"),"",40,-500.,500.,40.,-500.,500.);;
    TProfile2D* h2_zy_firstTime_all = new TProfile2D(Form("h2_zy_firstTime"),"",200,0.,3000.,40.,-500.,500.);;
    
    
    {
      std::string volName = "AllVol";
      
      // create the subdirectories for all the different SD (overkill): 
      if (VolumeDirset.find(volName) == VolumeDirset.end())
        VolumeDirset.insert(std::make_pair(volName,histotop->mkdir(volName.c_str())));
      
      VolumeDirset[volName]->cd();
      
      h2_xy_firstTime[volName] = new TProfile2D(Form("h2_xy_firstTime"),"",40,-500.,500.,40.,-500.,500.);
      h2_zy_firstTime[volName] = new TProfile2D(Form("h2_zy_firstTime"),"",200,0.,3000.,40.,-500.,500.);
      
      
      EvtEnergybyglobaltimelow = new TProfile("Eobs_vs_globalTimeLow",  "",nSlicesLow+2,minTimeLow-timeSliceLow,maxTimeLow+timeSliceLow);
      EvtEnergybyglobaltimemed = new TProfile("Eobs_vs_globalTimeMed",  "",nSlicesMed+2,minTimeMed-timeSliceMed,maxTimeMed+timeSliceMed);
      EvtEnergybyglobaltimehig = new TProfile("Eobs_vs_globalTimeHig",  "",nSlicesHig+2,minTimeHig-timeSliceHig,maxTimeHig+timeSliceHig);
      
      EvtEnergybylocaltime1low = new TProfile("Eobs_vs_localTime1Low",  "",nSlicesLow+2,minTimeLow-timeSliceLow,maxTimeLow+timeSliceLow);
      EvtEnergybylocaltime1med = new TProfile("Eobs_vs_localTime1Med",  "",nSlicesMed+2,minTimeMed-timeSliceMed,maxTimeMed+timeSliceMed);
      EvtEnergybylocaltime1hig = new TProfile("Eobs_vs_localTime1Hig",  "",nSlicesHig+2,minTimeHig-timeSliceHig,maxTimeHig+timeSliceHig);
      
      EvtEnergybylocaltime2low = new TProfile("Eobs_vs_localTime2Low",  "",nSlicesLow+2,minTimeLow-timeSliceLow,maxTimeLow+timeSliceLow);
      EvtEnergybylocaltime2med = new TProfile("Eobs_vs_localTime2Med",  "",nSlicesMed+2,minTimeMed-timeSliceMed,maxTimeMed+timeSliceMed);
      EvtEnergybylocaltime2hig = new TProfile("Eobs_vs_localTime2Hig",  "",nSlicesHig+2,minTimeHig-timeSliceHig,maxTimeHig+timeSliceHig);
      
      
      for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
      {
        std::string partName = particleList.at(partIt);
        
        std::string histName = Form("Eobs_%s_vs_globalTimeLow",partName.c_str());
        EvtEnergybypartandglobaltimelow[partName] = new TProfile(histName.c_str(),"",nSlicesLow+2,minTimeLow-timeSliceLow,maxTimeLow+timeSliceLow);
        histName = Form("Eobs_%s_vs_globalTimeMed",partName.c_str());
        EvtEnergybypartandglobaltimemed[partName] = new TProfile(histName.c_str(),"",nSlicesMed+2,minTimeMed-timeSliceMed,maxTimeMed+timeSliceMed);
        histName = Form("Eobs_%s_vs_globalTimeHig",partName.c_str());
        EvtEnergybypartandglobaltimehig[partName] = new TProfile(histName.c_str(),"",nSlicesHig+2,minTimeHig-timeSliceHig,maxTimeHig+timeSliceHig);
        
        histName = Form("Eobs_%s_vs_localTime1Low",partName.c_str());
        EvtEnergybypartandlocaltime1low[partName] = new TProfile(histName.c_str(),"",nSlicesLow+2,minTimeLow-timeSliceLow,maxTimeLow+timeSliceLow);
        histName = Form("Eobs_%s_vs_localTime1Med",partName.c_str());
        EvtEnergybypartandlocaltime1med[partName] = new TProfile(histName.c_str(),"",nSlicesMed+2,minTimeMed-timeSliceMed,maxTimeMed+timeSliceMed);
        histName = Form("Eobs_%s_vs_localTime1Hig",partName.c_str());
        EvtEnergybypartandlocaltime1hig[partName] = new TProfile(histName.c_str(),"",nSlicesHig+2,minTimeHig-timeSliceHig,maxTimeHig+timeSliceHig);
        
        histName = Form("Eobs_%s_vs_localTime2Low",partName.c_str());
        EvtEnergybypartandlocaltime2low[partName] = new TProfile(histName.c_str(),"",nSlicesLow+2,minTimeLow-timeSliceLow,maxTimeLow+timeSliceLow);
        histName = Form("Eobs_%s_vs_localTime2Med",partName.c_str());
        EvtEnergybypartandlocaltime2med[partName] = new TProfile(histName.c_str(),"",nSlicesMed+2,minTimeMed-timeSliceMed,maxTimeMed+timeSliceMed);
        histName = Form("Eobs_%s_vs_localTime2Hig",partName.c_str());
        EvtEnergybypartandlocaltime2hig[partName] = new TProfile(histName.c_str(),"",nSlicesHig+2,minTimeHig-timeSliceHig,maxTimeHig+timeSliceHig);
      }
    }
    
    
    for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
    {
      std::string volName = volumeList.at(volIt);
      
      // create the subdirectories for all the different SD (overkill): 
      if (VolumeDirset.find(volName) == VolumeDirset.end())
        VolumeDirset.insert(std::make_pair(volName,histotop->mkdir(volName.c_str())));
      
      VolumeDirset[volName]->cd();
      
      h2_xy_firstTime[volName] = new TProfile2D(Form("h2_xy_firstTime"),"",40,-500.,500.,40.,-500.,500.);
      h2_zy_firstTime[volName] = new TProfile2D(Form("h2_zy_firstTime"),"",200,0.,3000.,40.,-500.,500.);
      
      for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
      {
        std::string partName = particleList.at(partIt);
        
        std::string histName = Form("Eobs_%s_vs_globalTimeLow",partName.c_str());
        (EvtEnergybySDandpartandglobaltimelow[volName])[partName] = new TProfile(histName.c_str(),"",nSlicesLow+2,minTimeLow-timeSliceLow,maxTimeLow+timeSliceLow);
        histName = Form("Eobs_%s_vs_globalTimeMed",partName.c_str());
        (EvtEnergybySDandpartandglobaltimemed[volName])[partName] = new TProfile(histName.c_str(),"",nSlicesMed+2,minTimeMed-timeSliceMed,maxTimeMed+timeSliceMed);
        histName = Form("Eobs_%s_vs_globalTimeHig",partName.c_str());
        (EvtEnergybySDandpartandglobaltimehig[volName])[partName] = new TProfile(histName.c_str(),"",nSlicesHig+2,minTimeHig-timeSliceHig,maxTimeHig+timeSliceHig);
        
        histName = Form("Eobs_%s_vs_localTime1Low",partName.c_str());
        (EvtEnergybySDandpartandlocaltime1low[volName])[partName] = new TProfile(histName.c_str(),"",nSlicesLow+2,minTimeLow-timeSliceLow,maxTimeLow+timeSliceLow);
        histName = Form("Eobs_%s_vs_localTime1Med",partName.c_str());
        (EvtEnergybySDandpartandlocaltime1med[volName])[partName] = new TProfile(histName.c_str(),"",nSlicesMed+2,minTimeMed-timeSliceMed,maxTimeMed+timeSliceMed);
        histName = Form("Eobs_%s_vs_localTime1Hig",partName.c_str());
        (EvtEnergybySDandpartandlocaltime1hig[volName])[partName] = new TProfile(histName.c_str(),"",nSlicesHig+2,minTimeHig-timeSliceHig,maxTimeHig+timeSliceHig);
        
        histName = Form("Eobs_%s_vs_localTime2Low",partName.c_str());
        (EvtEnergybySDandpartandlocaltime2low[volName])[partName] = new TProfile(histName.c_str(),"",nSlicesLow+2,minTimeLow-timeSliceLow,maxTimeLow+timeSliceLow);
        histName = Form("Eobs_%s_vs_localTime2Med",partName.c_str());
        (EvtEnergybySDandpartandlocaltime2med[volName])[partName] = new TProfile(histName.c_str(),"",nSlicesMed+2,minTimeMed-timeSliceMed,maxTimeMed+timeSliceMed);
        histName = Form("Eobs_%s_vs_localTime2Hig",partName.c_str());
        (EvtEnergybySDandpartandlocaltime2hig[volName])[partName] = new TProfile(histName.c_str(),"",nSlicesHig+2,minTimeHig-timeSliceHig,maxTimeHig+timeSliceHig);
      }
    }
    
    
    
    //-----------------
    // loop over events
    
    Int_t nevent = fevtbranch->GetEntries();
    G4cout << " Nr. of Events:  " << nevent << " in input file "<< G4endl;
    for (Int_t i = 0; i < nevent; i++)
    {
      std::cout << ">>> processing event " << i << " / " << nevent << "\r" << std::flush;
      fevtbranch->GetEntry(i);
      
      double tot_en = 0;
      double obs_en = 0;
      int tot_nceren = 0;
      std::map<G4String, std::vector<G4VHit*> >* hcMap = event->GetHCMap();
      Edep = event->GetTotEnergy() / 1000.;
      Eobs = event->GetTotObsEnergy() / 1000.;
      RunTotEnergy = RunTotEnergy + Edep;
      RunTotObsEnergy = RunTotObsEnergy + Eobs;
      Ehisto2->Fill(Edep);
      Ehisto3->Fill(Eobs);
      RunTotNCeren = RunTotNCeren + event->GetNCeren();
      
      std::map<G4String, std::map<G4String, G4double> >* E_byPart = event->GetE_byParticle();
      std::map<G4String, std::map<G4String, G4double> >::iterator Eiter;
      for (Eiter = E_byPart->begin(); Eiter != E_byPart->end(); ++Eiter) {
        std::map<G4String, G4double>::iterator partiter;
        for (partiter = E_byPart->at((*Eiter).first).begin(); partiter != E_byPart->at((*Eiter).first).end(); partiter++) {
          E_byParticle[(*Eiter).first][(*partiter).first] += (*partiter).second;
        }
      }
      
      std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* E_byPartAndGlobTimeMed = event->GetE_byParticleAndGlobalTimeMed();
      std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >::iterator ETimeMediter;
      for (ETimeMediter = E_byPartAndGlobTimeMed->begin(); ETimeMediter != E_byPartAndGlobTimeMed->end(); ETimeMediter++) {
        std::map<G4int, std::map<G4String, G4double> >::iterator timeiter;        
        for (timeiter = ((*ETimeMediter).second).begin(); timeiter != ((*ETimeMediter).second).end(); timeiter++) {
          std::map<G4String, G4double>::iterator partiter;
          for (partiter = ((*timeiter).second).begin(); partiter != ((*timeiter).second).end(); ++partiter) {
            E_byParticleAndGlobalTimeMed[(*ETimeMediter).first][(*timeiter).first][(*partiter).first] += (*partiter).second;
          }
        }
      }
      
      
      //
      // histogram the relative energy contribution by particle:
      //        
      //G4cout << ">>> histogram the relative energy contribution by particle:" << G4endl;
      
      for (Eiter = E_byPart->begin(); Eiter != E_byPart->end(); Eiter++) {
        G4String volName = (*Eiter).first;
        VolumeDirset[volName]->cd();
        if (EvtEnergybySDandpart.find(volName) == EvtEnergybySDandpart.end()) {
          std::map<G4String, TH1F*> tmpmap;
          EvtEnergybySDandpart.insert(std::make_pair(volName, tmpmap));
          std::map<G4String, G4double>::iterator partiter;
          for (partiter = E_byPart->at(volName).begin(); partiter != E_byPart->at(volName).end(); partiter++) {
            const char* hisname = (*partiter).first.c_str();
            TH1F *histo = new TH1F(hisname, "energy by particle", 200, 0, 100);
            EvtEnergybySDandpart[volName].insert(std::make_pair((*partiter).first, histo));
            EvtEnergybySDandpart[volName].at((*partiter).first)->Fill(((*partiter).second * 0.1) / Ein);
          }
        } else {
          std::map<G4String, G4double>::iterator partiter;
          for (partiter = E_byPart->at(volName).begin(); partiter != E_byPart->at(volName).end(); partiter++) {
            G4String PName = (*partiter).first;
            if (EvtEnergybySDandpart[volName].find(PName) == EvtEnergybySDandpart[volName].end()) {
              const char* hisname = PName.c_str();
              TH1F *histo = new TH1F(hisname, "energy by particle", 200, 0, 100);
              EvtEnergybySDandpart[volName].insert(std::make_pair(PName, histo));
              EvtEnergybySDandpart[volName].at(PName)->Fill(((*partiter).second * 0.1) / Ein);
            } else {
              EvtEnergybySDandpart[volName].at(PName)->Fill(((*partiter).second * 0.1) / Ein);
            }
          }
        }
      }
      
      
      //------------------------------------------------------------------------------------------------------
      std::map<G4String, std::map<G4String, G4double> >* Eobs_byPart = event->GetEobs_byParticle();
      std::map<G4String, std::map<G4String, G4double> >::iterator Eobsiter;
      for (Eobsiter = Eobs_byPart->begin(); Eobsiter != Eobs_byPart->end(); Eobsiter++) {
        std::map<G4String, G4double>::iterator partobsiter;
        for (partobsiter = Eobs_byPart->at((*Eobsiter).first).begin(); partobsiter != Eobs_byPart->at((*Eobsiter).first).end(); partobsiter++) {
          Eobs_byParticle[(*Eobsiter).first][(*partobsiter).first] += (*partobsiter).second;
        }
      }
      
      
      std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byPartAndGlobTimeLow = event->GetEobs_byParticleAndGlobalTimeLow();
      std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byPartAndGlobTimeMed = event->GetEobs_byParticleAndGlobalTimeMed();
      std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byPartAndGlobTimeHig = event->GetEobs_byParticleAndGlobalTimeHig();
      
      std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byPartAndLocTime1Low = event->GetEobs_byParticleAndLocalTime1Low();
      std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byPartAndLocTime1Med = event->GetEobs_byParticleAndLocalTime1Med();
      std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byPartAndLocTime1Hig = event->GetEobs_byParticleAndLocalTime2Hig();
      
      std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byPartAndLocTime2Low = event->GetEobs_byParticleAndLocalTime2Low();
      std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byPartAndLocTime2Med = event->GetEobs_byParticleAndLocalTime2Med();
      std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byPartAndLocTime2Hig = event->GetEobs_byParticleAndLocalTime2Hig();
      
      
      for(int sliceLowIt = 0; sliceLowIt <= nSlicesLow+1; ++sliceLowIt)
      {
        G4double time = sliceLowIt * timeSliceLow - timeSliceLow;
        G4double sum_globTime = 0.;
        G4double sum_locTime1 = 0.;
        G4double sum_locTime2 = 0.;
        for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
        {
          std::string partName = particleList.at(partIt);
          for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
          {
            std::string volName = volumeList.at(volIt);
            sum_globTime += (*Eobs_byPartAndGlobTimeLow)[volName][sliceLowIt][partName];
            sum_locTime1 += (*Eobs_byPartAndLocTime1Low)[volName][sliceLowIt][partName];
            sum_locTime2 += (*Eobs_byPartAndLocTime2Low)[volName][sliceLowIt][partName];
          }
        }
        EvtEnergybyglobaltimelow -> Fill(time,sum_globTime);
        EvtEnergybylocaltime1low -> Fill(time,sum_locTime1);
        EvtEnergybylocaltime2low -> Fill(time,sum_locTime2);
      }
      for(int sliceMedIt = 0; sliceMedIt <= nSlicesMed+1; ++sliceMedIt)
      {
        G4double time = sliceMedIt * timeSliceMed - timeSliceMed;
        G4double sum_globTime = 0.;
        G4double sum_locTime1 = 0.;
        G4double sum_locTime2 = 0.;
        for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
        {
          std::string partName = particleList.at(partIt);
          for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
          {
            std::string volName = volumeList.at(volIt);
            sum_globTime += (*Eobs_byPartAndGlobTimeMed)[volName][sliceMedIt][partName];
            sum_locTime1 += (*Eobs_byPartAndLocTime1Med)[volName][sliceMedIt][partName];
            sum_locTime2 += (*Eobs_byPartAndLocTime2Med)[volName][sliceMedIt][partName];
          }
        }
        EvtEnergybyglobaltimemed -> Fill(time,sum_globTime);
        EvtEnergybylocaltime1med -> Fill(time,sum_locTime1);
        EvtEnergybylocaltime2med -> Fill(time,sum_locTime2);
      }  
      for(int sliceHigIt = 0; sliceHigIt <= nSlicesHig+1; ++sliceHigIt)
      {
        G4double time = sliceHigIt * timeSliceHig - timeSliceHig;
        G4double sum_globTime = 0.;
        G4double sum_locTime1 = 0.;
        G4double sum_locTime2 = 0.;
        for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
        {
          std::string partName = particleList.at(partIt);
          for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
          {
            std::string volName = volumeList.at(volIt);
            sum_globTime += (*Eobs_byPartAndGlobTimeHig)[volName][sliceHigIt][partName];
            sum_locTime1 += (*Eobs_byPartAndLocTime1Hig)[volName][sliceHigIt][partName];
            sum_locTime2 += (*Eobs_byPartAndLocTime2Hig)[volName][sliceHigIt][partName];
          }
        }
          EvtEnergybyglobaltimehig -> Fill(time,sum_globTime);
          EvtEnergybylocaltime1hig -> Fill(time,sum_locTime1);
          EvtEnergybylocaltime2hig -> Fill(time,sum_locTime2);
      }
      
      
      for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
      {
        std::string partName = particleList.at(partIt);
        //G4cout << "partName: " << partName << G4endl;
        
        for(int sliceLowIt = 0; sliceLowIt <= nSlicesLow+1; ++sliceLowIt)
        {
          G4double time = sliceLowIt * timeSliceLow - timeSliceLow;
          G4double sum_globTime = 0.;
          G4double sum_locTime1 = 0.;
          G4double sum_locTime2 = 0.;
          for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
          {
            std::string volName = volumeList.at(volIt);
            sum_globTime += (*Eobs_byPartAndGlobTimeLow)[volName][sliceLowIt][partName];
            sum_locTime1 += (*Eobs_byPartAndLocTime1Low)[volName][sliceLowIt][partName];
            sum_locTime2 += (*Eobs_byPartAndLocTime2Low)[volName][sliceLowIt][partName];
          }
          EvtEnergybypartandglobaltimelow[partName] -> Fill(time,sum_globTime);
          EvtEnergybypartandlocaltime1low[partName] -> Fill(time,sum_locTime1);
          EvtEnergybypartandlocaltime2low[partName] -> Fill(time,sum_locTime2);
        }
        for(int sliceMedIt = 0; sliceMedIt <= nSlicesMed+1; ++sliceMedIt)
        {
          G4double time = sliceMedIt * timeSliceMed - timeSliceMed;
          G4double sum_globTime = 0.;
          G4double sum_locTime1 = 0.;
          G4double sum_locTime2 = 0.;
          for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
          {
            std::string volName = volumeList.at(volIt);
            sum_globTime += (*Eobs_byPartAndGlobTimeMed)[volName][sliceMedIt][partName];
            sum_locTime1 += (*Eobs_byPartAndLocTime1Med)[volName][sliceMedIt][partName];
            sum_locTime2 += (*Eobs_byPartAndLocTime2Med)[volName][sliceMedIt][partName];
          }
          EvtEnergybypartandglobaltimemed[partName] -> Fill(time,sum_globTime);
          EvtEnergybypartandlocaltime1med[partName] -> Fill(time,sum_locTime1);
          EvtEnergybypartandlocaltime2med[partName] -> Fill(time,sum_locTime2);
        }
        
        for(int sliceHigIt = 0; sliceHigIt <= nSlicesHig+1; ++sliceHigIt)
        {
          G4double time = sliceHigIt * timeSliceHig - timeSliceHig;
          G4double sum_globTime = 0.;
          G4double sum_locTime1 = 0.;
          G4double sum_locTime2 = 0.;
          for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
          {
            std::string volName = volumeList.at(volIt);
            sum_globTime += (*Eobs_byPartAndGlobTimeHig)[volName][sliceHigIt][partName];
            sum_locTime1 += (*Eobs_byPartAndLocTime1Hig)[volName][sliceHigIt][partName];
            sum_locTime2 += (*Eobs_byPartAndLocTime2Hig)[volName][sliceHigIt][partName];
          }
          EvtEnergybypartandglobaltimehig[partName] -> Fill(time,sum_globTime);
          EvtEnergybypartandlocaltime1hig[partName] -> Fill(time,sum_locTime1);
          EvtEnergybypartandlocaltime2hig[partName] -> Fill(time,sum_locTime2);
        }
      }
      
      
      for(unsigned int volIt = 0; volIt < volumeList.size(); ++volIt)
      {
        std::string volName = volumeList.at(volIt);
        //G4cout << "volName: " << volName << G4endl;
        
        for(unsigned int partIt = 0; partIt < particleList.size(); ++partIt)
        {
          std::string partName = particleList.at(partIt);
          //G4cout << "partName: " << partName << G4endl;
          
          for(int sliceLowIt = 0; sliceLowIt <= nSlicesLow+1; ++sliceLowIt)
          {
            G4double time = sliceLowIt * timeSliceLow - timeSliceLow;
            EvtEnergybySDandpartandglobaltimelow[volName][partName] -> Fill(time,(*Eobs_byPartAndGlobTimeLow)[volName][sliceLowIt][partName]);
            EvtEnergybySDandpartandlocaltime1low[volName][partName] -> Fill(time,(*Eobs_byPartAndLocTime1Low)[volName][sliceLowIt][partName]);
            EvtEnergybySDandpartandlocaltime2low[volName][partName] -> Fill(time,(*Eobs_byPartAndLocTime2Low)[volName][sliceLowIt][partName]);
          }
          
          for(int sliceMedIt = 0; sliceMedIt <= nSlicesMed+1; ++sliceMedIt)
          {
            G4double time = sliceMedIt * timeSliceMed - timeSliceMed;
            EvtEnergybySDandpartandglobaltimemed[volName][partName] -> Fill(time,(*Eobs_byPartAndGlobTimeMed)[volName][sliceMedIt][partName]);
            EvtEnergybySDandpartandlocaltime1med[volName][partName] -> Fill(time,(*Eobs_byPartAndLocTime1Med)[volName][sliceMedIt][partName]);
            EvtEnergybySDandpartandlocaltime2med[volName][partName] -> Fill(time,(*Eobs_byPartAndLocTime2Med)[volName][sliceMedIt][partName]);
          }
          
          for(int sliceHigIt = 0; sliceHigIt <= nSlicesHig+1; ++sliceHigIt)
          {
            G4double time = sliceHigIt * timeSliceHig - timeSliceHig;
            EvtEnergybySDandpartandglobaltimehig[volName][partName] -> Fill(time,(*Eobs_byPartAndGlobTimeHig)[volName][sliceHigIt][partName]);
            EvtEnergybySDandpartandlocaltime1hig[volName][partName] -> Fill(time,(*Eobs_byPartAndLocTime1Hig)[volName][sliceHigIt][partName]);
            EvtEnergybySDandpartandlocaltime2hig[volName][partName] -> Fill(time,(*Eobs_byPartAndLocTime2Hig)[volName][sliceHigIt][partName]);
          }
        }
      }
      
      
      // std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >::iterator EobsTimeLowiter;
      // for (EobsTimeLowiter = Eobs_byPartAndGlobTimeLow->begin(); EobsTimeLowiter != Eobs_byPartAndGlobTimeLow->end(); EobsTimeLowiter++) {
      //   std::map<G4int, std::map<G4String, G4double> >::iterator timeiter;
      //   for (timeiter = ((*EobsTimeLowiter).second).begin(); timeiter != ((*EobsTimeLowiter).second).end(); timeiter++) {
      //     std::map<G4String, G4double>::iterator partobsiter;
      //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
      //       Eobs_byParticleAndGlobalTimeLow[(*EobsTimeLowiter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
      //       EvtEnergybySDandpartandglobaltimelow[(*EobsTimeLowiter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
      //     }
      //   }
      // }
      // std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >::iterator EobsTimeMediter;
      // for (EobsTimeMediter = Eobs_byPartAndGlobTimeMed->begin(); EobsTimeMediter != Eobs_byPartAndGlobTimeMed->end(); EobsTimeMediter++) {
      //   std::map<G4int, std::map<G4String, G4double> >::iterator timeiter;
      //   for (timeiter = ((*EobsTimeMediter).second).begin(); timeiter != ((*EobsTimeMediter).second).end(); timeiter++) {
      //     std::map<G4String, G4double>::iterator partobsiter;
      //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
      //       Eobs_byParticleAndGlobalTimeMed[(*EobsTimeMediter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
      //       EvtEnergybySDandpartandglobaltimemed[(*EobsTimeMediter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
      //     }
      //   }
      // }
      // std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >::iterator EobsTimeHigiter;
      // for (EobsTimeHigiter = Eobs_byPartAndGlobTimeHig->begin(); EobsTimeHigiter != Eobs_byPartAndGlobTimeHig->end(); EobsTimeHigiter++) {
      //   std::map<G4int, std::map<G4String, G4double> >::iterator timeiter;
      //   for (timeiter = ((*EobsTimeHigiter).second).begin(); timeiter != ((*EobsTimeHigiter).second).end(); timeiter++) {
      //     std::map<G4String, G4double>::iterator partobsiter;
      //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
      //       Eobs_byParticleAndGlobalTimeHig[(*EobsTimeHigiter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
      //       EvtEnergybySDandpartandglobaltimehig[(*EobsTimeHigiter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
      //     }
      //   }
      // }
      
      
      // for (EobsTimeLowiter = Eobs_byPartAndLocTime1Low->begin(); EobsTimeLowiter != Eobs_byPartAndLocTime1Low->end(); EobsTimeLowiter++) {
      //   std::map<G4int, std::map<G4String, G4double> >::iterator timeiter;
      //   for (timeiter = ((*EobsTimeLowiter).second).begin(); timeiter != ((*EobsTimeLowiter).second).end(); timeiter++) {
      //     std::map<G4String, G4double>::iterator partobsiter;
      //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
      //       Eobs_byParticleAndLocalTime1Low[(*EobsTimeLowiter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
      //       EvtEnergybySDandpartandlocaltime1low[(*EobsTimeLowiter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
      //     }
      //   }
      // }
      // for (EobsTimeMediter = Eobs_byPartAndLocTime1Med->begin(); EobsTimeMediter != Eobs_byPartAndLocTime1Med->end(); EobsTimeMediter++) {
      //   std::map<G4int, std::map<G4String, G4double> >::iterator timeiter;
      //   for (timeiter = ((*EobsTimeMediter).second).begin(); timeiter != ((*EobsTimeMediter).second).end(); timeiter++) {
      //     std::map<G4String, G4double>::iterator partobsiter;
      //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
      //       Eobs_byParticleAndLocalTime1Med[(*EobsTimeMediter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
      //       EvtEnergybySDandpartandlocaltime1med[(*EobsTimeMediter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
      //     }
      //   }
      // }
      // for (EobsTimeHigiter = Eobs_byPartAndLocTime1Hig->begin(); EobsTimeHigiter != Eobs_byPartAndLocTime1Hig->end(); EobsTimeHigiter++) {
      //   std::map<G4int, std::map<G4String, G4double> >::iterator timeiter;
      //   for (timeiter = ((*EobsTimeHigiter).second).begin(); timeiter != ((*EobsTimeHigiter).second).end(); timeiter++) {
      //     std::map<G4String, G4double>::iterator partobsiter;
      //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
      //       Eobs_byParticleAndLocalTime1Hig[(*EobsTimeHigiter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
      //       EvtEnergybySDandpartandlocaltime1hig[(*EobsTimeHigiter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
      //     }
      //   }
      // }
      
      
      // for (EobsTimeLowiter = Eobs_byPartAndLocTime2Low->begin(); EobsTimeLowiter != Eobs_byPartAndLocTime2Low->end(); EobsTimeLowiter++) {
      //   std::map<G4int, std::map<G4String, G4double> >::iterator timeiter;
      //   for (timeiter = ((*EobsTimeLowiter).second).begin(); timeiter != ((*EobsTimeLowiter).second).end(); timeiter++) {
      //     std::map<G4String, G4double>::iterator partobsiter;
      //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
      //       Eobs_byParticleAndLocalTime2Low[(*EobsTimeLowiter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
      //       EvtEnergybySDandpartandlocaltime2low[(*EobsTimeLowiter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
      //     }
      //   }
      // }
      // for (EobsTimeMediter = Eobs_byPartAndLocTime2Med->begin(); EobsTimeMediter != Eobs_byPartAndLocTime2Med->end(); EobsTimeMediter++) {
      //   std::map<G4int, std::map<G4String, G4double> >::iterator timeiter;
      //   for (timeiter = ((*EobsTimeMediter).second).begin(); timeiter != ((*EobsTimeMediter).second).end(); timeiter++) {
      //     std::map<G4String, G4double>::iterator partobsiter;
      //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
      //       Eobs_byParticleAndLocalTime2Med[(*EobsTimeMediter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
      //       EvtEnergybySDandpartandlocaltime2med[(*EobsTimeMediter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
      //     }
      //   }
      // }
      // for (EobsTimeHigiter = Eobs_byPartAndLocTime2Hig->begin(); EobsTimeHigiter != Eobs_byPartAndLocTime2Hig->end(); EobsTimeHigiter++) {
      //   std::map<G4int, std::map<G4String, G4double> >::iterator timeiter;
      //   for (timeiter = ((*EobsTimeHigiter).second).begin(); timeiter != ((*EobsTimeHigiter).second).end(); timeiter++) {
      //     std::map<G4String, G4double>::iterator partobsiter;
      //     for (partobsiter = ((*timeiter).second).begin(); partobsiter != ((*timeiter).second).end(); partobsiter++) {
      //       Eobs_byParticleAndLocalTime2Hig[(*EobsTimeHigiter).first][(*timeiter).first][(*partobsiter).first] += (*partobsiter).second;
      //       EvtEnergybySDandpartandlocaltime2hig[(*EobsTimeHigiter).first][(*partobsiter).first] -> Fill((*timeiter).first,(*partobsiter).second);            
      //     }
      //   }
      // }
      
      
      //
      // histogram the relative observed energy contribution by particle:
      //  
      //G4cout << ">>> histogram the relative observed energy contribution by particle:" << G4endl;
      
      for (Eobsiter = Eobs_byPart->begin(); Eobsiter != Eobs_byPart->end(); Eobsiter++) {
        G4String volName = (*Eobsiter).first;
        VolumeDirset[volName]->cd();
        if (EvtObsEnergybySDandpart.find(volName) == EvtObsEnergybySDandpart.end()) {
          std::map<G4String, TH1F*> tmpmap;
          EvtObsEnergybySDandpart.insert(std::make_pair(volName, tmpmap));
          std::map<G4String, G4double>::iterator partiter;
          for (partiter = Eobs_byPart->at(volName).begin(); partiter != Eobs_byPart->at(volName).end(); partiter++) {
            G4String hname = "Obs" + (*partiter).first;
            const char* hisname = hname.c_str();
            TH1F *histo = new TH1F(hisname, "observed energy contribution by particle", 200, 0, 100);
            EvtObsEnergybySDandpart[volName].insert(std::make_pair((*partiter).first, histo));
            EvtObsEnergybySDandpart[volName].at((*partiter).first)->Fill(((*partiter).second * 0.1) / Ein);
          }
        } else {
          std::map<G4String, G4double>::iterator partiter;
          for (partiter = Eobs_byPart->at(volName).begin(); partiter != Eobs_byPart->at(volName).end(); partiter++) {
            G4String PName = (*partiter).first;
            if (EvtObsEnergybySDandpart[volName].find(PName) == EvtObsEnergybySDandpart[volName].end()) {
              G4String hname = "Obs" + PName;
              const char* hisname = hname.c_str();
              TH1F *histo = new TH1F(hisname, "observed Energy contribution by particle", 200, 0, 100);
              EvtObsEnergybySDandpart[volName].insert(std::make_pair(PName, histo));
              EvtObsEnergybySDandpart[volName].at(PName)->Fill(((*partiter).second * 0.1) / Ein);
            } else {
              EvtObsEnergybySDandpart[volName].at(PName)->Fill(((*partiter).second * 0.1) / Ein);
            }
          }
        }
      }
      
      std::map<G4String, std::map<G4String, G4double> >* NCeren_byPart = event->GetNCeren_byParticle();
      std::map<G4String, std::map<G4String, G4double> >::iterator Niter;
      for (Niter = NCeren_byPart->begin(); Niter != NCeren_byPart->end(); Niter++) {
        if (NCeren_byParticle.find((*Niter).first) == NCeren_byParticle.end()) {
          NCeren_byParticle[(*Niter).first] = (*Niter).second;
        } else {
          std::map<G4String, G4double>::iterator partiter;
          for (partiter = NCeren_byPart->at((*Niter).first).begin(); partiter != NCeren_byPart->at((*Niter).first).end(); partiter++) {
            NCeren_byParticle[(*Niter).first][(*partiter).first] = NCeren_byParticle[(*Niter).first][(*partiter).first] + (*partiter).second;
          }
        }
      }
      
      //
      // histogram the relative NCeren contribution by particle:
      //  
      //G4cout << ">>> histogram the relative NCeren contribution by particle:" << G4endl;
      
      for (Niter = NCeren_byPart->begin(); Niter != NCeren_byPart->end(); Niter++) {
        G4String volName = (*Niter).first;
        VolumeDirset[volName]->cd();
        if (EvtNcerenbySDandpart.find(volName) == EvtNcerenbySDandpart.end()) {
          std::map<G4String, TH1F*> tmpmap;
          EvtNcerenbySDandpart.insert(std::make_pair(volName, tmpmap));
          std::map<G4String, G4double>::iterator partiter;
          for (partiter = NCeren_byPart->at(volName).begin(); partiter != NCeren_byPart->at(volName).end(); partiter++) {
            G4String hname = "NC" + (*partiter).first;
            const char* hisname = hname.c_str();
            TH1F *histo = new TH1F(hisname, "nceren contribution by particle", 200, 0, 100);
            EvtNcerenbySDandpart[volName].insert(std::make_pair((*partiter).first, histo));
            EvtNcerenbySDandpart[volName].at((*partiter).first)->Fill(100.*(*partiter).second / event->GetNCeren());
          }
        } else {
          std::map<G4String, G4double>::iterator partiter;
          for (partiter = NCeren_byPart->at(volName).begin(); partiter != NCeren_byPart->at(volName).end(); partiter++) {
            G4String PName = (*partiter).first;
            if (EvtNcerenbySDandpart[volName].find(PName) == EvtNcerenbySDandpart[volName].end()) {
              G4String hname = "NC" + PName;
              const char* hisname = hname.c_str();
              TH1F *histo = new TH1F(hisname, "nceren contribution by particle", 200, 0, 100);
              EvtNcerenbySDandpart[volName].insert(std::make_pair(PName, histo));
              EvtNcerenbySDandpart[volName].at(PName)->Fill(100.*(*partiter).second / event->GetNCeren());
            } else {
              EvtNcerenbySDandpart[volName].at(PName)->Fill(100.*(*partiter).second / event->GetNCeren());
            }
          }
        }
      }
      
      std::map<G4String, std::map<G4String, G4int > >*procmult = event->GetProcessMult();
      std::map<G4String, std::map<G4String, G4int> >::iterator prociter;
      for (prociter = procmult->begin(); prociter != procmult->end(); prociter++) {
        if (processMult.find((*prociter).first) == processMult.end()) {
          processMult[(*prociter).first] = (*prociter).second;
        } else {
          std::map<G4String, G4int>::iterator partiter;
          for (partiter = procmult->at((*prociter).first).begin(); partiter != procmult->at((*prociter).first).end(); partiter++) {
            processMult[(*prociter).first][(*partiter).first] = processMult[(*prociter).first][(*partiter).first] + (*partiter).second;
          }
        }
      }
      
      std::map<G4String, std::map<G4String, std::map<G4String, products > > >*procmap = event->GetProcessMap();
      std::map<G4String, std::map<G4String, std::map<G4String, products> > >::iterator SDiter; // iterator over sensitive detectors
      std::map<G4String, std::map<G4String, products> >::iterator procmapiter; // Iterator over processes
      std::map<G4String, products>::iterator partmpiter; // Iterator over particles
      
      for (SDiter = procmap->begin(); SDiter != procmap->end(); SDiter++) { // check if sensitive det already included 
        if (processMap.find((*SDiter).first) == processMap.end()) {
          processMap[(*SDiter).first] = (*SDiter).second;
        } else {
          for (procmapiter = (*SDiter).second.begin(); procmapiter != (*SDiter).second.end(); procmapiter++) {
            if (processMap[(*SDiter).first].find((*procmapiter).first) == processMap[(*SDiter).first].end()) // check if process already included  
            {
              processMap[(*SDiter).first][(*procmapiter).first] = (*procmapiter).second;
            } else {
              for (partmpiter = (*procmapiter).second.begin(); partmpiter != (*procmapiter).second.end(); partmpiter++) {
                if (processMap[(*SDiter).first][(*procmapiter).first].find((*partmpiter).first)
                    == processMap[(*SDiter).first][(*procmapiter).first].end()) {
                  processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first] = (*partmpiter).second;
                } else {
                  processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first].NParticles =
                    processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first].NParticles
                    + (*partmpiter).second.NParticles;
                  processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first].kinE =
                    processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first].kinE
                    + (*partmpiter).second.kinE;
                  processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first].totE =
                    processMap[(*SDiter).first][(*procmapiter).first][(*partmpiter).first].totE
                    + (*partmpiter).second.totE;
                }
              }
            }
          }
        }
      }
      
      //-------------------------------------------
      //
      // Now we deal with the Hit Collections.
      //
      //-------------------------------------------
      //G4cout << ">>> Now we deal with the Hit Collections." << G4endl;
      
      std::map<G4String, std::vector<G4VHit*> >::iterator hciter;
      for (hciter = hcMap->begin(); hciter != hcMap->end(); hciter++) {
        std::vector<G4VHit*> hits = (*hciter).second;
        G4int NbHits = hits.size();
        //G4cout << "NbHits: " << NbHits << G4endl;
        std::vector<std::string> y = split((*hciter).first, '_');
        std::string Classname = y[1];
        std::string Colltype = y[2];
        std::string volName = y[0] + "_" + y[1];
        EvtEnergy[volName] = 0.0;
        EvtObsEnergy[volName] = 0.0;
        EvtNceren[volName] = 0.0;
        
        if (RunEnergy.find(volName) == RunEnergy.end()) {
          RunEnergy[volName] = 0.0;
        }
        if (RunObsEnergy.find(volName) == RunObsEnergy.end()) {
          RunObsEnergy[volName] = 0.0;
        }
        if (RunNceren.find(volName) == RunNceren.end()) {
          RunNceren[volName] = 0.0;
        }
        if (Classname == "DRCalorimeter") {
          for (G4int ii = 0; ii < NbHits; ii++) {
            DRCalorimeterHit* DRHit = dynamic_cast<DRCalorimeterHit*> (hits[ii]);
            tot_en = tot_en + DRHit->GetEdep();
            tot_nceren = tot_nceren + DRHit->GetNCeren();
            EvtEnergy[volName] = EvtEnergy[volName] + DRHit->GetEdep();
            EvtNceren[volName] = EvtNceren[volName] + DRHit->GetNCeren();
          }
        } else if (Classname == "DRTSCalorimeter" && Colltype == "HC") {
          for (G4int ii = 0; ii < NbHits; ii++) {
            DRTSCalorimeterHit* DRHit = dynamic_cast<DRTSCalorimeterHit*> (hits[ii]);
            tot_en = tot_en + DRHit->GetEdep();
            obs_en = obs_en + DRHit->GetEobsbirks();
            tot_nceren = tot_nceren + DRHit->GetNCeren();
            EvtEnergy[volName] = EvtEnergy[volName] + DRHit->GetEdep();
            EvtObsEnergy[volName] = EvtObsEnergy[volName] + DRHit->GetEobsbirks();
            EvtNceren[volName] = EvtNceren[volName] + DRHit->GetNCeren();
            h2_xy_firstTime_all -> Fill(DRHit->GetPos().x(),DRHit->GetPos().y(),DRHit->GetGlobalTime());
            h2_zy_firstTime_all -> Fill(DRHit->GetPos().z(),DRHit->GetPos().y(),DRHit->GetGlobalTime());
            h2_xy_firstTime[volName] -> Fill(DRHit->GetPos().x(),DRHit->GetPos().y(),DRHit->GetGlobalTime());
            h2_zy_firstTime[volName] -> Fill(DRHit->GetPos().z(),DRHit->GetPos().y(),DRHit->GetGlobalTime());
            //std::vector<DRTimeSliceHit>* tshc = DRHit->GetDRGlobalTimeSliceHitCol();
            // for (unsigned int iii = 0; iii < tshc->size(); iii++) {
            //   tedephisto->Fill((tshc->at(iii)).GetSlice(), (tshc->at(iii)).GetEdep());
            //   tedepEMhisto->Fill((tshc->at(iii)).GetSlice(), (tshc->at(iii)).GetEdepEM());
            //   tedepnEMhisto->Fill((tshc->at(iii)).GetSlice(), (tshc->at(iii)).GetEdepnonEM());
            //   tNChisto->Fill((tshc->at(i)).GetSlice(), (tshc->at(iii)).GetNCeren());
            // }
          }
        } else if (Classname == "Calorimeter") {
          for (G4int ii = 0; ii < NbHits; ii++) {
            CalorimeterHit* Hit = dynamic_cast<CalorimeterHit*> (hits[ii]);
            tot_en = tot_en + Hit->GetEdep();
          }
        } else if (Classname == "Tracker") {
          for (G4int ii = 0; ii < NbHits; ii++) {
            TrackerHit* THit = dynamic_cast<TrackerHit*> (hits[ii]);
            tot_en = tot_en + THit->GetEdep();
          }
        } else if (Classname == "PhotonDetector") {
          tot_nceren = tot_nceren + NbHits;
          G4cout << "Photon Detector hits: " << NbHits << G4endl;
          for (G4int ii = 0; ii < NbHits; ii++) {
            PhotonHit* PHit = dynamic_cast<PhotonHit*> (hits[ii]);
            PHit->Print(); // HJW just for now 
          }
        } else {
          //G4cout << "SD type: " << Classname << " unknown" << G4endl;
        }
      } // end loop over Hit collections
      tot_en = tot_en * 0.001;
      obs_en = obs_en * 0.001;
      Ehisto->Fill(tot_en);
      Eobshisto->Fill(obs_en);
      Nhisto->Fill(tot_nceren);
      std::map<G4String, G4double>::iterator e_itere;
      for (e_itere = EvtEnergy.begin(); e_itere != EvtEnergy.end(); e_itere++) {
        G4String VolumeName = (*e_itere).first;
        if (EHistosbySD.find(VolumeName) == EHistosbySD.end()) {
          G4String hname = "E_" + (*e_itere).first;
          const char* hisname = hname.c_str();
          TH1F* tmphisto = new TH1F(hisname, "Total energy deposited in", 100, 0.0, Ein);
          EHistosbySD[VolumeName].insert(std::make_pair(VolumeName, tmphisto));
          EHistosbySD[VolumeName].at(VolumeName)->Fill(EvtEnergy[VolumeName]*0.001);
        } else {
          EHistosbySD[VolumeName].at(VolumeName)->Fill(EvtEnergy[VolumeName]*0.001);
        }
      }
      for (e_itere = EvtObsEnergy.begin(); e_itere != EvtObsEnergy.end(); e_itere++) {
        G4String VolumeName = (*e_itere).first;
        if (EObsHistosbySD.find(VolumeName) == EObsHistosbySD.end()) {
          G4String hname = "EObs_" + (*e_itere).first;
          const char* hisname = hname.c_str();
          TH1F* tmphisto = new TH1F(hisname, "Total observed energy deposited in", 100, 0.0, Ein);
          EObsHistosbySD[VolumeName].insert(std::make_pair(VolumeName, tmphisto));
          EObsHistosbySD[VolumeName].at(VolumeName)->Fill(EvtObsEnergy[VolumeName]*0.001);
        } else {
          EObsHistosbySD[VolumeName].at(VolumeName)->Fill(EvtObsEnergy[VolumeName]*0.001);
        }
      }
      
      for (e_itere = EvtNceren.begin(); e_itere != EvtNceren.end(); e_itere++) {
        G4String VolumeName = (*e_itere).first;
        if (NCHistosbySD.find(VolumeName) == NCHistosbySD.end()) {
          G4String hname = "NC_" + (*e_itere).first;
          const char* hisname = hname.c_str();
          TH1F* tmphisto = new TH1F(hisname, "Total Nr of Cerenkov photons in", 100, 0.0, 65000 * (Ein + 0.1 * Ein));
          NCHistosbySD[VolumeName].insert(std::make_pair(VolumeName, tmphisto));
          NCHistosbySD[VolumeName].at(VolumeName)->Fill(EvtObsEnergy[VolumeName]*0.001);
        } else {
          NCHistosbySD[VolumeName].at(VolumeName)->Fill(EvtObsEnergy[VolumeName]*0.001);
        }
      }
      std::map<G4String, G4double>::iterator e_iter;
      for (e_iter = EvtEnergy.begin(); e_iter != EvtEnergy.end(); e_iter++) {
        RunEnergy[(*e_iter).first] = RunEnergy[(*e_iter).first]+(*e_iter).second;
      }
      std::map<G4String, G4double>::iterator o_iter;
      for (o_iter = EvtObsEnergy.begin(); o_iter != EvtObsEnergy.end(); o_iter++) {
        RunObsEnergy[(*o_iter).first] = RunObsEnergy[(*o_iter).first]+(*o_iter).second;
      }
      std::map<G4String, G4double>::iterator c_iter;
      for (c_iter = EvtObsEnergy.begin(); o_iter != EvtObsEnergy.end(); o_iter++) {
        RunObsEnergy[(*o_iter).first] = RunObsEnergy[(*o_iter).first]+(*o_iter).second;
      }
      
    } // end loop over events
    
    
    
    G4cout << "==============================================================================" << G4endl;
    G4cout << "=========    Average deposited Energy/Event:   " << std::setw(12) << RunTotEnergy / double(nevent) << "  [GeV]" << G4endl;
    G4cout << "=========    Average observed Energy/Event:    " << std::setw(12) << RunTotObsEnergy / double(nevent) << "  [GeV]" << G4endl;
    G4cout << "==============================================================================" << G4endl;
    
    std::map<G4String, std::map<G4String, G4double> >::iterator Eiter;
    for (Eiter = E_byParticle.begin(); Eiter != E_byParticle.end(); Eiter++) {
      double sum = 0.0;
      G4cout << "-----------------------------------------------" << G4endl;
      G4cout << "Sensitive Volume:  " << (*Eiter).first << G4endl;
      G4cout << "Total Energy: (per Evt)     " << RunEnergy[(*Eiter).first] / double(nevent) << G4endl;
      G4cout << "Observed Energy: (per Evt)  " << RunObsEnergy[(*Eiter).first] / double(nevent) << G4endl;
      G4cout << "-----------------------------------------------" << G4endl;
      std::map<G4String, G4double> ::iterator p_iter;
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
    
    std::map<G4String, std::map<G4String, G4double> >::iterator Eobsiter;
    for (Eobsiter = Eobs_byParticle.begin(); Eobsiter != Eobs_byParticle.end(); Eobsiter++) {
      G4cout << "==============================================================================" << G4endl;
      G4cout << "Sensitive Volume:  " << (*Eobsiter).first << G4endl;
      G4cout << "==============================================================================" << G4endl;
      double sumObs = 0.0;
      std::map<G4String, G4double> ::iterator pobs_iter;
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
      G4cout << "Average number of Cerenkov Photons/Event:" << RunTotNCeren / double(nevent) << G4endl;
      double sumN = 0.0;
      for (Eiter = NCeren_byParticle.begin(); Eiter != NCeren_byParticle.end(); Eiter++) {
        G4cout << "Sensitive Volume:  " << (*Eiter).first << G4endl;
        std::map<G4String, G4double> ::iterator p_iter;
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
    
    std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >::iterator EobsTimeMediter;
    for (EobsTimeMediter = Eobs_byParticleAndGlobalTimeMed.begin(); EobsTimeMediter != Eobs_byParticleAndGlobalTimeMed.end(); ++EobsTimeMediter++) {
      G4cout << "==============================================================================" << G4endl;
      G4cout << "Sensitive Volume:  " << (*EobsTimeMediter).first << G4endl;
      G4cout << "==============================================================================" << G4endl;
      std::map<G4int, std::map<G4String, G4double> >::iterator timeiter;
      for (timeiter = ((*EobsTimeMediter).second).begin(); timeiter != ((*EobsTimeMediter).second).end(); ++timeiter)
      {
        G4cout << ">>> global time slice (med): " << (*timeiter).first << G4endl;
        double sumObs = 0.0;
        std::map<G4String, G4double> ::iterator pobs_iter;
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
    //     G4cout << "Average number of Cerenkov Photons/Event:" << RunTotNCeren / double(nevent) << G4endl;
    //     double sumN = 0.0;
    //     for (Eiter = NCeren_byParticle.begin(); Eiter != NCeren_byParticle.end(); Eiter++) {
    //         G4cout << "Sensitive Volume:  " << (*Eiter).first << G4endl;
    //         std::map<G4String, G4double> ::iterator p_iter;
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
        G4cout << std::setw(20) << (*npartiter).first << std::setw(20) << (*npartiter).second / double(nevent) << G4endl;
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
                        << " #: " << std::setw(10) << (*pariter1).second.NParticles / double(nevent)
                        << " kinE: " << std::setw(10) << (*pariter1).second.kinE / double(nevent)
                        << " totE: " << std::setw(10) << (*pariter1).second.totE / double(nevent)
                        << G4endl;
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
