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
#include "DRTSCalorimeterSD.hh"
#include "TrackInformation.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Poisson.hh"
#include "G4VVisManager.hh"
#include "G4VProcess.hh"
#include "G4UserSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "DRTSCalorimeterSDMessenger.hh"
#include "G4PhysicalConstants.hh"
//
#include "RootIO.hh"
#include "EventAction.hh"
#include "Cerenkov.hh"
#include<iomanip>
#ifdef G4ANALYSIS_USE
#include "Analysis.hh"
Analysis* DRTSCalorimeterSD::analysis;
#endif

DRTSCalorimeterSD* DRTSCalorimeterSD::instance = 0;

EventAction* DRTSCalorimeterSD::EvtAction;

G4bool verbosity = false;



DRTSCalorimeterSD::DRTSCalorimeterSD(G4String name):
  G4VSensitiveDetector(name),
  particleList( new std::vector<G4String> ),
  particleListCeren( new std::vector<G4String> ),
  particleListShort( new std::vector<G4String> )

{
  pMessenger = new DRTSCalorimeterSDMessenger(this);
  birksc1 = 1.29e-2;
  birksc2 = 9.59e-6;
  CerenGenerator = new(Cerenkov);
  G4String HCname = name + "_HC";
  collectionName.insert(HCname);
  G4cout << collectionName.size() << "   DRTSCalorimeterSD name:  " << name << " collection Name: " << HCname << G4endl;
  HCname = name + "_HC2";
  collectionName.insert(HCname);
  G4cout << collectionName.size() << "   DRTSCalorimeterSD name:  " << name << " collection Name: " << HCname << G4endl;
  HCID = -1;
#ifdef G4ANALYSIS_USE
  analysis = Analysis::getInstance();
#endif
  
  // particleList -> push_back("fragment");
  // particleList -> push_back("He3");
  // particleList -> push_back("alpha");
  // particleList -> push_back("deuteron");
  // particleList -> push_back("triton");
  // particleList -> push_back("p_ev");
  // particleList -> push_back("p_sp");
  // particleList -> push_back("p_he");
  // particleList -> push_back("neutron");
  // particleList -> push_back("e+");
  // particleList -> push_back("e-");
  // particleList -> push_back("pi+");
  // particleList -> push_back("pi-");
  // particleList -> push_back("gamma");
  // particleList -> push_back("mu+");
  // particleList -> push_back("mu-");
  // particleList -> push_back("sigma+");
  // particleList -> push_back("sigma-");
  // particleList -> push_back("kaon+");
  // particleList -> push_back("kaon-");
  // particleList -> push_back("kaon0L");
  // particleList -> push_back("kaon0S");
  // particleList -> push_back("lambda");
  // particleList -> push_back("xi-");
  // particleList -> push_back("anti_neutron");
  // particleList -> push_back("anti_sigma-");
  // particleList -> push_back("anti_proton");
  // particleList -> push_back("anti_xi-");
  // particleList -> push_back("anti_omega-");
  // particleList -> push_back("anti_sigma+");
  // particleList -> push_back("anti_lambda");
  // particleList -> push_back("anti_xi0");
  // particleList -> push_back("other"); //Just in case
  
  particleList -> push_back("em");
  particleList -> push_back("fragment");
  particleList -> push_back("proton");
  particleList -> push_back("neutron");
  particleList -> push_back("e");
  particleList -> push_back("gamma");
  particleList -> push_back("mu");
  particleList -> push_back("pi");
  particleList -> push_back("meson");
  particleList -> push_back("baryon");
  particleList -> push_back("other"); //Just in case
  
  particleListShort -> push_back("em");
  particleListShort -> push_back("fragment");
  particleListShort -> push_back("proton");
  particleListShort -> push_back("neutron");
  particleListShort -> push_back("e");
  particleListShort -> push_back("gamma");
  particleListShort -> push_back("mu");
  particleListShort -> push_back("pi");
  particleListShort -> push_back("meson");
  particleListShort -> push_back("baryon");
  particleListShort -> push_back("other"); //Just in case
  
  particleListCeren -> push_back("e+");
  particleListCeren -> push_back("e-");
  particleListCeren -> push_back("kaon+");
  particleListCeren -> push_back("kaon-");
  particleListCeren -> push_back("mu+");
  particleListCeren -> push_back("mu-");
  particleListCeren -> push_back("pi+");
  particleListCeren -> push_back("pi-");
  particleListCeren -> push_back("proton");
  particleListCeren -> push_back("deuteron");
  particleListCeren -> push_back("triton");
  particleListCeren -> push_back("He3");
  particleListCeren -> push_back("sigma-");
  particleListCeren -> push_back("sigma+");
  particleListCeren -> push_back("xi-");
  particleListCeren -> push_back("anti_sigma-");
  particleListCeren -> push_back("anti_proton");
  particleListCeren -> push_back("anti_xi-");
  particleListCeren -> push_back("anti_omega-");
  particleListCeren -> push_back("anti_sigma+");
  particleListCeren -> push_back("other"); //Just in case
  
  instance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



DRTSCalorimeterSD::~DRTSCalorimeterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void DRTSCalorimeterSD::Initialize(G4HCofThisEvent* HCE)
{
  
#ifdef G4ANALYSIS_USE
  analysis->SDdir->cd();
  G4String SensitiveDetectorDirName = SensitiveDetectorName + "/";
  if (!(analysis->SDdir->GetDirectory(SensitiveDetectorDirName.c_str())))
  {
    mydir = analysis->SDdir->mkdir(SensitiveDetectorDirName.c_str());
    mydir->cd();
    dEdxweighted = new TH1F("dEdxweighted", "dEdx", 1000, 0, 4000);
    dEdxunweighted = new TH1F("dEdxunweighted", "dEdxunweighted", 1000, 0, 4000.);
    Birksweighted = new TH1F("Birksweighted", "Birksweighted", 1000, 0, 100);
    Birksunweighted = new TH1F("Birksunweighted", "Birksunweighted", 1000, 0, 100.);
    globalPositionX = new TH1F("globalPositionX","globalPositionX",20000,-2500.,2500.);
    globalPositionY = new TH1F("globalPositionY","globalPositionY",20000,-2500.,2500.);
    globalPositionZ = new TH1F("globalPositionZ","globalPositionZ",20000,-2500.,2500.);
    localPositionX = new TH1F("localPositionX","localPositionX",10000,-100.,100.);
    localPositionY = new TH1F("localPositionY","localPositionY",10000,-100.,100.);
    localPositionZ = new TH1F("localPositionZ","localPositionZ",10000,-100.,100.);
  }
#endif   
  
  drtscalorimeterCollection = new DRTSCalorimeterHitsCollection(SensitiveDetectorName, collectionName[0]);
  //if (HCID < 0) {
  HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  //}
  HCE->AddHitsCollection(HCID, drtscalorimeterCollection);
  
  drtscalorimeterCollection2 = new DRTSCalorimeterHits2Collection(SensitiveDetectorName, collectionName[1]);
  //if (HCID < 0) {
  HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[1]);
  //}
  HCE->AddHitsCollection(HCID, drtscalorimeterCollection2);
  
  
  
  EvtAction = EventAction::GetInstance();
  Event* CaTSEvt = EvtAction->GetEvent();
  
  std::map<G4String, std::map<G4String, G4int > > * processMult = CaTSEvt -> GetProcessMult();
  std::map<G4String, G4int> tmpmult;
  processMult -> insert(std::make_pair(SensitiveDetectorName, tmpmult));
  
  std::map<G4String, std::map<G4String, std::map<G4String, products> > > * processMap = CaTSEvt -> GetProcessMap();
  std::map<G4String, std::map<G4String, products> >tmpprocmap;
  processMap -> insert(std::make_pair(SensitiveDetectorName, tmpprocmap));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



G4bool DRTSCalorimeterSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  if( verbosity )
  {
    G4cout << ">>> DRTSCalorimeterSD::ProcessHits <<<" << G4endl;
  }
  
  EvtAction = EventAction::GetInstance();
  Event* CaTSEvt = EvtAction->GetEvent();
  
  const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
  const G4VSolid* solid = touch -> GetSolid();
  G4ThreeVector globalPosition = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector localPosition = touch->GetHistory()->GetTopTransform().TransformPoint(globalPosition);
  G4ThreeVector cellPosition = touch->GetTranslation();
  G4ThreeVector zAxisNeg(0.,0.,-1.);
  G4ThreeVector zAxisPos(0.,0.,+1.);
  G4double offsetPosition = 0.5 * fabs( solid->DistanceToOut(localPosition,zAxisNeg)+solid->DistanceToOut(localPosition,zAxisPos) );
  
  const G4double globalTime = aStep->GetPreStepPoint()->GetGlobalTime();
  const G4double localTime1 = globalTime - cellPosition.z()/c_light + offsetPosition/c_light;
  const G4double localTime2 = globalTime - sqrt(pow(cellPosition.x(),2)+pow(cellPosition.y(),2)+pow(cellPosition.z(),2))/c_light  + offsetPosition/c_light;
  
  G4int globalSliceLow = int(globalTime/timeslicelow) + 1;
  if( (globalTime-mintimelow) <  mintimelow ) globalSliceLow = 0;
  if( (globalTime-mintimelow) >= maxtimelow ) globalSliceLow = int((maxtimelow-mintimelow)/timeslicelow) + 1;
  G4int localSlice1Low = int(localTime1/timeslicelow) + 1;
  if( (localTime1-mintimelow) <  mintimelow ) localSlice1Low = 0;
  if( (localTime1-mintimelow) >= maxtimelow ) localSlice1Low = int((maxtimelow-mintimelow)/timeslicelow) + 1;
  G4int localSlice2Low = int(localTime2/timeslicelow) + 1;
  if( (localTime2-mintimelow) <  mintimelow ) localSlice2Low = 0;
  if( (localTime2-mintimelow) >= maxtimelow ) localSlice2Low = int((maxtimelow-mintimelow)/timeslicelow) + 1;
  
  G4int globalSliceMed = int(globalTime/timeslicemed) + 1;
  if( (globalTime-mintimemed) <  mintimemed ) globalSliceMed = 0;
  if( (globalTime-mintimemed) >= maxtimemed ) globalSliceMed = int((maxtimemed-mintimemed)/timeslicemed) + 1;
  G4int localSlice1Med = int(localTime1/timeslicemed) + 1;
  if( (localTime1-mintimemed) <  mintimemed ) localSlice1Med = 0;
  if( (localTime1-mintimemed) >= maxtimemed ) localSlice1Med = int((maxtimemed-mintimemed)/timeslicemed) + 1;
  G4int localSlice2Med = int(localTime2/timeslicemed) + 1;
  if( (localTime2-mintimemed) <  mintimemed ) localSlice2Med = 0;
  if( (localTime2-mintimemed) >= maxtimemed ) localSlice2Med = int((maxtimemed-mintimemed)/timeslicemed) + 1;
  
  G4int globalSliceHig = int(globalTime/timeslicehig) + 1;
  if( (globalTime-mintimehig) <  mintimehig ) globalSliceHig = 0;
  if( (globalTime-mintimehig) >= maxtimehig ) globalSliceHig = int((maxtimehig-mintimehig)/timeslicehig) + 1;
  G4int localSlice1Hig = int(localTime1/timeslicehig) + 1;
  if( (localTime1-mintimehig) <  mintimehig ) localSlice1Hig = 0;
  if( (localTime1-mintimehig) >= maxtimehig ) localSlice1Hig = int((maxtimehig-mintimehig)/timeslicehig) + 1;
  G4int localSlice2Hig = int(localTime2/timeslicehig) + 1;
  if( (localTime2-mintimehig) <  mintimehig ) localSlice2Hig = 0;
  if( (localTime2-mintimehig) >= maxtimehig ) localSlice2Hig = int((maxtimehig-mintimehig)/timeslicehig) + 1;
  
  globalPositionX -> Fill(globalPosition.x());
  globalPositionY -> Fill(globalPosition.y());
  globalPositionZ -> Fill(globalPosition.z());
  localPositionX -> Fill(localPosition.x());
  localPositionY -> Fill(localPosition.y());
  localPositionZ -> Fill(localPosition.z());

  
  
  //-------------------
  // GLOBAL INFORMATION
  //-------------------
  if( verbosity )
  {
    G4cout << "//------------------------" << G4endl;
    G4cout << "// GLOBAL INFORMATION"      << G4endl;
    G4cout << "//------------------------" << G4endl;
  }
  
  //---------------------------------------------------------------------------------------------------
  // Birks' Law 
  // ===========
  //    
  // In the case of Scintillator as active medium, we can describe the quenching effects with the Birks' law,
  // using the expression and the coefficients taken from the paper NIM 80 (1970) 239-244 for the organic
  // scintillator NE-102:
  //                      S*dE/dr
  //          dL/dr = --------------
  //                  1 + C1*(dE/dr)
  // with:
  //          S=1
  //          C1 = 1.29 x 10^-2  g*cm^-2*MeV^-1
  //          C2 = 9.59 x 10^-6  g^2*cm^-4*MeV^-2
  // These are the same values used by ATLAS TileCal and CMS HCAL (and also the default in Geant3).
  // You can try different values for these parameters, to have an idea on the uncertainties due to them,
  // by uncommenting one of the lines below. To get the "dE/dr" that appears in the formula,
  // which has the dimensions 
  //          [ dE/dr ] = MeV * cm^2 / g
  // we have to divide the energy deposit in MeV by the product of the step length (in cm) and the density
  // of the scintillator:
  //          rho_scintillator = 1.032  g/cm3 .
  // Of course, in our case we use only the denominator of the Birks' formula above, because we do not have
  // digitization, i.e. we only have energy deposit but not conversion to photons.
  // Birks should not be applied in the case of gamma energy depositions (which happens only for the
  // photoelectric process), because in this case the step length is related to the photoelectric cross
  // section, and not to the distance in which the energy is actually deposited, that is what is
  // needed in order to determine dE/dx which enters
  // in the Birks' law.
  // Similarly, for neutron energy depositions (which have been introduced in Geant4 8.1 as an effective
  // way to describe the elastic nuclei recoils below a certain kinetic threshold, set by default to
  // 100 keV), the step length is related to the neutron elastic cross-section, and not to the real ionization
  // range of the recoiled nuclei, which should instead be considered for the dE/dx in the Birks' law.
  // In the case of neutron energy depositions, the most correct way to apply the Birks quench would
  // be to eliminate them by setting the kinetic threshold to 0.0 (variable  "edepLimit"  in the
  // file  "G4HadronElasticPhysics.cc"  in  "geant4/physics_lists/hadronic/Packaging/src/" ),
  // so that the recoiled nuclei tracks are always generated explicitly. This, of course, costs in
  // term of CPU performance.
  
  G4double edep = aStep->GetTotalEnergyDeposit() / CLHEP::MeV;
  if (edep == 0.) return false;
  G4double rho = aStep->GetPreStepPoint()->GetMaterial()->GetDensity() / (CLHEP::g / CLHEP::cm3);
  G4double stepLength = aStep->GetStepLength() / CLHEP::cm;
  G4double birksFactor = 1.0;
  G4double dedx;
  G4double vel;
  if (aStep->GetTrack()->GetDefinition() == G4Proton::ProtonDefinition()) {
    vel = aStep->GetTrack()->GetVelocity();
  }
  //Do not apply Birks for gamma deposits!
  if (stepLength > 1.0e-8)//Check, cut if necessary.
  {
    dedx = edep / (rho * stepLength); //[MeV*cm^2/g]
    birksFactor = 1.0 / (1.0 + birksc1 * dedx + birksc2 * dedx * dedx);
    if (aStep->GetTrack()->GetDefinition() == G4Gamma::GammaDefinition()) {
      birksFactor = 1.0;
    }
    if (aStep->GetTrack()->GetDefinition() == G4Neutron::NeutronDefinition()) {
      birksFactor = 1.0;
    }
  }
  if (aStep->GetTrack()->GetDefinition() == G4Proton::ProtonDefinition()) {
    if (aStep->GetTrack()->GetTrackStatus() == fAlive) {
      dEdxunweighted->Fill(aStep->GetTrack()->GetMomentum().mag());
      dEdxweighted->Fill(aStep->GetTrack()->GetMomentum().mag(), dedx);
      Birksunweighted->Fill(aStep->GetTrack()->GetKineticEnergy());
      Birksweighted->Fill(aStep->GetTrack()->GetKineticEnergy(), birksFactor);
    }
  }
  G4double eobsbirks = edep*birksFactor;
  CaTSEvt->SetTotEnergy(CaTSEvt->GetTotEnergy() + edep);
  CaTSEvt->SetTotObsEnergy(CaTSEvt->GetTotObsEnergy() + eobsbirks);
  
  if( verbosity )
  {
    TrackInformation* aTrackInfo = (TrackInformation*)( aStep->GetTrack()->GetUserInformation() );
    G4cout << "isEM: " << aTrackInfo->GetParticleIsEM() << G4endl;
  }
  
  
  //------------------------
  // INFORMATION BY PARTICLE
  //------------------------
  if( verbosity )
  {
    G4cout << "//------------------------" << G4endl;
    G4cout << "// INFORMATION BY PARTICLE" << G4endl;
    G4cout << "//------------------------" << G4endl;
  }
  
  //-------------------
  // energy by particle
  
  G4Track* theTrack = aStep->GetTrack();
  G4String particleName = theTrack->GetDefinition()->GetParticleName();
  G4String fragment = "fragment";
  if (theTrack->GetParticleDefinition()->GetParticleType() == "nucleus" && theTrack->GetParticleDefinition()->GetParticleSubType() == "generic") {
    particleName = fragment;
  }
  if (particleName == "proton") {
    if (theTrack->GetKineticEnergy()/CLHEP::MeV < 10.) {
      particleName = "p_ev";
    } else if (theTrack->GetKineticEnergy()/CLHEP::MeV < 500.) {
      particleName = "p_sp";
    } else {
      particleName = "p_he";
    }
  }
  
  std::map<G4String, std::map<G4String, G4double> >* E_byParticle = CaTSEvt->GetE_byParticle();
  if (std::find(particleList->begin(),particleList->end(),particleName) == particleList->end()) {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   particleName: " << "other"
             << "   Edep_old: " << (*E_byParticle)[SensitiveDetectorName]["other"]
             << "   Edep_now: " << edep
             << G4endl;
    (*E_byParticle)[SensitiveDetectorName]["other"] += edep;
  } else {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   particleName: " << particleName
             << "   Edep_old: " << (*E_byParticle)[SensitiveDetectorName][particleName]
             << "   Edep_now: " << edep
             << G4endl;
    (*E_byParticle)[SensitiveDetectorName][particleName] += edep;
  }
  std::map<G4String, std::map<G4String, G4double> >* Eobs_byParticle = CaTSEvt->GetEobs_byParticle();
  if (std::find(particleList->begin(),particleList->end(),particleName) == particleList->end()) {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   particleName: " << "other"
             << "   Eobs_old: " << (*Eobs_byParticle)[SensitiveDetectorName]["other"]
             << "   Eobs_now: " << eobsbirks
             << G4endl;
    (*Eobs_byParticle)[SensitiveDetectorName]["other"] += eobsbirks;
  } else {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   particleName: " << particleName
             << "   Eobs_old: " << (*Eobs_byParticle)[SensitiveDetectorName][particleName]
             << "   Eobs_now: " << eobsbirks
             << G4endl;
    (*Eobs_byParticle)[SensitiveDetectorName][particleName] += eobsbirks;
  }
  
  //-----------------------------
  // cerenkov photons by particle
  
  particleName = theTrack->GetDefinition()->GetParticleName();
  G4int NCerenPhotons = 0;
  G4StepPoint* pPreStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep->GetPostStepPoint();
  const G4double charge = theTrack->GetDefinition()->GetPDGCharge();
  G4double beta = 0.5 * (pPreStepPoint->GetBeta() + pPostStepPoint->GetBeta());
  G4String thematerial = touch->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName();
  G4double MeanNumberOfPhotons = CerenGenerator->GetAverageNumberOfPhotons(charge, beta, thematerial);
  if (MeanNumberOfPhotons > 0.0) {
    G4double step_length = aStep->GetStepLength();
    MeanNumberOfPhotons = MeanNumberOfPhotons * step_length;
    NCerenPhotons = (G4int) G4Poisson(MeanNumberOfPhotons);
  } else {
    NCerenPhotons = 0;
  }
  if (NCerenPhotons > 0)
  {
    CaTSEvt->SetTotNCeren(CaTSEvt->GetNCeren() + double(NCerenPhotons));
    std::map<G4String, std::map<G4String, G4double > >* NCeren_byParticle = CaTSEvt->GetNCeren_byParticle();
    if (std::find(particleListCeren->begin(),particleListCeren->end(),particleName) == particleListCeren->end()) {
      (*NCeren_byParticle)[SensitiveDetectorName]["other"] += double(NCerenPhotons);
    }
    else
    {
      (*NCeren_byParticle)[SensitiveDetectorName][particleName] += double(NCerenPhotons);
    }
  }
  
  
  
  //---------------------------------
  // INFORMATION BY PARTICLE AND TIME
  //---------------------------------
  if( verbosity )
  {
    G4cout << "//---------------------------------" << G4endl;
    G4cout << "// INFORMATION BY PARTICLE AND TIME" << G4endl;
    G4cout << "//---------------------------------" << G4endl;
  }
  
  //----------------------------
  // energy by particle and time
  
  // if (theTrack->GetParticleDefinition()->GetParticleType() == "nucleus" && theTrack->GetParticleDefinition()->GetParticleSubType() == "generic") {
  //   particleName = fragment;
  // }
  // if (particleName == "proton") {
  //   if (theTrack->GetKineticEnergy()/CLHEP::MeV < 10.) {
  //     particleName = "p_ev";
  //   } else if (theTrack->GetKineticEnergy()/CLHEP::MeV < 500.) {
  //     particleName = "p_sp";
  //   } else {
  //     particleName = "p_he";
  //   }
  // }

  TrackInformation* aTrackInfo = (TrackInformation*)( aStep->GetTrack()->GetUserInformation() );  
  G4String partName = particleName;
  if ( theTrack->GetParticleDefinition()->GetParticleType() == "nucleus" ) partName = fragment;
  if ( theTrack->GetParticleDefinition()->GetParticleType() == "meson" )   partName = "meson";
  if ( theTrack->GetParticleDefinition()->GetParticleType() == "baryon" )  partName = "baryon";
  if ( particleName == "proton"  || particleName == "anti_proton" )        partName = "proton";
  if ( particleName == "neutron" || particleName == "anti_neutron" )       partName = "neutron";
  if ( particleName == "e-"      || particleName == "e+" )                 partName = "e";
  if ( particleName == "mu-"     || particleName == "mu+" )                partName = "mu";
  if ( particleName == "pi+"     || particleName == "pi-" )                partName = "pi";
  if( aTrackInfo->GetParticleIsEM() ) partName = "em";
  
  if ( std::find(particleListShort->begin(),particleListShort->end(),partName) == particleList->end() )
    partName = "other";
  
  particleName = partName;
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* E_byParticleAndGlobalTimeLow = CaTSEvt->GetE_byParticleAndGlobalTimeLow();
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* E_byParticleAndGlobalTimeMed = CaTSEvt->GetE_byParticleAndGlobalTimeMed();
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* E_byParticleAndGlobalTimeHig = CaTSEvt->GetE_byParticleAndGlobalTimeHig();
  if (std::find(particleList->begin(),particleList->end(),particleName) == particleList->end()) {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   globalTime (ns): " << globalTime/CLHEP::ns
             << "   globalTimeSliceLow: " << globalSliceLow
             << "   globalTimeSliceMed: " << globalSliceMed
             << "   globalTimeSliceHig: " << globalSliceHig
             << "   particleName: " << "other"
             << "   Edep_old (low): " << (*E_byParticleAndGlobalTimeLow)[SensitiveDetectorName][globalSliceLow]["other"]
             << "   Edep_old (med): " << (*E_byParticleAndGlobalTimeMed)[SensitiveDetectorName][globalSliceMed]["other"]
             << "   Edep_old (hig): " << (*E_byParticleAndGlobalTimeHig)[SensitiveDetectorName][globalSliceHig]["other"]
             << "   Edep_now: " << edep
             << G4endl;
    (*E_byParticleAndGlobalTimeLow)[SensitiveDetectorName][globalSliceLow]["other"] += edep;
    (*E_byParticleAndGlobalTimeMed)[SensitiveDetectorName][globalSliceMed]["other"] += edep;
    (*E_byParticleAndGlobalTimeHig)[SensitiveDetectorName][globalSliceHig]["other"] += edep;
  } else {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   globalTime (ns): " << globalTime/CLHEP::ns
             << "   globalTimeSliceLow: " << globalSliceLow
             << "   globalTimeSliceMed: " << globalSliceMed
             << "   globalTimeSliceHig: " << globalSliceHig
             << "   particleName: " << particleName
             << "   Edep_old (low): " << (*E_byParticleAndGlobalTimeLow)[SensitiveDetectorName][globalSliceLow][particleName]
             << "   Edep_old (med): " << (*E_byParticleAndGlobalTimeMed)[SensitiveDetectorName][globalSliceMed][particleName]
             << "   Edep_old (hig): " << (*E_byParticleAndGlobalTimeHig)[SensitiveDetectorName][globalSliceHig][particleName]
             << "   Edep_now: " << edep
             << G4endl;
    (*E_byParticleAndGlobalTimeLow)[SensitiveDetectorName][globalSliceLow][particleName] += edep;
    (*E_byParticleAndGlobalTimeMed)[SensitiveDetectorName][globalSliceMed][particleName] += edep;
    (*E_byParticleAndGlobalTimeHig)[SensitiveDetectorName][globalSliceHig][particleName] += edep;
  }
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byParticleAndGlobalTimeLow = CaTSEvt->GetEobs_byParticleAndGlobalTimeLow();
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byParticleAndGlobalTimeMed = CaTSEvt->GetEobs_byParticleAndGlobalTimeMed();
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byParticleAndGlobalTimeHig = CaTSEvt->GetEobs_byParticleAndGlobalTimeHig();
  if (std::find(particleList->begin(),particleList->end(),particleName) == particleList->end()) {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   globalTime (ns): " << globalTime/CLHEP::ns
             << "   globalTimeSliceLow: " << globalSliceLow
             << "   globalTimeSliceMed: " << globalSliceMed
             << "   globalTimeSliceHig: " << globalSliceHig
             << "   particleName: " << "other"
             << "   Eobs_old (low): " << (*Eobs_byParticleAndGlobalTimeLow)[SensitiveDetectorName][globalSliceLow]["other"]
             << "   Eobs_old (med): " << (*Eobs_byParticleAndGlobalTimeMed)[SensitiveDetectorName][globalSliceMed]["other"]
             << "   Eobs_old (hig): " << (*Eobs_byParticleAndGlobalTimeHig)[SensitiveDetectorName][globalSliceHig]["other"]
             << "   Eobs_now: " << eobsbirks
             << G4endl;
    (*Eobs_byParticleAndGlobalTimeLow)[SensitiveDetectorName][globalSliceLow]["other"] += eobsbirks;
    (*Eobs_byParticleAndGlobalTimeMed)[SensitiveDetectorName][globalSliceMed]["other"] += eobsbirks;
    (*Eobs_byParticleAndGlobalTimeHig)[SensitiveDetectorName][globalSliceHig]["other"] += eobsbirks;
  } else {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   globalTime (ns): " << globalTime/CLHEP::ns
             << "   globalTimeSliceLow: " << globalSliceLow
             << "   globalTimeSliceMed: " << globalSliceMed
             << "   globalTimeSliceHig: " << globalSliceHig
             << "   particleName: " << particleName
             << "   Eobs_old (low): " << (*Eobs_byParticleAndGlobalTimeLow)[SensitiveDetectorName][globalSliceLow][particleName]
             << "   Eobs_old (med): " << (*Eobs_byParticleAndGlobalTimeMed)[SensitiveDetectorName][globalSliceMed][particleName]
             << "   Eobs_old (hig): " << (*Eobs_byParticleAndGlobalTimeHig)[SensitiveDetectorName][globalSliceHig][particleName]
             << "   Eobs_now: " << eobsbirks
             << G4endl;
    (*Eobs_byParticleAndGlobalTimeLow)[SensitiveDetectorName][globalSliceLow][particleName] += eobsbirks;
    (*Eobs_byParticleAndGlobalTimeMed)[SensitiveDetectorName][globalSliceMed][particleName] += eobsbirks;
    (*Eobs_byParticleAndGlobalTimeHig)[SensitiveDetectorName][globalSliceHig][particleName] += eobsbirks;
  }
  
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* E_byParticleAndLocalTime1Low = CaTSEvt->GetE_byParticleAndLocalTime1Low();
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* E_byParticleAndLocalTime1Med = CaTSEvt->GetE_byParticleAndLocalTime1Med();
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* E_byParticleAndLocalTime1Hig = CaTSEvt->GetE_byParticleAndLocalTime1Hig();
  if (std::find(particleList->begin(),particleList->end(),particleName) == particleList->end()) {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   localTime1 (ns): " << localTime1/CLHEP::ns
             << "   localTime1SliceLow: " << localSlice1Low
             << "   localTime1SliceMed: " << localSlice1Med
             << "   localTime1SliceHig: " << localSlice1Hig
             << "   particleName: " << "other"
             << "   Edep_old (low): " << (*E_byParticleAndLocalTime1Low)[SensitiveDetectorName][localSlice1Low]["other"]
             << "   Edep_old (med): " << (*E_byParticleAndLocalTime1Med)[SensitiveDetectorName][localSlice1Med]["other"]
             << "   Edep_old (hig): " << (*E_byParticleAndLocalTime1Hig)[SensitiveDetectorName][localSlice1Hig]["other"]
             << "   Edep_now: " << edep
             << G4endl;
    (*E_byParticleAndLocalTime1Low)[SensitiveDetectorName][localSlice1Low]["other"] += edep;
    (*E_byParticleAndLocalTime1Med)[SensitiveDetectorName][localSlice1Med]["other"] += edep;
    (*E_byParticleAndLocalTime1Hig)[SensitiveDetectorName][localSlice1Hig]["other"] += edep;
  } else {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   localTime1 (ns): " << localTime1/CLHEP::ns
             << "   localTime1SliceLow: " << localSlice1Low
             << "   localTime1SliceMed: " << localSlice1Med
             << "   localTime1SliceHig: " << localSlice1Hig
             << "   particleName: " << particleName
             << "   Edep_old (low): " << (*E_byParticleAndLocalTime1Low)[SensitiveDetectorName][localSlice1Low][particleName]
             << "   Edep_old (med): " << (*E_byParticleAndLocalTime1Med)[SensitiveDetectorName][localSlice1Med][particleName]
             << "   Edep_old (hig): " << (*E_byParticleAndLocalTime1Hig)[SensitiveDetectorName][localSlice1Hig][particleName]
             << "   Edep_now: " << edep
             << G4endl;
    (*E_byParticleAndLocalTime1Low)[SensitiveDetectorName][localSlice1Low][particleName] += edep;
    (*E_byParticleAndLocalTime1Med)[SensitiveDetectorName][localSlice1Med][particleName] += edep;
    (*E_byParticleAndLocalTime1Hig)[SensitiveDetectorName][localSlice1Hig][particleName] += edep;
  }
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byParticleAndLocalTime1Low = CaTSEvt->GetEobs_byParticleAndLocalTime1Low();
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byParticleAndLocalTime1Med = CaTSEvt->GetEobs_byParticleAndLocalTime1Med();
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byParticleAndLocalTime1Hig = CaTSEvt->GetEobs_byParticleAndLocalTime1Hig();
  if (std::find(particleList->begin(),particleList->end(),particleName) == particleList->end()) {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   localTime1 (ns): " << localTime1/CLHEP::ns
             << "   localTime1SliceLow: " << localSlice1Low
             << "   localTime1SliceMed: " << localSlice1Med
             << "   localTime1SliceHig: " << localSlice1Hig
             << "   particleName: " << "other"
             << "   Eobs_old (low): " << (*Eobs_byParticleAndLocalTime1Low)[SensitiveDetectorName][localSlice1Low]["other"]
             << "   Eobs_old (med): " << (*Eobs_byParticleAndLocalTime1Med)[SensitiveDetectorName][localSlice1Med]["other"]
             << "   Eobs_old (hig): " << (*Eobs_byParticleAndLocalTime1Hig)[SensitiveDetectorName][localSlice1Hig]["other"]
             << "   Eobs_now: " << eobsbirks
             << G4endl;
    (*Eobs_byParticleAndLocalTime1Low)[SensitiveDetectorName][localSlice1Low]["other"] += eobsbirks;
    (*Eobs_byParticleAndLocalTime1Med)[SensitiveDetectorName][localSlice1Med]["other"] += eobsbirks;
    (*Eobs_byParticleAndLocalTime1Hig)[SensitiveDetectorName][localSlice1Hig]["other"] += eobsbirks;
  } else {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   localTime1 (ns): " << localTime1/CLHEP::ns
             << "   localTime1SliceLow: " << localSlice1Low
             << "   localTime1SliceMed: " << localSlice1Med
             << "   localTime1SliceHig: " << localSlice1Hig
             << "   particleName: " << particleName
             << "   Eobs_old (low): " << (*Eobs_byParticleAndLocalTime1Low)[SensitiveDetectorName][localSlice1Low][particleName]
             << "   Eobs_old (med): " << (*Eobs_byParticleAndLocalTime1Med)[SensitiveDetectorName][localSlice1Med][particleName]
             << "   Eobs_old (hig): " << (*Eobs_byParticleAndLocalTime1Hig)[SensitiveDetectorName][localSlice1Hig][particleName]
             << "   Eobs_now: " << eobsbirks
             << G4endl;
    (*Eobs_byParticleAndLocalTime1Low)[SensitiveDetectorName][localSlice1Low][particleName] += eobsbirks;
    (*Eobs_byParticleAndLocalTime1Med)[SensitiveDetectorName][localSlice1Med][particleName] += eobsbirks;
    (*Eobs_byParticleAndLocalTime1Hig)[SensitiveDetectorName][localSlice1Hig][particleName] += eobsbirks;
  }
  
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* E_byParticleAndLocalTime2Low = CaTSEvt->GetE_byParticleAndLocalTime2Low();
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* E_byParticleAndLocalTime2Med = CaTSEvt->GetE_byParticleAndLocalTime2Med();
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* E_byParticleAndLocalTime2Hig = CaTSEvt->GetE_byParticleAndLocalTime2Hig();
  if (std::find(particleList->begin(),particleList->end(),particleName) == particleList->end()) {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   localTime2 (ns): " << localTime2/CLHEP::ns
             << "   localTime2SliceLow: " << localSlice2Low
             << "   localTime2SliceMed: " << localSlice2Med
             << "   localTime2SliceHig: " << localSlice2Hig
             << "   particleName: " << "other"
             << "   Edep_old (low): " << (*E_byParticleAndLocalTime2Low)[SensitiveDetectorName][localSlice2Low]["other"]
             << "   Edep_old (med): " << (*E_byParticleAndLocalTime2Med)[SensitiveDetectorName][localSlice2Med]["other"]
             << "   Edep_old (hig): " << (*E_byParticleAndLocalTime2Hig)[SensitiveDetectorName][localSlice2Hig]["other"]
             << "   Edep_now: " << edep
             << G4endl;
    (*E_byParticleAndLocalTime2Low)[SensitiveDetectorName][localSlice2Low]["other"] += edep;
    (*E_byParticleAndLocalTime2Med)[SensitiveDetectorName][localSlice2Med]["other"] += edep;
    (*E_byParticleAndLocalTime2Hig)[SensitiveDetectorName][localSlice2Hig]["other"] += edep;
  } else {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   localTime2 (ns): " << localTime2/CLHEP::ns
             << "   localTime2SliceLow: " << localSlice2Low
             << "   localTime2SliceMed: " << localSlice2Med
             << "   localTime2SliceHig: " << localSlice2Hig
             << "   particleName: " << particleName
             << "   Edep_old (low): " << (*E_byParticleAndLocalTime2Low)[SensitiveDetectorName][localSlice2Low][particleName]
             << "   Edep_old (med): " << (*E_byParticleAndLocalTime2Med)[SensitiveDetectorName][localSlice2Med][particleName]
             << "   Edep_old (hig): " << (*E_byParticleAndLocalTime2Hig)[SensitiveDetectorName][localSlice2Hig][particleName]
             << "   Edep_now: " << edep
             << G4endl;
    (*E_byParticleAndLocalTime2Low)[SensitiveDetectorName][localSlice2Low][particleName] += edep;
    (*E_byParticleAndLocalTime2Med)[SensitiveDetectorName][localSlice2Med][particleName] += edep;
    (*E_byParticleAndLocalTime2Hig)[SensitiveDetectorName][localSlice2Hig][particleName] += edep;
  }
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byParticleAndLocalTime2Low = CaTSEvt->GetEobs_byParticleAndLocalTime2Low();
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byParticleAndLocalTime2Med = CaTSEvt->GetEobs_byParticleAndLocalTime2Med();
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* Eobs_byParticleAndLocalTime2Hig = CaTSEvt->GetEobs_byParticleAndLocalTime2Hig();
  if (std::find(particleList->begin(),particleList->end(),particleName) == particleList->end()) {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   localTime2 (ns): " << localTime2/CLHEP::ns
             << "   localTime2SliceLow: " << localSlice2Low
             << "   localTime2SliceMed: " << localSlice2Med
             << "   localTime2SliceHig: " << localSlice2Hig
             << "   particleName: " << "other"
             << "   Eobs_old (low): " << (*Eobs_byParticleAndLocalTime2Low)[SensitiveDetectorName][localSlice2Low]["other"]
             << "   Eobs_old (med): " << (*Eobs_byParticleAndLocalTime2Med)[SensitiveDetectorName][localSlice2Med]["other"]
             << "   Eobs_old (hig): " << (*Eobs_byParticleAndLocalTime2Hig)[SensitiveDetectorName][localSlice2Hig]["other"]
             << "   Eobs_now: " << eobsbirks
             << G4endl;
    (*Eobs_byParticleAndLocalTime2Low)[SensitiveDetectorName][localSlice2Low]["other"] += eobsbirks;
    (*Eobs_byParticleAndLocalTime2Med)[SensitiveDetectorName][localSlice2Med]["other"] += eobsbirks;
    (*Eobs_byParticleAndLocalTime2Hig)[SensitiveDetectorName][localSlice2Hig]["other"] += eobsbirks;
  } else {
    if( verbosity )
      G4cout << "SensitiveDetectorName: " << SensitiveDetectorName
             << "   localTime2 (ns): " << localTime2/CLHEP::ns
             << "   localTime2SliceLow: " << localSlice2Low
             << "   localTime2SliceMed: " << localSlice2Med
             << "   localTime2SliceHig: " << localSlice2Hig
             << "   particleName: " << particleName
             << "   Eobs_old (low): " << (*Eobs_byParticleAndLocalTime2Low)[SensitiveDetectorName][localSlice2Low][particleName]
             << "   Eobs_old (med): " << (*Eobs_byParticleAndLocalTime2Med)[SensitiveDetectorName][localSlice2Med][particleName]
             << "   Eobs_old (hig): " << (*Eobs_byParticleAndLocalTime2Hig)[SensitiveDetectorName][localSlice2Hig][particleName]
             << "   Eobs_now: " << eobsbirks
             << G4endl;
    (*Eobs_byParticleAndLocalTime2Low)[SensitiveDetectorName][localSlice2Low][particleName] += eobsbirks;
    (*Eobs_byParticleAndLocalTime2Med)[SensitiveDetectorName][localSlice2Med][particleName] += eobsbirks;
    (*Eobs_byParticleAndLocalTime2Hig)[SensitiveDetectorName][localSlice2Hig][particleName] += eobsbirks;
  }
  
  
  //--------------------------------------
  // cerenkov photons by particle and time
  
  particleName = theTrack->GetDefinition()->GetParticleName();
  
  if (NCerenPhotons > 0)
  {
    std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* NCeren_byParticleAndGlobalTimeLow = CaTSEvt->GetNCeren_byParticleAndGlobalTimeLow();
    std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* NCeren_byParticleAndGlobalTimeMed = CaTSEvt->GetNCeren_byParticleAndGlobalTimeMed();
    std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* NCeren_byParticleAndGlobalTimeHig = CaTSEvt->GetNCeren_byParticleAndGlobalTimeHig();
    if (std::find(particleListCeren->begin(),particleListCeren->end(),particleName) == particleListCeren->end()) {
      (*NCeren_byParticleAndGlobalTimeLow)[SensitiveDetectorName][globalSliceLow]["other"] += double(NCerenPhotons);
      (*NCeren_byParticleAndGlobalTimeMed)[SensitiveDetectorName][globalSliceMed]["other"] += double(NCerenPhotons);
      (*NCeren_byParticleAndGlobalTimeHig)[SensitiveDetectorName][globalSliceHig]["other"] += double(NCerenPhotons);
    } else {
      (*NCeren_byParticleAndGlobalTimeLow)[SensitiveDetectorName][globalSliceLow][particleName] += double(NCerenPhotons);
      (*NCeren_byParticleAndGlobalTimeMed)[SensitiveDetectorName][globalSliceMed][particleName] += double(NCerenPhotons);
      (*NCeren_byParticleAndGlobalTimeHig)[SensitiveDetectorName][globalSliceHig][particleName] += double(NCerenPhotons);
    }
    
    std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* NCeren_byParticleAndLocalTime1Low = CaTSEvt->GetNCeren_byParticleAndLocalTime1Low();
    std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* NCeren_byParticleAndLocalTime1Med = CaTSEvt->GetNCeren_byParticleAndLocalTime1Med();
    std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* NCeren_byParticleAndLocalTime1Hig = CaTSEvt->GetNCeren_byParticleAndLocalTime1Hig();
    if (std::find(particleListCeren->begin(),particleListCeren->end(),particleName) == particleListCeren->end()) {
      (*NCeren_byParticleAndLocalTime1Low)[SensitiveDetectorName][localSlice1Low]["other"] += double(NCerenPhotons);
      (*NCeren_byParticleAndLocalTime1Med)[SensitiveDetectorName][localSlice1Med]["other"] += double(NCerenPhotons);
      (*NCeren_byParticleAndLocalTime1Hig)[SensitiveDetectorName][localSlice1Hig]["other"] += double(NCerenPhotons);
    } else {
      (*NCeren_byParticleAndLocalTime1Low)[SensitiveDetectorName][localSlice1Low][particleName] += double(NCerenPhotons);
      (*NCeren_byParticleAndLocalTime1Med)[SensitiveDetectorName][localSlice1Med][particleName] += double(NCerenPhotons);
      (*NCeren_byParticleAndLocalTime1Hig)[SensitiveDetectorName][localSlice1Hig][particleName] += double(NCerenPhotons);
    }
    
    std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* NCeren_byParticleAndLocalTime2Low = CaTSEvt->GetNCeren_byParticleAndLocalTime2Low();
    std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* NCeren_byParticleAndLocalTime2Med = CaTSEvt->GetNCeren_byParticleAndLocalTime2Med();
    std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* NCeren_byParticleAndLocalTime2Hig = CaTSEvt->GetNCeren_byParticleAndLocalTime2Hig();
    if (std::find(particleListCeren->begin(),particleListCeren->end(),particleName) == particleListCeren->end()) {
      (*NCeren_byParticleAndLocalTime2Low)[SensitiveDetectorName][localSlice2Low]["other"] += double(NCerenPhotons);
      (*NCeren_byParticleAndLocalTime2Med)[SensitiveDetectorName][localSlice2Med]["other"] += double(NCerenPhotons);
      (*NCeren_byParticleAndLocalTime2Hig)[SensitiveDetectorName][localSlice2Hig]["other"] += double(NCerenPhotons);
    } else {
      (*NCeren_byParticleAndLocalTime2Low)[SensitiveDetectorName][localSlice2Low][particleName] += double(NCerenPhotons);
      (*NCeren_byParticleAndLocalTime2Med)[SensitiveDetectorName][localSlice2Med][particleName] += double(NCerenPhotons);
      (*NCeren_byParticleAndLocalTime2Hig)[SensitiveDetectorName][localSlice2Hig][particleName] += double(NCerenPhotons);
    }
  }
  
  
  
  //-----------------------
  // INFORMATION BY PROCESS
  //-----------------------
  if( verbosity )
  {
    G4cout << "//-----------------------" << G4endl;
    G4cout << "// INFORMATION BY PROCESS" << G4endl;
    G4cout << "//-----------------------" << G4endl;
  }
  
  const G4VProcess * process = aStep -> GetPostStepPoint() -> GetProcessDefinedStep();
  G4String processName = process -> GetProcessName();
  G4String NameofdecayingParticle = "unknown";
  
  G4SteppingManager* fpSteppingManager = G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager();
  const G4TrackVector * secondary = fpSteppingManager -> GetSecondary();
  G4int totsec = fpSteppingManager -> GetfN2ndariesAlongStepDoIt() +
                 fpSteppingManager -> GetfN2ndariesAtRestDoIt() +
                 fpSteppingManager -> GetfN2ndariesPostStepDoIt();
  
  if (processName == "hadElastic")
  {
    NameofdecayingParticle = fpSteppingManager->GetStep()->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
    if (NameofdecayingParticle != "neutron")
    {
      for (size_t lp = (*secondary).size() - totsec; lp < (*secondary).size(); lp++) {
        G4String pname = (*secondary)[lp] -> GetParticleDefinition() -> GetParticleName();
      }
    }
  }
  
  if (processName == "Decay")
  {
    NameofdecayingParticle = fpSteppingManager->GetStep()->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
    for (size_t lp = (*secondary).size() - totsec; lp < (*secondary).size(); lp++) {
      G4String pname = (*secondary)[lp] -> GetParticleDefinition() -> GetParticleName();
    }
  }
  
  if (!((processName == "Decay")&&(NameofdecayingParticle == "triton"))) // triton decay is just there to get rid of long-living tritons
  {
    std::map<G4String, std::map<G4String, G4int > > * processMult = CaTSEvt -> GetProcessMult();
    processMult->at(SensitiveDetectorName)[processName]++;
    if( verbosity )
      G4cout << "process multiplicity[" << processName << "] = " << processMult->at(SensitiveDetectorName)[processName] << G4endl;
    
    // check if process already occurred:
    std::map<G4String, std::map<G4String, std::map<G4String, products> > >* processMap = CaTSEvt -> GetProcessMap();
    if (processMap->at(SensitiveDetectorName).find(processName) == processMap->at(SensitiveDetectorName).end())
    { // new process
      for (size_t lp = (*secondary).size() - totsec; lp < (*secondary).size(); lp++)
      { // loop over secondaries
        std::map<G4String, products> tmpmap;
        products tmpprod;
        G4String pname = (*secondary)[lp] -> GetParticleDefinition() -> GetParticleName();
        if ((*secondary)[lp] -> GetParticleDefinition() -> GetParticleType() == "nucleus"
            && (*secondary)[lp] -> GetParticleDefinition() -> GetParticleSubType() == "generic") {
          pname = fragment;
        }
        
        // check if particle already occurred:
        if (processMap->at(SensitiveDetectorName)[processName].find(pname) == processMap->at(SensitiveDetectorName)[processName].end()){
          processMap->at(SensitiveDetectorName)[processName][pname].NParticles = 1;
          processMap->at(SensitiveDetectorName)[processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
          processMap->at(SensitiveDetectorName)[processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
        } else {
          processMap->at(SensitiveDetectorName)[processName][pname].NParticles++;
          processMap->at(SensitiveDetectorName)[processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
          processMap->at(SensitiveDetectorName)[processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();
        }
        
        processMap->at(SensitiveDetectorName).insert(std::make_pair(processName, tmpmap));
      } // end loop over secondaries
    }// end new process
    else
    { // process already there
      for (size_t lp = (*secondary).size() - totsec; lp < (*secondary).size(); lp++)
      {
        G4String pname = (*secondary)[lp] -> GetParticleDefinition() -> GetParticleName();
        if ((*secondary)[lp] -> GetParticleDefinition() -> GetParticleType() == "nucleus"
            && (*secondary)[lp] -> GetParticleDefinition() -> GetParticleSubType() == "generic") {
          pname = fragment;
        }
        
        // check if particle already occurred:
        if (processMap->at(SensitiveDetectorName)[processName].find(pname) == processMap->at(SensitiveDetectorName)[processName].end()){
          processMap->at(SensitiveDetectorName)[processName][pname].NParticles = 1;
          processMap->at(SensitiveDetectorName)[processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
          processMap->at(SensitiveDetectorName)[processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
        } else {
          processMap->at(SensitiveDetectorName)[processName][pname].NParticles++;
          processMap->at(SensitiveDetectorName)[processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
          processMap->at(SensitiveDetectorName)[processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();
        }
      }
    } // end of process already there
    
    if( verbosity )
    {
      std::map<G4String, products>::const_iterator mapIt;
      for(mapIt = (*processMap)[SensitiveDetectorName][processName].begin(); mapIt != (*processMap)[SensitiveDetectorName][processName].end(); ++mapIt)
      {
        G4cout << "(*processMap)[" << SensitiveDetectorName << "][" << processName << "][" << mapIt->first << "].NParticles: " << (mapIt->second).NParticles << G4endl;
        G4cout << "(*processMap)[" << SensitiveDetectorName << "][" << processName << "][" << mapIt->first << "].kinE: "       << (mapIt->second).kinE << G4endl;
        G4cout << "(*processMap)[" << SensitiveDetectorName << "][" << processName << "][" << mapIt->first << "].totE: "       << (mapIt->second).totE << G4endl;
      }
    }
  }
  
  
  
  //--------------------------------
  // INFORMATION BY PROCESS AND TIME
  //--------------------------------
  if( verbosity )
  {
    G4cout << "//--------------------------------" << G4endl;
    G4cout << "// INFORMATION BY PROCESS AND TIME" << G4endl;
    G4cout << "//--------------------------------" << G4endl;
  }
  
  if (!((processName == "Decay")&&(NameofdecayingParticle == "triton"))) // triton decay is just there to get rid of long-living tritons
  {
    std::map<G4String, std::map<int, std::map<G4String, G4int > > >* processAndGlobalTimeLowMult = CaTSEvt -> GetProcessAndGlobalTimeLowMult();
    std::map<G4String, std::map<int, std::map<G4String, G4int > > >* processAndGlobalTimeMedMult = CaTSEvt -> GetProcessAndGlobalTimeMedMult();
    std::map<G4String, std::map<int, std::map<G4String, G4int > > >* processAndGlobalTimeHigMult = CaTSEvt -> GetProcessAndGlobalTimeHigMult();
    (*processAndGlobalTimeLowMult)[SensitiveDetectorName][globalSliceLow][processName] += 1;
    (*processAndGlobalTimeMedMult)[SensitiveDetectorName][globalSliceMed][processName] += 1;
    (*processAndGlobalTimeHigMult)[SensitiveDetectorName][globalSliceHig][processName] += 1;
    if( verbosity )
    {
      G4cout << "process multiplicity[" << processName << "][" << globalSliceLow << "] = " << (*processAndGlobalTimeLowMult)[SensitiveDetectorName][globalSliceLow][processName] << G4endl;
      G4cout << "process multiplicity[" << processName << "][" << globalSliceMed << "] = " << (*processAndGlobalTimeMedMult)[SensitiveDetectorName][globalSliceMed][processName] << G4endl;
      G4cout << "process multiplicity[" << processName << "][" << globalSliceHig << "] = " << (*processAndGlobalTimeHigMult)[SensitiveDetectorName][globalSliceHig][processName] << G4endl;
    }
    
    std::map<G4String, std::map<int, std::map<G4String, G4int > > >* processAndLocalTime1LowMult = CaTSEvt -> GetProcessAndLocalTime1LowMult();
    std::map<G4String, std::map<int, std::map<G4String, G4int > > >* processAndLocalTime1MedMult = CaTSEvt -> GetProcessAndLocalTime1MedMult();
    std::map<G4String, std::map<int, std::map<G4String, G4int > > >* processAndLocalTime1HigMult = CaTSEvt -> GetProcessAndLocalTime1HigMult();
    (*processAndLocalTime1LowMult)[SensitiveDetectorName][localSlice1Low][processName] += 1;
    (*processAndLocalTime1MedMult)[SensitiveDetectorName][localSlice1Med][processName] += 1;
    (*processAndLocalTime1HigMult)[SensitiveDetectorName][localSlice1Hig][processName] += 1;
    if( verbosity )
    {
      G4cout << "process multiplicity[" << processName << "][" << localSlice1Low << "] = " << (*processAndLocalTime1LowMult)[SensitiveDetectorName][localSlice1Low][processName] << G4endl;
      G4cout << "process multiplicity[" << processName << "][" << localSlice1Med << "] = " << (*processAndLocalTime1MedMult)[SensitiveDetectorName][localSlice1Med][processName] << G4endl;
      G4cout << "process multiplicity[" << processName << "][" << localSlice1Hig << "] = " << (*processAndLocalTime1HigMult)[SensitiveDetectorName][localSlice1Hig][processName] << G4endl;
    }
    
    std::map<G4String, std::map<int, std::map<G4String, G4int > > >* processAndLocalTime2LowMult = CaTSEvt -> GetProcessAndLocalTime2LowMult();
    std::map<G4String, std::map<int, std::map<G4String, G4int > > >* processAndLocalTime2MedMult = CaTSEvt -> GetProcessAndLocalTime2MedMult();
    std::map<G4String, std::map<int, std::map<G4String, G4int > > >* processAndLocalTime2HigMult = CaTSEvt -> GetProcessAndLocalTime2HigMult();
    (*processAndLocalTime2LowMult)[SensitiveDetectorName][localSlice2Low][processName] += 1;
    (*processAndLocalTime2MedMult)[SensitiveDetectorName][localSlice2Med][processName] += 1;
    (*processAndLocalTime2HigMult)[SensitiveDetectorName][localSlice2Hig][processName] += 1;
    if( verbosity )
    {
      G4cout << "process multiplicity[" << processName << "][" << localSlice2Low << "] = " << (*processAndLocalTime2LowMult)[SensitiveDetectorName][localSlice2Low][processName] << G4endl;
      G4cout << "process multiplicity[" << processName << "][" << localSlice2Med << "] = " << (*processAndLocalTime2MedMult)[SensitiveDetectorName][localSlice2Med][processName] << G4endl;
      G4cout << "process multiplicity[" << processName << "][" << localSlice2Hig << "] = " << (*processAndLocalTime2HigMult)[SensitiveDetectorName][localSlice2Hig][processName] << G4endl;
    }
    
    
    std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* processAndGlobalTimeLowMap = CaTSEvt -> GetProcessAndGlobalTimeLowMap();
    std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* processAndGlobalTimeMedMap = CaTSEvt -> GetProcessAndGlobalTimeMedMap();
    std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* processAndGlobalTimeHigMap = CaTSEvt -> GetProcessAndGlobalTimeHigMap();
    for (size_t lp = (*secondary).size() - totsec; lp < (*secondary).size(); lp++)
    {
      G4String pname = (*secondary)[lp] -> GetParticleDefinition() -> GetParticleName();
      if ((*secondary)[lp] -> GetParticleDefinition() -> GetParticleType() == "nucleus"
          && (*secondary)[lp] -> GetParticleDefinition() -> GetParticleSubType() == "generic") {
        pname = fragment;
      }
      
      if ((*processAndGlobalTimeLowMap)[SensitiveDetectorName][globalSliceLow][processName].find(pname) == (*processAndGlobalTimeLowMap)[SensitiveDetectorName][globalSliceLow][processName].end()) {
        (*processAndGlobalTimeLowMap)[SensitiveDetectorName][globalSliceLow][processName][pname].NParticles = 1;
        (*processAndGlobalTimeLowMap)[SensitiveDetectorName][globalSliceLow][processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
        (*processAndGlobalTimeLowMap)[SensitiveDetectorName][globalSliceLow][processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
      }
      else
      {
        (*processAndGlobalTimeLowMap)[SensitiveDetectorName][globalSliceLow][processName][pname].NParticles++;
        (*processAndGlobalTimeLowMap)[SensitiveDetectorName][globalSliceLow][processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
        (*processAndGlobalTimeLowMap)[SensitiveDetectorName][globalSliceLow][processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();      
      }
      if ((*processAndGlobalTimeMedMap)[SensitiveDetectorName][globalSliceMed][processName].find(pname) == (*processAndGlobalTimeMedMap)[SensitiveDetectorName][globalSliceMed][processName].end()) {
        (*processAndGlobalTimeMedMap)[SensitiveDetectorName][globalSliceMed][processName][pname].NParticles = 1;
        (*processAndGlobalTimeMedMap)[SensitiveDetectorName][globalSliceMed][processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
        (*processAndGlobalTimeMedMap)[SensitiveDetectorName][globalSliceMed][processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
      }
      else
      {
        (*processAndGlobalTimeMedMap)[SensitiveDetectorName][globalSliceMed][processName][pname].NParticles++;
        (*processAndGlobalTimeMedMap)[SensitiveDetectorName][globalSliceMed][processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
        (*processAndGlobalTimeMedMap)[SensitiveDetectorName][globalSliceMed][processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();      
      }
      if ((*processAndGlobalTimeHigMap)[SensitiveDetectorName][globalSliceHig][processName].find(pname) == (*processAndGlobalTimeHigMap)[SensitiveDetectorName][globalSliceHig][processName].end()) {
        (*processAndGlobalTimeHigMap)[SensitiveDetectorName][globalSliceHig][processName][pname].NParticles = 1;
        (*processAndGlobalTimeHigMap)[SensitiveDetectorName][globalSliceHig][processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
        (*processAndGlobalTimeHigMap)[SensitiveDetectorName][globalSliceHig][processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
      }
      else
      {
        (*processAndGlobalTimeHigMap)[SensitiveDetectorName][globalSliceHig][processName][pname].NParticles++;
        (*processAndGlobalTimeHigMap)[SensitiveDetectorName][globalSliceHig][processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
        (*processAndGlobalTimeHigMap)[SensitiveDetectorName][globalSliceHig][processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();      
      }
    }
    
    std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* processAndLocalTime1LowMap = CaTSEvt -> GetProcessAndLocalTime1LowMap();
    std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* processAndLocalTime1MedMap = CaTSEvt -> GetProcessAndLocalTime1MedMap();
    std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* processAndLocalTime1HigMap = CaTSEvt -> GetProcessAndLocalTime1HigMap();
    for (size_t lp = (*secondary).size() - totsec; lp < (*secondary).size(); lp++)
    {
      G4String pname = (*secondary)[lp] -> GetParticleDefinition() -> GetParticleName();
      if ((*secondary)[lp] -> GetParticleDefinition() -> GetParticleType() == "nucleus"
          && (*secondary)[lp] -> GetParticleDefinition() -> GetParticleSubType() == "generic") {
        pname = fragment;
      }
      
      if ((*processAndLocalTime1LowMap)[SensitiveDetectorName][localSlice1Low][processName].find(pname) == (*processAndLocalTime1LowMap)[SensitiveDetectorName][localSlice1Low][processName].end()) {
        (*processAndLocalTime1LowMap)[SensitiveDetectorName][localSlice1Low][processName][pname].NParticles = 1;
        (*processAndLocalTime1LowMap)[SensitiveDetectorName][localSlice1Low][processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
        (*processAndLocalTime1LowMap)[SensitiveDetectorName][localSlice1Low][processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
      }
      else
      {
        (*processAndLocalTime1LowMap)[SensitiveDetectorName][localSlice1Low][processName][pname].NParticles++;
        (*processAndLocalTime1LowMap)[SensitiveDetectorName][localSlice1Low][processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
        (*processAndLocalTime1LowMap)[SensitiveDetectorName][localSlice1Low][processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();      
      }
      if ((*processAndLocalTime1MedMap)[SensitiveDetectorName][localSlice1Med][processName].find(pname) == (*processAndLocalTime1MedMap)[SensitiveDetectorName][localSlice1Med][processName].end()) {
        (*processAndLocalTime1MedMap)[SensitiveDetectorName][localSlice1Med][processName][pname].NParticles = 1;
        (*processAndLocalTime1MedMap)[SensitiveDetectorName][localSlice1Med][processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
        (*processAndLocalTime1MedMap)[SensitiveDetectorName][localSlice1Med][processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
      }
      else
      {
        (*processAndLocalTime1MedMap)[SensitiveDetectorName][localSlice1Med][processName][pname].NParticles++;
        (*processAndLocalTime1MedMap)[SensitiveDetectorName][localSlice1Med][processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
        (*processAndLocalTime1MedMap)[SensitiveDetectorName][localSlice1Med][processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();      
      }
      if ((*processAndLocalTime1HigMap)[SensitiveDetectorName][localSlice1Hig][processName].find(pname) == (*processAndLocalTime1HigMap)[SensitiveDetectorName][localSlice1Hig][processName].end()) {
        (*processAndLocalTime1HigMap)[SensitiveDetectorName][localSlice1Hig][processName][pname].NParticles = 1;
        (*processAndLocalTime1HigMap)[SensitiveDetectorName][localSlice1Hig][processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
        (*processAndLocalTime1HigMap)[SensitiveDetectorName][localSlice1Hig][processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
      }
      else
      {
        (*processAndLocalTime1HigMap)[SensitiveDetectorName][localSlice1Hig][processName][pname].NParticles++;
        (*processAndLocalTime1HigMap)[SensitiveDetectorName][localSlice1Hig][processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
        (*processAndLocalTime1HigMap)[SensitiveDetectorName][localSlice1Hig][processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();      
      }
    }
    
    std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* processAndLocalTime2LowMap = CaTSEvt -> GetProcessAndLocalTime2LowMap();
    std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* processAndLocalTime2MedMap = CaTSEvt -> GetProcessAndLocalTime2MedMap();
    std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* processAndLocalTime2HigMap = CaTSEvt -> GetProcessAndLocalTime2HigMap();
    for (size_t lp = (*secondary).size() - totsec; lp < (*secondary).size(); lp++)
    {
      G4String pname = (*secondary)[lp] -> GetParticleDefinition() -> GetParticleName();
      if ((*secondary)[lp] -> GetParticleDefinition() -> GetParticleType() == "nucleus"
          && (*secondary)[lp] -> GetParticleDefinition() -> GetParticleSubType() == "generic") {
        pname = fragment;
      }
      
      if ((*processAndLocalTime2LowMap)[SensitiveDetectorName][localSlice2Low][processName].find(pname) == (*processAndLocalTime2LowMap)[SensitiveDetectorName][localSlice2Low][processName].end()) {
        (*processAndLocalTime2LowMap)[SensitiveDetectorName][localSlice2Low][processName][pname].NParticles = 1;
        (*processAndLocalTime2LowMap)[SensitiveDetectorName][localSlice2Low][processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
        (*processAndLocalTime2LowMap)[SensitiveDetectorName][localSlice2Low][processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
      }
      else
      {
        (*processAndLocalTime2LowMap)[SensitiveDetectorName][localSlice2Low][processName][pname].NParticles++;
        (*processAndLocalTime2LowMap)[SensitiveDetectorName][localSlice2Low][processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
        (*processAndLocalTime2LowMap)[SensitiveDetectorName][localSlice2Low][processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();      
      }
      if ((*processAndLocalTime2MedMap)[SensitiveDetectorName][localSlice2Med][processName].find(pname) == (*processAndLocalTime2MedMap)[SensitiveDetectorName][localSlice2Med][processName].end()) {
        (*processAndLocalTime2MedMap)[SensitiveDetectorName][localSlice2Med][processName][pname].NParticles = 1;
        (*processAndLocalTime2MedMap)[SensitiveDetectorName][localSlice2Med][processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
        (*processAndLocalTime2MedMap)[SensitiveDetectorName][localSlice2Med][processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
      }
      else
      {
        (*processAndLocalTime2MedMap)[SensitiveDetectorName][localSlice2Med][processName][pname].NParticles++;
        (*processAndLocalTime2MedMap)[SensitiveDetectorName][localSlice2Med][processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
        (*processAndLocalTime2MedMap)[SensitiveDetectorName][localSlice2Med][processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();      
      }
      if ((*processAndLocalTime2HigMap)[SensitiveDetectorName][localSlice2Hig][processName].find(pname) == (*processAndLocalTime2HigMap)[SensitiveDetectorName][localSlice2Hig][processName].end()) {
        (*processAndLocalTime2HigMap)[SensitiveDetectorName][localSlice2Hig][processName][pname].NParticles = 1;
        (*processAndLocalTime2HigMap)[SensitiveDetectorName][localSlice2Hig][processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
        (*processAndLocalTime2HigMap)[SensitiveDetectorName][localSlice2Hig][processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
      }
      else
      {
        (*processAndLocalTime2HigMap)[SensitiveDetectorName][localSlice2Hig][processName][pname].NParticles++;
        (*processAndLocalTime2HigMap)[SensitiveDetectorName][localSlice2Hig][processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
        (*processAndLocalTime2HigMap)[SensitiveDetectorName][localSlice2Hig][processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();      
      }
    }
    
    if( verbosity )
    {
      std::map<G4String, products>::const_iterator mapIt;
      
      for(mapIt = (*processAndGlobalTimeLowMap)[SensitiveDetectorName][globalSliceLow][processName].begin(); mapIt != (*processAndGlobalTimeLowMap)[SensitiveDetectorName][globalSliceLow][processName].end(); ++mapIt)
      {
        G4cout << "(*processAndGlobalTimeLowMap)[" << SensitiveDetectorName << "][" << globalSliceLow << "][" << processName << "][" << mapIt->first << "].NParticles: " << (mapIt->second).NParticles << G4endl;
        G4cout << "(*processAndGlobalTimeLowMap)[" << SensitiveDetectorName << "][" << globalSliceLow << "][" << processName << "][" << mapIt->first << "].kinE: "       << (mapIt->second).kinE << G4endl;
        G4cout << "(*processAndGlobalTimeLowMap)[" << SensitiveDetectorName << "][" << globalSliceLow << "][" << processName << "][" << mapIt->first << "].totE: "       << (mapIt->second).totE << G4endl;
      }
      for(mapIt = (*processAndGlobalTimeMedMap)[SensitiveDetectorName][globalSliceMed][processName].begin(); mapIt != (*processAndGlobalTimeMedMap)[SensitiveDetectorName][globalSliceMed][processName].end(); ++mapIt)
      {
        G4cout << "(*processAndGlobalTimeMedMap)[" << SensitiveDetectorName << "][" << globalSliceMed << "][" << processName << "][" << mapIt->first << "].NParticles: " << (mapIt->second).NParticles << G4endl;
        G4cout << "(*processAndGlobalTimeMedMap)[" << SensitiveDetectorName << "][" << globalSliceMed << "][" << processName << "][" << mapIt->first << "].kinE: "       << (mapIt->second).kinE << G4endl;
        G4cout << "(*processAndGlobalTimeMedMap)[" << SensitiveDetectorName << "][" << globalSliceMed << "][" << processName << "][" << mapIt->first << "].totE: "       << (mapIt->second).totE << G4endl;
      }
      for(mapIt = (*processAndGlobalTimeHigMap)[SensitiveDetectorName][globalSliceHig][processName].begin(); mapIt != (*processAndGlobalTimeHigMap)[SensitiveDetectorName][globalSliceHig][processName].end(); ++mapIt)
      {
        G4cout << "(*processAndGlobalTimeHigMap)[" << SensitiveDetectorName << "][" << globalSliceHig << "][" << processName << "][" << mapIt->first << "].NParticles: " << (mapIt->second).NParticles << G4endl;
        G4cout << "(*processAndGlobalTimeHigMap)[" << SensitiveDetectorName << "][" << globalSliceHig << "][" << processName << "][" << mapIt->first << "].kinE: "       << (mapIt->second).kinE << G4endl;
        G4cout << "(*processAndGlobalTimeHigMap)[" << SensitiveDetectorName << "][" << globalSliceHig << "][" << processName << "][" << mapIt->first << "].totE: "       << (mapIt->second).totE << G4endl;
      }
    }
  }
  
  
  
  //--------------------
  // INFORMATION BY CELL
  //--------------------
  if( verbosity )
  {
    G4cout << "//--------------------" << G4endl;
    G4cout << "// INFORMATION BY CELL" << G4endl;
    G4cout << "//--------------------" << G4endl;
  }
 
  partName = particleName;
  if ( theTrack->GetParticleDefinition()->GetParticleType() == "nucleus" ) partName = fragment;
  if ( theTrack->GetParticleDefinition()->GetParticleType() == "meson" )   partName = "meson";
  if ( theTrack->GetParticleDefinition()->GetParticleType() == "baryon" )  partName = "baryon";
  if ( particleName == "proton"  || particleName == "anti_proton" )        partName = "proton";
  if ( particleName == "neutron" || particleName == "anti_neutron" )       partName = "neutron";
  if ( particleName == "e-"      || particleName == "e+" )                 partName = "e";
  if ( particleName == "mu-"     || particleName == "mu+" )                partName = "mu";
  if ( particleName == "pi+"     || particleName == "pi-" )                partName = "pi";
  if( aTrackInfo->GetParticleIsEM() )                                      partName = "em";
  
  if ( std::find(particleListShort->begin(),particleListShort->end(),partName) == particleList->end() )
    partName = "other";
  
  //-------------------------------------------
  // this cell already exists in hit collection
  
  G4bool hitExists = false;
  for (G4int j = 0; j < drtscalorimeterCollection2->entries(); j++)
  {
    DRTSCalorimeterHit2* aPreviousHit = (*drtscalorimeterCollection2)[j];
    if (cellPosition == aPreviousHit->GetPos() &&
        partName == aPreviousHit->GetParticleName() && 
        globalSliceLow == aPreviousHit->GetGlobalSliceLow() &&
        globalSliceMed == aPreviousHit->GetGlobalSliceMed() &&
        globalSliceHig == aPreviousHit->GetGlobalSliceHig() )
    {
      aPreviousHit->SetEdep(edep + aPreviousHit->GetEdep());
      aPreviousHit->SetEobsbirks(eobsbirks + aPreviousHit->GetEobsbirks());
      aPreviousHit->SetNCeren(NCerenPhotons + aPreviousHit->GetNCeren());
      
      if( verbosity )
        G4cout << "hit:   cellPosition: "   << cellPosition.x() << "," << cellPosition.y() << "," << cellPosition.z()
               << "   Edep: "               << aPreviousHit->GetEdep()
               << "   Eobs: "               << aPreviousHit->GetEobsbirks()
               << "   globalTimeSliceLow: " << globalSliceMed
               << "   globalTimeSliceMed: " << globalSliceMed
               << "   globalTimeSliceHig: " << globalSliceMed
               << G4endl;
      
      if( verbosity )
        G4cout << G4endl;
      
      hitExists = true;
      break;
    }
  }
  
  if( hitExists == false )
  {
    DRTSCalorimeterHit2* newHit = new DRTSCalorimeterHit2(partName,edep,eobsbirks,NCerenPhotons,cellPosition,globalSliceLow,globalSliceMed,globalSliceHig);
    if( verbosity )
      G4cout << "new hit: "           << drtscalorimeterCollection2->GetSize()
             << "new hit: "           << drtscalorimeterCollection2->entries()
             << "   cellPosition: "   << cellPosition.x() << "," << cellPosition.y() << "," << cellPosition.z()
             << "   Edep: "           << edep
             << "   Eobs: "           << eobsbirks
             << "   globalSliceLow: " << globalSliceLow
             << "   globalSliceMed: " << globalSliceMed
             << "   globalSliceHig: " << globalSliceHig
             << G4endl;
    drtscalorimeterCollection2->insert(newHit);
  }
  
  
  
  //--------------------
  // INFORMATION BY CELL
  //--------------------
  if( verbosity )
  {
    G4cout << "//--------------------" << G4endl;
    G4cout << "// INFORMATION BY CELL" << G4endl;
    G4cout << "//--------------------" << G4endl;
  }
  
  //-------------------------------------------
  // this cell already exists in hit collection
  
  for (G4int j = 0; j < drtscalorimeterCollection->entries(); j++)
  {
    DRTSCalorimeterHit* aPreviousHit = (*drtscalorimeterCollection)[j];
    if (cellPosition == aPreviousHit->GetPos())
    {
      aPreviousHit->SetEdep(edep + aPreviousHit->GetEdep());
      aPreviousHit->SetEobsbirks(eobsbirks + aPreviousHit->GetEobsbirks());
      aPreviousHit->SetNCeren(NCerenPhotons + aPreviousHit->GetNCeren());
      if ((particleName == "e+") || (particleName == "gamma") || (particleName == "e-")) {
        aPreviousHit->SetEdepEM(edep + aPreviousHit->GetEdepEM());
      } else {
        aPreviousHit->SetEdepnonEM(edep + aPreviousHit->GetEdepnonEM());
      }
      if (globalTime < aPreviousHit->GetGlobalTime()) {
        aPreviousHit->SetGlobalTime(globalTime/CLHEP::ns);
      }
      
      if( verbosity )
        G4cout << "hit:   cellPosition: " << cellPosition.x() << "," << cellPosition.y() << "," << cellPosition.z()
               << "   Edep(em): "     << aPreviousHit->GetEdepEM()
               << "   Edep(non-em): " << aPreviousHit->GetEdepnonEM()
               << "   Eobs: "         << aPreviousHit->GetEobsbirks()
               << "   globalTimeSliceMed: " << int(globalTime/timeslicemed)
               << "   localTimeSlice1Med: " << int(localTime1/timeslicemed)
               << "   localTimeSlice2Med: " << int(localTime2/timeslicemed)
               << G4endl;
      
      bool foundslice = false;
      std::vector<DRTimeSliceHit>* tslhtcol = aPreviousHit->GetDRGlobalTimeSliceHitCol();
      for (unsigned int jj = 0; jj < tslhtcol->size(); jj++)
      {
        DRTimeSliceHit tslht = tslhtcol->at(jj);
        if (globalSliceMed == tslht.GetSlice()) // TimeSliceHit for this time slice exists
        {
          foundslice = true;
          if ((particleName == "e+") || (particleName == "gamma") || (particleName == "e-")) {
            tslht.SetEdepEM(edep + tslht.GetEdepEM());
          } else {
            tslht.SetEdepnonEM(edep + tslht.GetEdepnonEM());
          }
          tslht.SetEdep(edep + tslht.GetEdep());
          std::swap(tslhtcol->at(jj), tslht);
          break;
        }
      }
      if (!foundslice) {
        if ((particleName == "e+") || (particleName == "gamma") || (particleName == "e-")) {
          aPreviousHit->AddDRTimeSliceHit(0, globalSliceMed, edep, edep, 0.0, NCerenPhotons);
        } else {
          aPreviousHit->AddDRTimeSliceHit(0, globalSliceMed, edep, 0.0, edep, NCerenPhotons);
        }
      }
      
      foundslice = false;
      tslhtcol = aPreviousHit->GetDRLocalTime1SliceHitCol();
      for (unsigned int jj = 0; jj < tslhtcol->size(); jj++)
      {
        DRTimeSliceHit tslht = tslhtcol->at(jj);
        if (localSlice1Med == tslht.GetSlice()) // TimeSliceHit for this time slice exists
        {
          foundslice = true;
          if ((particleName == "e+") || (particleName == "gamma") || (particleName == "e-")) {
            tslht.SetEdepEM(edep + tslht.GetEdepEM());
          } else {
            tslht.SetEdepnonEM(edep + tslht.GetEdepnonEM());
          }
          tslht.SetEdep(edep + tslht.GetEdep());
          std::swap(tslhtcol->at(jj), tslht);
          break;
        }
      }
      if (!foundslice) {
        if ((particleName == "e+") || (particleName == "gamma") || (particleName == "e-")) {
          aPreviousHit->AddDRTimeSliceHit(1, localSlice1Med, edep, edep, 0.0, NCerenPhotons);
        } else {
          aPreviousHit->AddDRTimeSliceHit(1, localSlice1Med, edep, 0.0, edep, NCerenPhotons);
        }
      }
      
      foundslice = false;
      tslhtcol = aPreviousHit->GetDRLocalTime1SliceHitCol();
      for (unsigned int jj = 0; jj < tslhtcol->size(); jj++)
      {
        DRTimeSliceHit tslht = tslhtcol->at(jj);
        if (localSlice2Med == tslht.GetSlice()) // TimeSliceHit for this time slice exists
        {
          foundslice = true;
          if ((particleName == "e+") || (particleName == "gamma") || (particleName == "e-")) {
            tslht.SetEdepEM(edep + tslht.GetEdepEM());
          } else {
            tslht.SetEdepnonEM(edep + tslht.GetEdepnonEM());
          }
          tslht.SetEdep(edep + tslht.GetEdep());
          std::swap(tslhtcol->at(jj), tslht);
          break;
        }
      }
      if (!foundslice) {
        if ((particleName == "e+") || (particleName == "gamma") || (particleName == "e-")) {
          aPreviousHit->AddDRTimeSliceHit(2, localSlice2Med, edep, edep, 0.0, NCerenPhotons);
        } else {
          aPreviousHit->AddDRTimeSliceHit(2, localSlice2Med, edep, 0.0, edep, NCerenPhotons);
        }
      }

      if( verbosity )
        G4cout << G4endl;
      
      return true;
    }
  }
  
  //------------------------------------------
  // this cell doesn't exist in hit collection
  
  DRTSCalorimeterHit* newHit;
  if ((particleName == "e+") || (particleName == "gamma") || (particleName == "e-")) {
    newHit = new DRTSCalorimeterHit(edep, edep, 0.0, eobsbirks, NCerenPhotons, cellPosition, globalTime, localTime1, localTime2, timeslicemed, dedx, vel);
    if( verbosity )
      G4cout << "new hit:   cellPosition: " << cellPosition.x() << "," << cellPosition.y() << "," << cellPosition.z()
             << "   Edep(em): "     << edep
             << "   Edep(non-em): " << 0.
             << "   Eobs: "         << eobsbirks
             << "   globalTimeSliceMed: " << globalSliceMed
             << "   localTimeSlice1Med: " << localSlice1Med
             << "   localTimeSlice2Med: " << localSlice1Med
             << G4endl;
  } else {
    newHit = new DRTSCalorimeterHit(edep, 0.0, edep, eobsbirks, NCerenPhotons, cellPosition, globalTime, localTime1, localTime2, timeslicemed, dedx, vel);
    if( verbosity )
      G4cout << "new hit:   cellPosition: " << cellPosition.x() << "," << cellPosition.y() << "," << cellPosition.z()
             << "   Edep(em): "     << 0.
             << "   Edep(non-em): " << edep
             << "   Eobs: "         << eobsbirks
             << "   globalTimeSliceMed: " << globalSliceMed
             << "   localTimeSlice1Med: " << localSlice1Med
             << "   localTimeSlice2Med: " << localSlice2Med
             << G4endl;
  }
  drtscalorimeterCollection->insert(newHit);
  
  
  
  if( verbosity )
    G4cout << G4endl;
    
  return true;
}
