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
  G4VSensitiveDetector(name)
{
#ifdef G4ANALYSIS_USE
  analysis = Analysis::getInstance();
#endif
  
  pMessenger = new DRTSCalorimeterSDMessenger(this);
  birksc1 = 1.29e-2;
  birksc2 = 9.59e-6;
  CerenGenerator = new(Cerenkov);
  
  timeSliceTypes.push_back("Low");
  timeSliceTypes.push_back("Med");
  timeSliceTypes.push_back("Hig");
  
  particleList = NULL;
  processList = NULL;
  
  instance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



DRTSCalorimeterSD::~DRTSCalorimeterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void DRTSCalorimeterSD::Initialize(G4HCofThisEvent*)
{
  if( particleList == NULL )
  {
    RunAction* RA = RunAction::getInstance();
    
    timeSliceSizes = RA -> GetTimeSliceSizes();
    minTimes       = RA -> GetMinTimes();
    maxTimes       = RA -> GetMaxTimes();
    
    particleList = RA -> GetParticleList();
    processList  = RA -> GetProcessList();
  }
  
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
  G4float offsetPosition = 0.5 * fabs( solid->DistanceToOut(localPosition,zAxisNeg)+solid->DistanceToOut(localPosition,zAxisPos) );
  
  const G4float globalTime = aStep->GetPreStepPoint()->GetGlobalTime();
  const G4float localTime1 = globalTime - cellPosition.z()/c_light + offsetPosition/c_light;
  const G4float localTime2 = globalTime - sqrt(pow(cellPosition.x(),2)+pow(cellPosition.y(),2)+pow(cellPosition.z(),2))/c_light  + offsetPosition/c_light;
  
  times["globalTime"] = globalTime;
  times["localTime1"] = localTime1;
  //times["localTime2"] = localTime2;
  
  globalPositionX -> Fill(globalPosition.x());
  globalPositionY -> Fill(globalPosition.y());
  globalPositionZ -> Fill(globalPosition.z());
  localPositionX -> Fill(localPosition.x());
  localPositionY -> Fill(localPosition.y());
  localPositionZ -> Fill(localPosition.z());
  
  
  
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
  
  G4float edep = aStep->GetTotalEnergyDeposit() / CLHEP::GeV;
  if (edep == 0.) return false;
  G4float rho = aStep->GetPreStepPoint()->GetMaterial()->GetDensity() / (CLHEP::g / CLHEP::cm3);
  G4float stepLength = aStep->GetStepLength() / CLHEP::cm;
  G4float birksFactor = 1.0;
  G4float dedx;
  
  //Do not apply Birks for gamma deposits!
  if (stepLength > 1.0e-8)//Check, cut if necessary.
  {
    dedx = edep/CLHEP::MeV / (rho * stepLength); //[MeV*cm^2/g]
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
  G4float eobs = edep*birksFactor;
  
  G4Track* theTrack = aStep->GetTrack();
  G4int NCerenPhotons = 0;
  G4String particleName = GetParticleName(theTrack);
  if( std::find(particleList->begin(),particleList->end(),particleName) == particleList->end() ) particleList->push_back(particleName);
  G4String processName = (aStep->GetPostStepPoint()->GetProcessDefinedStep())->GetProcessName();
  if( std::find(processList->begin(),processList->end(),processName) == processList->end() ) processList->push_back(processName);
  G4StepPoint* pPreStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep->GetPostStepPoint();
  const G4float charge = theTrack->GetDefinition()->GetPDGCharge();
  G4float beta = 0.5 * (pPreStepPoint->GetBeta() + pPostStepPoint->GetBeta());
  G4String thematerial = touch->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName();
  G4float MeanNumberOfPhotons = CerenGenerator->GetAverageNumberOfPhotons(charge, beta, thematerial);
  if (MeanNumberOfPhotons > 0.0) {
    G4float step_length = aStep->GetStepLength();
    MeanNumberOfPhotons = MeanNumberOfPhotons * step_length;
    NCerenPhotons = (G4int) G4Poisson(MeanNumberOfPhotons);
  } else {
    NCerenPhotons = 0;
  }

  if( verbosity )
  {
    TrackInformation* aTrackInfo = (TrackInformation*)( aStep->GetTrack()->GetUserInformation() );
    G4cout << "isEM: " << aTrackInfo->GetParticleIsEM() << G4endl;
  }
  
  
  
  
  

  //-------------------
  // GLOBAL INFORMATION
  //-------------------
  if( verbosity )
  {
    G4cout << "//------------------------" << G4endl;
    G4cout << "// GLOBAL INFORMATION"      << G4endl;
    G4cout << "//------------------------" << G4endl;
  }
  
  
  //--------------------
  // global event energy
  
  CaTSEvt->SetTotDepEnergy(CaTSEvt->GetTotDepEnergy() + edep);
  CaTSEvt->SetTotObsEnergy(CaTSEvt->GetTotObsEnergy() + eobs);
  CaTSEvt->SetTotNCeren(CaTSEvt->GetNCeren() + G4float(NCerenPhotons));
  
  
  //----------------------------
  // event energy per sub-volume
  
  std::map<G4String,G4float>* m_Edep = CaTSEvt->GetEdepMap();
  std::map<G4String,G4float>* m_Eobs = CaTSEvt->GetEobsMap();
  std::map<G4String,G4float>* m_NCeren = CaTSEvt->GetNCerenMap();
  (*m_Edep)[SensitiveDetectorName] += edep;
  (*m_Eobs)[SensitiveDetectorName] += eobs;
  (*m_NCeren)[SensitiveDetectorName] += NCerenPhotons;
  
  
  //-------------------------------------------------
  // event energy per sub-volume and particle/process
  
  std::map<G4String, std::map<G4String, G4float> >* m_Edep_byParticle = CaTSEvt->GetEdepByParticleMap();
  std::map<G4String, std::map<G4String, G4float> >* m_Eobs_byParticle = CaTSEvt->GetEobsByParticleMap();
  std::map<G4String, std::map<G4String, G4float > >* m_NCeren_byParticle = CaTSEvt->GetNCerenByParticleMap();
  (*m_Edep_byParticle)[SensitiveDetectorName][particleName] += edep;
  (*m_Eobs_byParticle)[SensitiveDetectorName][particleName] += eobs;
  (*m_NCeren_byParticle)[SensitiveDetectorName][particleName] += G4float(NCerenPhotons);
  
  std::map<G4String, std::map<G4String, G4float> >* m_Edep_byProcess = CaTSEvt->GetEdepByProcessMap();
  std::map<G4String, std::map<G4String, G4float> >* m_Eobs_byProcess = CaTSEvt->GetEobsByProcessMap();
  std::map<G4String, std::map<G4String, G4float > >* m_NCeren_byProcess = CaTSEvt->GetNCerenByProcessMap();
  (*m_Edep_byProcess)[SensitiveDetectorName][processName] += edep;
  (*m_Eobs_byProcess)[SensitiveDetectorName][processName] += eobs;
  (*m_NCeren_byProcess)[SensitiveDetectorName][processName] += G4float(NCerenPhotons);
  
  if( verbosity )
  {
    G4cout << ">>> SensitiveVolume: " << SensitiveDetectorName
           << "   particleName: " << particleName
           << "   processName: " << processName
           << "   Edep: " << edep/CLHEP::MeV
           << G4endl;
  }
  
  
  // //-------------------------------------------------------
  // // event energy per sub-volume and detailed process infos
  
  // const G4VProcess * process = aStep -> GetPostStepPoint() -> GetProcessDefinedStep();
  // G4String processName = process -> GetProcessName();
  // if( std::find(processList->begin(),processList->end(),processName) == processList->end() )processList->push_back(processName);
  
  // G4SteppingManager* fpSteppingManager = G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager();
  // G4String NameofdecayingParticle = "unknown";
  // if( processName == "hadElastic" || processName == "Decay" )
  //   NameofdecayingParticle = fpSteppingManager->GetStep()->GetTrack()->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
  
  // const G4TrackVector* secondary = fpSteppingManager -> GetSecondary();
  // G4int totsec = fpSteppingManager -> GetfN2ndariesAlongStepDoIt() +
  //                fpSteppingManager -> GetfN2ndariesAtRestDoIt() +
  //                fpSteppingManager -> GetfN2ndariesPostStepDoIt();
  
  // std::map<G4String,std::map<G4String,G4int > >* processMult = CaTSEvt -> GetProcessMult();
  // std::map<G4String,std::map<G4String,std::map<G4String, products> > >* processMap = CaTSEvt -> GetProcessMap();
  
  // if( !(processName == "Decay" && NameofdecayingParticle == "triton") ) // triton decay is just there to get rid of long-living tritons
  // {
  //   (*processMult)[SensitiveDetectorName][processName] += 1;
    
  //   for(size_t lp = (*secondary).size()-totsec; lp < (*secondary).size(); ++lp) // loop over secondaries
  //   {
  //     G4String pname = (*secondary)[lp] -> GetParticleDefinition() -> GetParticleName();
  //     if ((*secondary)[lp] -> GetParticleDefinition() -> GetParticleType() == "nucleus"
  //         && (*secondary)[lp] -> GetParticleDefinition() -> GetParticleSubType() == "generic") {
  //       pname = "fragment";
  //     }
      
  //     // check if particle already occurred:
  //     if( (*processMap)[SensitiveDetectorName][processName].find(pname) == (*processMap)[SensitiveDetectorName][processName].end() ){
  //       (*processMap)[SensitiveDetectorName][processName][pname].NParticles = 1;
  //       (*processMap)[SensitiveDetectorName][processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
  //       (*processMap)[SensitiveDetectorName][processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
  //     } else {
  //       (*processMap)[SensitiveDetectorName][processName][pname].NParticles++;
  //       (*processMap)[SensitiveDetectorName][processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
  //       (*processMap)[SensitiveDetectorName][processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();
  //     }
  //   } // end loop over secondaries
  // }
  
  // if( verbosity )
  // {
  //   G4cout << ">>> SensitiveVolume: " << SensitiveDetectorName
  //          << "   processName: " << processName
  //          << "   Edep: " << edep/CLHEP::MeV
  //          << G4endl;
    
  //   G4cout << ">>> loop over secondaries <<<" << G4endl;
  //   for(size_t lp = (*secondary).size()-totsec; lp < (*secondary).size(); ++lp) // loop over secondaries
  //   {
  //     G4String pname = (*secondary)[lp] -> GetParticleDefinition() -> GetParticleName();
  //     if ((*secondary)[lp] -> GetParticleDefinition() -> GetParticleType() == "nucleus"
  //         && (*secondary)[lp] -> GetParticleDefinition() -> GetParticleSubType() == "generic") {
  //       pname = "fragment";
  //     }
      
  //     G4cout << ">>>>>> processName: " << processName
  //            << "   pname: " << pname
  //            << "   num: " << (*processMap)[SensitiveDetectorName][processName][pname].NParticles
  //            << "   totE: " << (*processMap)[SensitiveDetectorName][processName][pname].totE/CLHEP::MeV
  //            << G4endl;
  //   }
  // }
  
  
  
  //--------------------
  // INFORMATION BY TIME
  //--------------------
  if( verbosity )
  {
    G4cout << "//------------------------" << G4endl;
    G4cout << "// INFORMATION BY TIME"     << G4endl;
    G4cout << "//------------------------" << G4endl;
  }
  
  for(std::map<G4String,G4float>::const_iterator mapIt = times.begin(); mapIt != times.end(); ++mapIt)
  {
    G4String timeType = mapIt -> first;
    
    for(unsigned int timeSliceTypeIt = 0; timeSliceTypeIt < timeSliceTypes.size(); ++timeSliceTypeIt)
    {
      G4String timeSliceType = timeSliceTypes.at(timeSliceTypeIt);
      G4String timeSliceName = timeType+timeSliceType;
      
      G4int timeSlice = int(mapIt->second/timeSliceSizes[timeSliceType]) + 1;
      if( (mapIt->second-minTimes[timeSliceType]) <  minTimes[timeSliceType] ) timeSlice = 0;
      if( (mapIt->second-minTimes[timeSliceType]) >= maxTimes[timeSliceType] ) timeSlice = int((maxTimes[timeSliceType]-minTimes[timeSliceType])/timeSliceSizes[timeSliceType]) + 1;
      
      
      //----------------------------
      // event energy per sub-volume
      
      std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* m_Edep_byTime = CaTSEvt->GetEdepByTimeMap();
      std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* m_Eobs_byTime = CaTSEvt->GetEobsByTimeMap();
      std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* m_NCeren_byTime = CaTSEvt->GetNCerenByTimeMap();
      (*m_Edep_byTime)[SensitiveDetectorName][timeSliceName][timeSlice] += edep;
      (*m_Eobs_byTime)[SensitiveDetectorName][timeSliceName][timeSlice] += eobs;
      (*m_NCeren_byTime)[SensitiveDetectorName][timeSliceName][timeSlice] += NCerenPhotons;
      
      
      //-------------------------------------------------
      // event energy per sub-volume and particle/process
      
      std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String, G4float> > > >* m_Edep_byParticleAndTime = CaTSEvt->GetEdepByParticleAndTimeMap();
      std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String, G4float> > > >* m_Eobs_byParticleAndTime = CaTSEvt->GetEobsByParticleAndTimeMap();
      std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String, G4float> > > >* m_NCeren_byParticleAndTime = CaTSEvt->GetNCerenByParticleAndTimeMap();
      (*m_Edep_byParticleAndTime)[SensitiveDetectorName][timeSliceName][timeSlice][particleName] += edep;
      (*m_Eobs_byParticleAndTime)[SensitiveDetectorName][timeSliceName][timeSlice][particleName] += eobs;
      (*m_NCeren_byParticleAndTime)[SensitiveDetectorName][timeSliceName][timeSlice][particleName] += G4float(NCerenPhotons);
      
      std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String, G4float> > > >* m_Edep_byProcessAndTime = CaTSEvt->GetEdepByProcessAndTimeMap();
      std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String, G4float> > > >* m_Eobs_byProcessAndTime = CaTSEvt->GetEobsByProcessAndTimeMap();
      std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String, G4float> > > >* m_NCeren_byProcessAndTime = CaTSEvt->GetNCerenByProcessAndTimeMap();
      (*m_Edep_byProcessAndTime)[SensitiveDetectorName][timeSliceName][timeSlice][processName] += edep;
      (*m_Eobs_byProcessAndTime)[SensitiveDetectorName][timeSliceName][timeSlice][processName] += eobs;
      (*m_NCeren_byProcessAndTime)[SensitiveDetectorName][timeSliceName][timeSlice][processName] += G4float(NCerenPhotons);
      
      
      // //------------------------------------------------------
      // // event energy per sub-volume and detailed process info
      
      // std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4int > > > >* processAndTimeMult = CaTSEvt -> GetProcessAndTimeMult();
      // std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,std::map<G4String, products> > > > >* processAndTimeMap = CaTSEvt -> GetProcessAndTimeMap();
      
      // if( !(processName == "Decay" && NameofdecayingParticle == "triton") ) // triton decay is just there to get rid of long-living tritons
      // {
      //   (*processAndTimeMult)[SensitiveDetectorName][timeSliceName][timeSlice][processName] += 1;
        
      //   for(size_t lp = (*secondary).size()-totsec; lp < (*secondary).size(); ++lp) // loop over secondaries
      //   {
      //     G4String pname = (*secondary)[lp] -> GetParticleDefinition() -> GetParticleName();
      //     if ((*secondary)[lp] -> GetParticleDefinition() -> GetParticleType() == "nucleus"
      //         && (*secondary)[lp] -> GetParticleDefinition() -> GetParticleSubType() == "generic") {
      //       pname = "fragment";
      //     }
          
      //     // check if particle already occurred:
      //     if( (*processAndTimeMap)[SensitiveDetectorName][timeSliceName][timeSlice][processName].find(pname) == (*processAndTimeMap)[SensitiveDetectorName][timeSliceName][timeSlice][processName].end() ){
      //       (*processAndTimeMap)[SensitiveDetectorName][timeSliceName][timeSlice][processName][pname].NParticles = 1;
      //       (*processAndTimeMap)[SensitiveDetectorName][timeSliceName][timeSlice][processName][pname].kinE = (*secondary)[lp] -> GetKineticEnergy();
      //       (*processAndTimeMap)[SensitiveDetectorName][timeSliceName][timeSlice][processName][pname].totE = (*secondary)[lp] -> GetTotalEnergy();
      //     } else {
      //       (*processAndTimeMap)[SensitiveDetectorName][timeSliceName][timeSlice][processName][pname].NParticles++;
      //       (*processAndTimeMap)[SensitiveDetectorName][timeSliceName][timeSlice][processName][pname].kinE += (*secondary)[lp] -> GetKineticEnergy();
      //       (*processAndTimeMap)[SensitiveDetectorName][timeSliceName][timeSlice][processName][pname].totE += (*secondary)[lp] -> GetTotalEnergy();
      //     }
      //   } // end loop over secondaries
      // }
      
      
      //------------------------------------
      // cell-based and particle information
      
      std::map<G4String,std::map<G4String,std::map<G4ThreeVector,std::vector<G4VHit*> > > >* HCMap = CaTSEvt->GetHCMap();
      std::vector<G4VHit*> previousHits = (*HCMap)[SensitiveDetectorName][timeSliceName][cellPosition];
      G4bool hitExists = false;
      for(unsigned int j = 0; j < previousHits.size(); ++j)
      {
        DRTSCalorimeterHit2* aPreviousHit = (DRTSCalorimeterHit2*)(previousHits.at(j));
        if( particleName == aPreviousHit->GetParticleName() &&
            processName  == aPreviousHit->GetProcessName() &&
            timeSlice == aPreviousHit->GetTimeSlice() )
        {
          aPreviousHit->SetEdep(edep + aPreviousHit->GetEdep());
          aPreviousHit->SetEobsbirks(eobs + aPreviousHit->GetEobsbirks());
          aPreviousHit->SetNCeren(NCerenPhotons + aPreviousHit->GetNCeren());
          
          hitExists = true;
          break;
        }
      }
      
      if( hitExists == false )
      {
        DRTSCalorimeterHit2* newHit = new DRTSCalorimeterHit2(particleName,processName,edep,eobs,NCerenPhotons,cellPosition,timeSlice);
        (*HCMap)[SensitiveDetectorName][timeSliceName][cellPosition].push_back(newHit);
      }
    }
  }
  
  
  return true;
}



G4String GetParticleName(G4Track* aTrack)
{
  G4String particleName = aTrack->GetDefinition()->GetParticleName();
  G4String partName = "";
  
  //if( aTrack->GetParticleDefinition()->GetParticleType() == "nucleus" && 
  //    aTrack->GetParticleDefinition()->GetParticleSubType() == "generic" )
  //  partName = "fragment";
  
  if( aTrack->GetParticleDefinition()->GetParticleType() == "nucleus" ) partName = "fragment";
  if( aTrack->GetParticleDefinition()->GetParticleType() == "meson" )   partName = "meson";
  if( aTrack->GetParticleDefinition()->GetParticleType() == "baryon" )  partName = "baryon";
  if( particleName == "proton"  || particleName == "anti_proton" )      partName = "proton";
  if( particleName == "neutron" || particleName == "anti_neutron" )     partName = "neutron";
  if( particleName == "mu-"     || particleName == "mu+" )              partName = "mu";
  if( particleName == "pi+"     || particleName == "pi-" )              partName = "pi";
  if( particleName == "e-"      || particleName == "e+" )               partName = particleName;
  if( particleName == "gamma" )                                         partName = particleName;
  
  if( partName == "proton") {
    if (aTrack->GetKineticEnergy()/CLHEP::MeV < 10.) {
      partName = "proton_ev";
    } else if (aTrack->GetKineticEnergy()/CLHEP::MeV < 500.) {
      partName = "proton_sp";
    } else {
      partName = "proton_he";
    }
  }
  
  if( partName == "" ) partName = "other";
  
  if( partName == "other" ) G4cout << "!!! " << particleName << " saved as 'other'" << G4endl;
  
  return partName;
}
