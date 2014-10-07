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
 Created on May 2, 2012, 8:44 AM
-------------------------------------------------------------------------*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#ifndef Event_HH
#define	Event_HH
#include <vector>
#include <map>

#include "TH1F.h"

#include "products.hh"
#include "G4Types.hh"
#include "G4VHit.hh"

class Event
{
public:
  products evtprod;
  
private:
  G4int fEvtNum;
  G4double TotEnergy; // total energy deposited/Evt
  G4double TotNCeren; // total number of Cerenkov photons/Evt
  G4double TotObsEnergy; // total observed/Evt
  
  // Hit Maps:
  std::map<G4String, std::vector<G4VHit*> > hcmap; // map of Hit Collections
  
  // some info what's happening in an event:
  std::map<G4String, std::map<G4String, G4double> > E_byParticle; // Energy deposited by particle type in a given detector: map<Detector,<particle,energy> >
  std::map<G4String, std::map<G4String, G4double> > Eobs_byParticle; // Birks suppressed energy by particle type in a given detector: map<Detector,<particle,energy> >
  std::map<G4String, std::map<G4String, G4double> >NCeren_byParticle; // Cerenkov photons by particle type in a given detector: map<Detector,<particle,NCeren> >
  std::map<G4String, std::map<G4String, std::map<G4String, products> > > processMap; // keep track of particles produced by a process in a given detector: // map<Detector,map <process,map < particle products> > >
  std::map<G4String, std::map<G4String, G4int > > processMult; // Keeps track how often a specific process occurs  in a given map <detector,map<processes,mult> >
  
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndGlobalTimeLow; // Energy deposited by particle type in a given detector: map<Detector,map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndGlobalTimeMed; // Energy deposited by particle type in a given detector: map<Detector,map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndGlobalTimeHig; // Energy deposited by particle type in a given detector: map<Detector,map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndGlobalTimeLow; // Birks suppressed energy by particle type in a given detector: map<Detector, map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndGlobalTimeMed; // Birks suppressed energy by particle type in a given detector: map<Detector, map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndGlobalTimeHig; // Birks suppressed energy by particle type in a given detector: map<Detector, map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndGlobalTimeLow; // Cerenkov photons by particle type in a given detector: map<Detector, map<timeslice, <particle,NCeren> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndGlobalTimeMed; // Cerenkov photons by particle type in a given detector: map<Detector, map<timeslice, <particle,NCeren> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndGlobalTimeHig; // Cerenkov photons by particle type in a given detector: map<Detector, map<timeslice, <particle,NCeren> > >
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndLocalTime1Low; // Energy deposited by particle type in a given detector: map<Detector,map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndLocalTime1Med; // Energy deposited by particle type in a given detector: map<Detector,map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndLocalTime1Hig; // Energy deposited by particle type in a given detector: map<Detector,map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndLocalTime1Low; // Birks suppressed energy by particle type in a given detector: map<Detector, map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndLocalTime1Med; // Birks suppressed energy by particle type in a given detector: map<Detector, map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndLocalTime1Hig; // Birks suppressed energy by particle type in a given detector: map<Detector, map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndLocalTime1Low; // Cerenkov photons by particle type in a given detector: map<Detector, map<timeslice, <particle,NCeren> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndLocalTime1Med; // Cerenkov photons by particle type in a given detector: map<Detector, map<timeslice, <particle,NCeren> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndLocalTime1Hig; // Cerenkov photons by particle type in a given detector: map<Detector, map<timeslice, <particle,NCeren> > >
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndLocalTime2Low; // Energy deposited by particle type in a given detector: map<Detector,map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndLocalTime2Med; // Energy deposited by particle type in a given detector: map<Detector,map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > E_byParticleAndLocalTime2Hig; // Energy deposited by particle type in a given detector: map<Detector,map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndLocalTime2Low; // Birks suppressed energy by particle type in a given detector: map<Detector, map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndLocalTime2Med; // Birks suppressed energy by particle type in a given detector: map<Detector, map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > Eobs_byParticleAndLocalTime2Hig; // Birks suppressed energy by particle type in a given detector: map<Detector, map<timeslice, <particle,energy> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndLocalTime2Low; // Cerenkov photons by particle type in a given detector: map<Detector, map<timeslice, <particle,NCeren> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndLocalTime2Med; // Cerenkov photons by particle type in a given detector: map<Detector, map<timeslice, <particle,NCeren> > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > > NCeren_byParticleAndLocalTime2Hig; // Cerenkov photons by particle type in a given detector: map<Detector, map<timeslice, <particle,NCeren> > >
  
  
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > processAndGlobalTimeLowMap; // keep track of particles produced by a process in a given detector: map<Detector, map< timeslice, map <process,map < particle products> > > >
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > processAndGlobalTimeMedMap; // keep track of particles produced by a process in a given detector: map<Detector, map< timeslice, map <process,map < particle products> > > >
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > processAndGlobalTimeHigMap; // keep track of particles produced by a process in a given detector: map<Detector, map< timeslice, map <process,map < particle products> > > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > > processAndGlobalTimeLowMult; // Keeps track how often a specific process occurs in a given map <detector,map<processes,mult> >
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > > processAndGlobalTimeMedMult; // Keeps track how often a specific process occurs in a given map <detector,map<processes,mult> >
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > > processAndGlobalTimeHigMult; // Keeps track how often a specific process occurs in a given map <detector,map<processes,mult> >
  
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > processAndLocalTime1LowMap; // keep track of particles produced by a process in a given detector: map<Detector, map< timeslice, map <process,map < particle products> > > >
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > processAndLocalTime1MedMap; // keep track of particles produced by a process in a given detector: map<Detector, map< timeslice, map <process,map < particle products> > > >
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > processAndLocalTime1HigMap; // keep track of particles produced by a process in a given detector: map<Detector, map< timeslice, map <process,map < particle products> > > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > > processAndLocalTime1LowMult; // Keeps track how often a specific process occurs in a given map <detector,map<processes,mult> >
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > > processAndLocalTime1MedMult; // Keeps track how often a specific process occurs in a given map <detector,map<processes,mult> >
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > > processAndLocalTime1HigMult; // Keeps track how often a specific process occurs in a given map <detector,map<processes,mult> >
  
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > processAndLocalTime2LowMap; // keep track of particles produced by a process in a given detector: map<Detector, map< timeslice, map <process,map < particle products> > > >
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > processAndLocalTime2MedMap; // keep track of particles produced by a process in a given detector: map<Detector, map< timeslice, map <process,map < particle products> > > >
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > > processAndLocalTime2HigMap; // keep track of particles produced by a process in a given detector: map<Detector, map< timeslice, map <process,map < particle products> > > >
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > > processAndLocalTime2LowMult; // Keeps track how often a specific process occurs in a given map <detector,map<processes,mult> >
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > > processAndLocalTime2MedMult; // Keeps track how often a specific process occurs in a given map <detector,map<processes,mult> >
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > > processAndLocalTime2HigMult; // Keeps track how often a specific process occurs in a given map <detector,map<processes,mult> >
  
  
public:
  
  Event() : fEvtNum(0), TotEnergy(0), TotNCeren(0), TotObsEnergy(0) {}
  
  virtual ~Event() {}
  
  
  void SetEventNr(G4int i) {
    fEvtNum = i;
  }
  void SetTotEnergy(G4double e) {
    TotEnergy = e;
  }
  void SetTotNCeren(G4double e) {
    TotNCeren = e;
  }
  void SetTotObsEnergy(G4double e) {
    TotObsEnergy = e;
  }
  G4int GetEventNumber() const {
    return fEvtNum;
  }
  G4double GetTotEnergy() const {
    return TotEnergy;
  }
  G4double GetNCeren() const {
    return TotNCeren;
  }
  G4double GetTotObsEnergy() const {
    return TotObsEnergy;
  }
  
  
  std::map<G4String, std::vector<G4VHit* > >* GetHCMap() {
    return &hcmap;
  }
  std::map<G4String, std::map<G4String, G4double> >* GetE_byParticle() {
    return &E_byParticle;
  }
  std::map<G4String, std::map<G4String, G4double> >* GetEobs_byParticle() {
    return &Eobs_byParticle;
  }
  std::map<G4String, std::map<G4String, G4double > >* GetNCeren_byParticle() {
    return &NCeren_byParticle;
  }
  
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetE_byParticleAndGlobalTimeLow() {
    return &E_byParticleAndGlobalTimeLow;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetE_byParticleAndGlobalTimeMed() {
    return &E_byParticleAndGlobalTimeMed;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetE_byParticleAndGlobalTimeHig() {
    return &E_byParticleAndGlobalTimeHig;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetEobs_byParticleAndGlobalTimeLow() {
    return &Eobs_byParticleAndGlobalTimeLow;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetEobs_byParticleAndGlobalTimeMed() {
    return &Eobs_byParticleAndGlobalTimeMed;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetEobs_byParticleAndGlobalTimeHig() {
    return &Eobs_byParticleAndGlobalTimeHig;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* GetNCeren_byParticleAndGlobalTimeLow() {
    return &NCeren_byParticleAndGlobalTimeLow;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* GetNCeren_byParticleAndGlobalTimeMed() {
    return &NCeren_byParticleAndGlobalTimeMed;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* GetNCeren_byParticleAndGlobalTimeHig() {
    return &NCeren_byParticleAndGlobalTimeHig;
  }
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetE_byParticleAndLocalTime1Low() {
    return &E_byParticleAndLocalTime1Low;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetE_byParticleAndLocalTime1Med() {
    return &E_byParticleAndLocalTime1Med;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetE_byParticleAndLocalTime1Hig() {
    return &E_byParticleAndLocalTime1Hig;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetEobs_byParticleAndLocalTime1Low() {
    return &Eobs_byParticleAndLocalTime1Low;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetEobs_byParticleAndLocalTime1Med() {
    return &Eobs_byParticleAndLocalTime1Med;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetEobs_byParticleAndLocalTime1Hig() {
    return &Eobs_byParticleAndLocalTime1Hig;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* GetNCeren_byParticleAndLocalTime1Low() {
    return &NCeren_byParticleAndLocalTime1Low;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* GetNCeren_byParticleAndLocalTime1Med() {
    return &NCeren_byParticleAndLocalTime1Med;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* GetNCeren_byParticleAndLocalTime1Hig() {
    return &NCeren_byParticleAndLocalTime1Hig;
  }
  
  
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetE_byParticleAndLocalTime2Low() {
    return &E_byParticleAndLocalTime2Low;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetE_byParticleAndLocalTime2Med() {
    return &E_byParticleAndLocalTime2Med;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetE_byParticleAndLocalTime2Hig() {
    return &E_byParticleAndLocalTime2Hig;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetEobs_byParticleAndLocalTime2Low() {
    return &Eobs_byParticleAndLocalTime2Low;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetEobs_byParticleAndLocalTime2Med() {
    return &Eobs_byParticleAndLocalTime2Med;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double> > >* GetEobs_byParticleAndLocalTime2Hig() {
    return &Eobs_byParticleAndLocalTime2Hig;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* GetNCeren_byParticleAndLocalTime2Low() {
    return &NCeren_byParticleAndLocalTime2Low;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* GetNCeren_byParticleAndLocalTime2Med() {
    return &NCeren_byParticleAndLocalTime2Med;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4double > > >* GetNCeren_byParticleAndLocalTime2Hig() {
    return &NCeren_byParticleAndLocalTime2Hig;
  }
  
  
  std::map<G4String, std::map<G4String, std::map<G4String, products> > >* GetProcessMap() {
    return &processMap;
  }
  std::map<G4String, std::map<G4String, G4int > >* GetProcessMult() {
    return &processMult;
  }
  
  
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* GetProcessAndGlobalTimeLowMap() {
    return &processAndGlobalTimeLowMap;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* GetProcessAndGlobalTimeMedMap() {
    return &processAndGlobalTimeMedMap;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* GetProcessAndGlobalTimeHigMap() {
    return &processAndGlobalTimeHigMap;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > >* GetProcessAndGlobalTimeLowMult() {
    return &processAndGlobalTimeLowMult;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > >* GetProcessAndGlobalTimeMedMult() {
    return &processAndGlobalTimeMedMult;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > >* GetProcessAndGlobalTimeHigMult() {
    return &processAndGlobalTimeHigMult;
  }
  
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* GetProcessAndLocalTime1LowMap() {
    return &processAndLocalTime1LowMap;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* GetProcessAndLocalTime1MedMap() {
    return &processAndLocalTime1MedMap;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* GetProcessAndLocalTime1HigMap() {
    return &processAndLocalTime1HigMap;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > >* GetProcessAndLocalTime1LowMult() {
    return &processAndLocalTime1LowMult;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > >* GetProcessAndLocalTime1MedMult() {
    return &processAndLocalTime1MedMult;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > >* GetProcessAndLocalTime1HigMult() {
    return &processAndLocalTime1HigMult;
  }
  
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* GetProcessAndLocalTime2LowMap() {
    return &processAndLocalTime2LowMap;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* GetProcessAndLocalTime2MedMap() {
    return &processAndLocalTime2MedMap;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, std::map<G4String, products> > > >* GetProcessAndLocalTime2HigMap() {
    return &processAndLocalTime2HigMap;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > >* GetProcessAndLocalTime2LowMult() {
    return &processAndLocalTime2LowMult;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > >* GetProcessAndLocalTime2MedMult() {
    return &processAndLocalTime2MedMult;
  }
  std::map<G4String, std::map<G4int, std::map<G4String, G4int > > >* GetProcessAndLocalTime2HigMult() {
    return &processAndLocalTime2HigMult;
  }
  
  
  void Reset();
};
#endif	/* Event_HH */
