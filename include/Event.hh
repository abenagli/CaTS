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

#ifndef Event_HH
#define	Event_HH

#include <vector>
#include <map>

#include "TH1F.h"

#include "products.hh"
#include "G4Types.hh"
#include "G4VHit.hh"
#include "G4ThreeVector.hh"



class Event
{
public:
  products evtprod;
  
private:
  
  // global variables
  G4int fEvtNum;
  G4float TotDepAndEscAndLosEnergy; // total energy deposited+escaped+lost/Evt
  G4float TotDepAndEscEnergy;       // total energy deposited+escaped/Evt
  G4float TotDepEnergy;             // total energy deposited/Evt
  G4float TotObsEnergy;             // total energy observed/Evt
  G4float TotNCeren;                // total number of Cerenkov photons/Evt
  
  
  std::map<G4String,G4float> m_EdepAndEscAndLos; // map<Detector,energy>
  std::map<G4String,G4float> m_EdepAndEsc;       // map<Detector,energy>
  std::map<G4String,G4float> m_Edep;             // map<Detector,energy>
  std::map<G4String,G4float> m_Eobs;             // map<Detector,energy>
  std::map<G4String,G4float> m_NCeren;           // map<Detector,energy>
  
  std::map<G4String,std::map<G4String,std::map<G4int,G4float> > > m_Edep_byTime;       // map<Detector,map<timeSliceType,map<timeSlice,energy> > >
  std::map<G4String,std::map<G4String,std::map<G4int,G4float> > > m_Eobs_byTime;       // map<Detector,map<timeSliceType,map<timeSlice,energy> > >
  std::map<G4String,std::map<G4String,std::map<G4int,G4float> > > m_NCeren_byTime;     // map<Detector,map<timeSliceType,map<timeSlice,energy> > >
  
  
  std::map<G4String,std::map<G4ThreeVector,G4float> > m_Edep_byPos;   // map<Detector,map<position,energy> >
  std::map<G4String,std::map<G4ThreeVector,G4float> > m_Eobs_byPos;   // map<Detector,map<position,energy> >
  std::map<G4String,std::map<G4ThreeVector,G4float> > m_NCeren_byPos; // map<Detector,map<position,energy> >
  
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4ThreeVector,G4float> > > > m_Edep_byPosAndTime;   // map<Detector,map<timeSliceType,map<timeSlice,map<position,energy> > > >
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4ThreeVector,G4float> > > > m_Eobs_byPosAndTime;   // map<Detector,map<timeSliceType,map<timeSlice,map<position,energy> > > >
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4ThreeVector,G4float> > > > m_NCeren_byPosAndTime; // map<Detector,map<timeSliceType,map<timeSlice,map<position,energy> > > >
  
  
  std::map<G4String,std::map<G4String,G4float> > m_Edep_byParticle;   // map<Detector,map<particle,energy> >
  std::map<G4String,std::map<G4String,G4float> > m_Eobs_byParticle;   // map<Detector,map<particle,energy> >
  std::map<G4String,std::map<G4String,G4float> > m_NCeren_byParticle; // map<Detector,map<particle,NCeren> >
  std::map<G4String,std::map<G4String,G4float> > m_Edep_byParticleType;   // map<Detector,map<particleType,energy> >
  std::map<G4String,std::map<G4String,G4float> > m_Eobs_byParticleType;   // map<Detector,map<particleType,energy> >
  std::map<G4String,std::map<G4String,G4float> > m_NCeren_byParticleType; // map<Detector,map<particleType,NCeren> >
  
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > > m_Edep_byParticleAndTime;   // map<Detector,map<timeSliceType,map<timeSlice,map<particle,energy> > > >
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > > m_Eobs_byParticleAndTime;   // map<Detector,map<timeSliceType,map<timeSlice,map<particle,energy> > > >
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > > m_NCeren_byParticleAndTime; // map<Detector,map<timeSliceType,map<timeSlice,map<particle,energy> > > >
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > > m_Edep_byParticleTypeAndTime;   // map<Detector,map<timeSliceType,map<timeSlice,map<particleType,energy> > > >
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > > m_Eobs_byParticleTypeAndTime;   // map<Detector,map<timeSliceType,map<timeSlice,map<particleType,energy> > > >
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > > m_NCeren_byParticleTypeAndTime; // map<Detector,map<timeSliceType,map<timeSlice,map<particleType,energy> > > >
  
  
  std::map<G4String,std::map<G4String,G4int> > m_particleMult;   // map<Detector,map<particle,mult> >
  std::map<G4String,std::map<G4String,TH1F*> > m_particleHist;   // map<Detector,map<particle,TH1F*> >
  
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4int> > > > m_particleAndTimeMult;   // map<Detector,map<timeSliceType,map<timeSlice,map<particle,mult> > > >
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,TH1F*> > > > m_particleAndTimeHist;   // map<Detector,map<timeSliceType,map<timeSlice,map<particle,mult> > > >
  
  
  std::map<G4String,std::map<G4String,G4int> > processMult;                       // keep track how often a specific process occurs: map<Detector,map<processes,mult> >
  std::map<G4String,std::map<G4String,std::map<G4String,products> > > processMap; // keep track of particles produced by a processr: map<Detector,map<process,map< particle,products> > >
  
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4int> > > > processAndTimeMult;                       // map<Detector,map<timeSliceType,map<timeSlice,map<processes,mult> > > >
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,std::map<G4String,products> > > > > processAndTimeMap; // map<Detector,map<timeSliceType,map<timeSlice,map<process,map< particle,products> > > > >
  
  
  std::map<G4String,std::map<G4String,std::map<G4ThreeVector,std::vector<G4VHit*> > > > HCMap; // map<Detector,map<timeSliceType,map<position,hit> > >
  
  
public:
  
  //! ctor
  Event();
  
  //! dtor
  virtual ~Event() {};
  
  
  //! methods
  void SetEventNr(G4int i) { fEvtNum = i; };
  void SetTotDepAndEscAndLosEnergy(G4float e) { TotDepAndEscAndLosEnergy = e; };
  void SetTotDepAndEscEnergy(G4float e)       { TotDepAndEscEnergy = e;       };
  void SetTotDepEnergy(G4float e)             { TotDepEnergy = e;             };
  void SetTotObsEnergy(G4float e)             { TotObsEnergy = e;             };
  void SetTotNCeren(G4float n)                { TotNCeren = n;                };
  
  G4int GetEventNumber() const { return fEvtNum; }
  G4float GetTotDepAndEscAndLosEnergy() const { return TotDepAndEscAndLosEnergy; };
  G4float GetTotDepAndEscEnergy()       const { return TotDepAndEscEnergy;       };
  G4float GetTotDepEnergy()             const { return TotDepEnergy;             };
  G4float GetTotObsEnergy()             const { return TotObsEnergy;             };
  G4float GetNCeren()                   const { return TotNCeren;                };
  
  
  std::map<G4String,G4float>* GetEdepAndEscAndLosMap() { return &m_EdepAndEscAndLos; };
  std::map<G4String,G4float>* GetEdepAndEscMap()       { return &m_EdepAndEsc;       };
  std::map<G4String,G4float>* GetEdepMap()             { return &m_Edep;             };
  std::map<G4String,G4float>* GetEobsMap()             { return &m_Eobs;             };
  std::map<G4String,G4float>* GetNCerenMap()           { return &m_NCeren;           };
    
  std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* GetEdepByTimeMap()   { return &m_Edep_byTime;   };
  std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* GetEobsByTimeMap()   { return &m_Eobs_byTime;   };
  std::map<G4String,std::map<G4String,std::map<G4int,G4float> > >* GetNCerenByTimeMap() { return &m_NCeren_byTime; };
  
  
  std::map<G4String,std::map<G4ThreeVector,G4float> >* GetEdepByPosMap()   { return &m_Edep_byPos;   };
  std::map<G4String,std::map<G4ThreeVector,G4float> >* GetEobsByPosMap()   { return &m_Eobs_byPos;   };
  std::map<G4String,std::map<G4ThreeVector,G4float> >* GetNCerenByPosMap() { return &m_NCeren_byPos; };
  
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4ThreeVector,G4float> > > >* GetEdepByPosAndTimeMap()   { return &m_Edep_byPosAndTime;   };
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4ThreeVector,G4float> > > >* GetEobsByPosAndTimeMap()   { return &m_Eobs_byPosAndTime;   };
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4ThreeVector,G4float> > > >* GetNCerenByPosAndTimeMap() { return &m_NCeren_byPosAndTime; };
  
  
  std::map<G4String,std::map<G4String,G4float> >* GetEdepByParticleMap()   { return &m_Edep_byParticle;   };
  std::map<G4String,std::map<G4String,G4float> >* GetEobsByParticleMap()   { return &m_Eobs_byParticle;   };
  std::map<G4String,std::map<G4String,G4float> >* GetNCerenByParticleMap() { return &m_NCeren_byParticle; };
  std::map<G4String,std::map<G4String,G4float> >* GetEdepByParticleTypeMap()   { return &m_Edep_byParticleType;   };
  std::map<G4String,std::map<G4String,G4float> >* GetEobsByParticleTypeMap()   { return &m_Eobs_byParticleType;   };
  std::map<G4String,std::map<G4String,G4float> >* GetNCerenByParticleTypeMap() { return &m_NCeren_byParticleType; };
  
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > >* GetEdepByParticleAndTimeMap()   { return &m_Edep_byParticleAndTime;   };
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > >* GetEobsByParticleAndTimeMap()   { return &m_Eobs_byParticleAndTime;   };
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > >* GetNCerenByParticleAndTimeMap() { return &m_NCeren_byParticleAndTime; };
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > >* GetEdepByParticleTypeAndTimeMap()   { return &m_Edep_byParticleTypeAndTime;   };
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > >* GetEobsByParticleTypeAndTimeMap()   { return &m_Eobs_byParticleTypeAndTime;   };
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4float> > > >* GetNCerenByParticleTypeAndTimeMap() { return &m_NCeren_byParticleTypeAndTime; };
  
  
  std::map<G4String,std::map<G4String,G4int> >* GetParticleMult() { return &m_particleMult; };
  std::map<G4String,std::map<G4String,TH1F*> >* GetParticleHist() { return &m_particleHist; };
  
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4int> > > >* GetParticleAndTimeMult() { return &m_particleAndTimeMult; };
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,TH1F*> > > >* GetParticleAndTimeHist() { return &m_particleAndTimeHist; };
  
  
  std::map<G4String,std::map<G4String,G4int> >* GetProcessMult() { return &processMult; };
  std::map<G4String,std::map<G4String,std::map<G4String, products> > >* GetProcessMap() { return &processMap; };
  
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,G4int> > > >* GetProcessAndTimeMult() { return &processAndTimeMult; };
  std::map<G4String,std::map<G4String,std::map<G4int,std::map<G4String,std::map<G4String, products> > > > >* GetProcessAndTimeMap() { return &processAndTimeMap; };
  
  
  std::map<G4String,std::map<G4String,std::map<G4ThreeVector,std::vector<G4VHit*> > > >* GetHCMap() { return &HCMap; };
  
  
  void Reset();
};
#endif
