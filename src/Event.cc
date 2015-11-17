
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
#include "Event.hh"



Event::Event():
  fEvtNum(0),
  TotDepAndEscAndLosEnergy(0.),
  TotDepAndEscEnergy(0.),
  TotDepEnergy(0.),
  TotObsEnergy(0.),
  TotNCeren(0.)
{}



void Event::Reset()
{
  TotDepAndEscAndLosEnergy = 0.0;
  TotDepAndEscEnergy = 0.0;
  TotDepEnergy = 0.0;
  TotObsEnergy = 0.0;
  TotNCeren = 0.0;
  
  
  m_EdepAndEscAndLos.clear();
  m_EdepAndEsc.clear();
  m_Edep.clear();
  m_Eobs.clear();
  m_NCeren.clear();
  
  m_Edep_byTime.clear();
  m_Eobs_byTime.clear();
  m_NCeren_byTime.clear();
  
  
  m_Edep_byPos.clear();
  m_Eobs_byPos.clear();
  m_NCeren_byPos.clear();
  
  m_Edep_byPosAndTime.clear();
  m_Eobs_byPosAndTime.clear();
  m_NCeren_byPosAndTime.clear();
  
  
  m_Edep_byParticle.clear();
  m_Eobs_byParticle.clear();
  m_NCeren_byParticle.clear();
  m_Edep_byParticleType.clear();
  m_Eobs_byParticleType.clear();
  m_NCeren_byParticleType.clear();
  
  m_Edep_byParticleAndTime.clear();
  m_Eobs_byParticleAndTime.clear();
  m_NCeren_byParticleAndTime.clear();
  m_Edep_byParticleTypeAndTime.clear();
  m_Eobs_byParticleTypeAndTime.clear();
  m_NCeren_byParticleTypeAndTime.clear();
  
  
  m_particleMult.clear();
  
  m_particleAndTimeMult.clear();
  
  
  processMult.clear();
  processMap.clear();
  
  processAndTimeMult.clear();
  processAndTimeMap.clear();
  
  
  HCMap.clear();
}
