
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
  TotDepEnergy(0.),
  TotObsEnergy(0.),
  TotNCeren(0.)
{}



void Event::Reset()
{
  TotDepEnergy = 0.0;
  TotObsEnergy = 0.0;
  TotNCeren = 0.0;
  
  
  m_Edep.clear();
  m_Eobs.clear();
  m_NCeren.clear();
  
  m_Edep_byTime.clear();
  m_Eobs_byTime.clear();
  m_NCeren_byTime.clear();
  
  
  m_Edep_byParticle.clear();
  m_Eobs_byParticle.clear();
  m_NCeren_byParticle.clear();
  m_Edep_byProcess.clear();
  m_Eobs_byProcess.clear();
  m_NCeren_byProcess.clear();
  
  m_Edep_byParticleAndTime.clear();
  m_Eobs_byParticleAndTime.clear();
  m_NCeren_byParticleAndTime.clear();
  m_Edep_byProcessAndTime.clear();
  m_Eobs_byProcessAndTime.clear();
  m_NCeren_byProcessAndTime.clear();
  
  
  processMult.clear();
  processMap.clear();
  
  processAndTimeMult.clear();
  processAndTimeMap.clear();
  
  
  HCMap.clear();
}
