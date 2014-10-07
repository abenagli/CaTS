
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

void Event::Reset()
{
  TotEnergy = 0.0;
  TotObsEnergy = 0.0;
  TotNCeren = 0.0;
  
  hcmap.clear();
  
  E_byParticle.clear();
  Eobs_byParticle.clear();
  NCeren_byParticle.clear();
  E_byParticleAndGlobalTimeLow.clear();
  E_byParticleAndGlobalTimeMed.clear();
  E_byParticleAndGlobalTimeHig.clear();
  E_byParticleAndLocalTime1Low.clear();
  E_byParticleAndLocalTime1Med.clear();
  E_byParticleAndLocalTime1Hig.clear();
  E_byParticleAndLocalTime2Low.clear();
  E_byParticleAndLocalTime2Med.clear();
  E_byParticleAndLocalTime2Hig.clear();
  Eobs_byParticleAndGlobalTimeLow.clear();
  Eobs_byParticleAndGlobalTimeMed.clear();
  Eobs_byParticleAndGlobalTimeHig.clear();
  Eobs_byParticleAndLocalTime1Low.clear();
  Eobs_byParticleAndLocalTime1Med.clear();
  Eobs_byParticleAndLocalTime1Hig.clear();
  Eobs_byParticleAndLocalTime2Low.clear();
  Eobs_byParticleAndLocalTime2Med.clear();
  Eobs_byParticleAndLocalTime2Hig.clear();
  NCeren_byParticleAndGlobalTimeLow.clear();
  NCeren_byParticleAndGlobalTimeMed.clear();
  NCeren_byParticleAndGlobalTimeHig.clear();
  NCeren_byParticleAndLocalTime1Low.clear();
  NCeren_byParticleAndLocalTime1Med.clear();
  NCeren_byParticleAndLocalTime1Hig.clear();
  NCeren_byParticleAndLocalTime2Low.clear();
  NCeren_byParticleAndLocalTime2Med.clear();
  NCeren_byParticleAndLocalTime2Hig.clear();
  
  processMap.clear();
  processMult.clear();
  processAndGlobalTimeLowMap.clear();
  processAndGlobalTimeMedMap.clear();
  processAndGlobalTimeHigMap.clear();
  processAndGlobalTimeLowMult.clear();
  processAndGlobalTimeMedMult.clear();
  processAndGlobalTimeHigMult.clear();
  
  // std::map<G4String, std::map<G4String, G4double> >::iterator Eiter;
  // for (Eiter = E_byParticle.begin(); Eiter != E_byParticle.end(); Eiter++)
  // {
  //   G4cout << (*Eiter).first << G4endl;
  //   std::map<G4String, G4double> partmap = (*Eiter).second;
  //   std::map<G4String, G4double>::iterator partiter;
  //   for (partiter = partmap.begin(); partiter != partmap.end(); partiter++)
  //   {
  //     G4cout << (*partiter).first << "   " << (*partiter).second << G4endl;
  //     (*partiter).second = 0.0;
  //   }
  // }
}
