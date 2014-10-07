/* 
 * File:   DRTSCalorimeterSDMessenger.hh
 * Author: andread
 *
 * Created on July 10, 2013, 4:38 PM
 */

#ifndef DRTSCALORIMETERSDMESSENGER_HH
#define	DRTSCALORIMETERSDMESSENGER_HH
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

#include "G4UImessenger.hh"
#include "globals.hh"

class DRTSCalorimeterSD;
class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DRTSCalorimeterSDMessenger : public G4UImessenger {

public:
  DRTSCalorimeterSDMessenger(DRTSCalorimeterSD*);
  ~DRTSCalorimeterSDMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  G4String GetCurrentValue(G4UIcommand*);
  
private:
  DRTSCalorimeterSD* pDRTSCalorimeterSD;
  G4UIdirectory* AnalysisDir;
  G4UIcmdWithADouble* psetc1Cmd;
  G4UIcmdWithADouble* psetc2Cmd;
  G4UIcmdWithADoubleAndUnit* psetc3Cmd;
  G4UIcmdWithADoubleAndUnit* psetc3MinCmd;
  G4UIcmdWithADoubleAndUnit* psetc3MaxCmd;
  G4UIcmdWithADoubleAndUnit* psetc4Cmd;
  G4UIcmdWithADoubleAndUnit* psetc4MinCmd;
  G4UIcmdWithADoubleAndUnit* psetc4MaxCmd;
  G4UIcmdWithADoubleAndUnit* psetc5Cmd;
  G4UIcmdWithADoubleAndUnit* psetc5MinCmd;
  G4UIcmdWithADoubleAndUnit* psetc5MaxCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif	/* DRTSCALORIMETERSDMESSENGER_HH */
