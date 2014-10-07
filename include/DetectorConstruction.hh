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

#ifndef _DETECTORCONSTRUCTION_H_
#define _DETECTORCONSTRUCTION_H_

#include "G4VUserDetectorConstruction.hh"
#include "G4GDMLParser.hh"



class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  DetectorConstruction(G4String fname);
  
  G4VPhysicalVolume *Construct()
  {
    return ConstructDetector();
  }
  
  std::vector<G4String>* GetVolumes() { return Volumes; };
  
  static DetectorConstruction* getInstance() { return instance; };
  
private:
  
  G4VPhysicalVolume *World;
  G4VPhysicalVolume * ConstructDetector();
  G4String gdmlFile;
  G4GDMLParser *parser;
  
  // Extended reader
  //
  G4GDMLReadStructure* fReader;
  
  std::vector<G4String>* Volumes;
  
  static DetectorConstruction* instance;
};

#endif
