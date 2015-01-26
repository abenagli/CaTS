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
  
  std::vector<G4String>* GetSolids()             { return Solids;            };
  std::vector<G4float>* GetSolidsXHalfLength()  { return SolidsXHalfLength; };
  std::vector<G4float>* GetSolidsYHalfLength()  { return SolidsYHalfLength; };
  std::vector<G4float>* GetSolidsZHalfLength()  { return SolidsZHalfLength; };
  std::vector<G4String>* GetVolumes()            { return Volumes;           };

  G4int GetXCellNum()  { return numXCells;  };
  G4int GetYCellNum()  { return numYCells;  };
  G4int GetZLayerNum() { return numZLayers; };
  
  static DetectorConstruction* getInstance() { return instance; };
  
private:
  
  G4VPhysicalVolume *World;
  G4VPhysicalVolume * ConstructDetector();
  G4String gdmlFile;
  G4GDMLParser *parser;
  
  // Extended reader
  //
  G4GDMLReadStructure* fReader;
  
  std::vector<G4String>* Solids;
  std::vector<G4float>* SolidsXHalfLength;
  std::vector<G4float>* SolidsYHalfLength;
  std::vector<G4float>* SolidsZHalfLength;
  std::vector<G4String>* Volumes;
  
  G4int numXCells;
  G4int numYCells;
  G4int numZLayers;
  
  static DetectorConstruction* instance;
};

#endif
