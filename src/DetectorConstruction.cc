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
///////////////////////////////////////////////////////////////////////////
// File: DetectorConstruction.cc
// Description: Performs the detector construction
///////////////////////////////////////////////////////////////////////////

#include "DetectorConstruction.hh"
#include "G4SDManager.hh"
#include "TrackerSD.hh"
#include "PhotonSD.hh"
#include "CalorimeterSD.hh"
#include "DRCalorimeterSD.hh"
#include "DRTSCalorimeterSD.hh"
#include "StoppingCalorimeterSD.hh"
#include "ColorReader.hh"
//#include "G4Colour.hh"

DetectorConstruction* DetectorConstruction::instance = 0;



DetectorConstruction::DetectorConstruction(G4String fname):
  gdmlFile(fname),
  Solids(new std::vector<G4String>),
  SolidsXHalfLength(new std::vector<G4float>),
  SolidsYHalfLength(new std::vector<G4float>),
  SolidsZHalfLength(new std::vector<G4float>),
  Volumes(new std::vector<G4String>)
{
  instance = this;
}



G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{
  fReader = new ColorReader;
  //  fWriter = new G03ColorWriter;
  //fParser = new G4GDMLParser(fReader, fWriter);
  parser = new G4GDMLParser(fReader);
  //    parser.Read(gdmlFile);
  //   World = parser.GetWorldVolume();
  parser->Read(gdmlFile,true);
  World = parser->GetWorldVolume();
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  const G4GDMLAuxMapType* auxmap = parser->GetAuxMap();
  G4cout << "\nFound " << auxmap->size()
         << " volume(s) with auxiliary information."
         << G4endl;
  for (G4GDMLAuxMapType::const_iterator iter = auxmap->begin(); iter != auxmap->end(); iter++)
  {
    G4cout << "\nVolume " << ((*iter).first)->GetName()
           << " has the following list of auxiliary information: "
           << G4endl;
    
    for (G4GDMLAuxListType::const_iterator vit = (*iter).second.begin(); vit != (*iter).second.end(); vit++)
    {
      G4cout << "--> Type: " << (*vit).type
             << " Value: " << (*vit).value << G4endl;
      
      G4String name = ((*iter).first)->GetName();
      
      if ((*vit).type == "SensDet")
      {
        if ((*vit).value == "DRCalorimeter") {
          name += "_DRCalorimeter";
          DRCalorimeterSD* aDRCalorimeterSD = new DRCalorimeterSD(name);
          SDman->AddNewDetector(aDRCalorimeterSD);
          ((*iter).first)->SetSensitiveDetector(aDRCalorimeterSD);
        }
        else if ((*vit).value == "DRTSCalorimeter") {
          name += "_DRTSCalorimeter";
          DRTSCalorimeterSD* aDRTSCalorimeterSD = new DRTSCalorimeterSD(name);
          SDman->AddNewDetector(aDRTSCalorimeterSD);
          ((*iter).first)->SetSensitiveDetector(aDRTSCalorimeterSD);
        }
        else if ((*vit).value == "Calorimeter") {
          name += "_Calorimeter";
          CalorimeterSD* aCalorimeterSD = new CalorimeterSD(name);
          SDman->AddNewDetector(aCalorimeterSD);
          ((*iter).first)->SetSensitiveDetector(aCalorimeterSD);
        }
        else if ((*vit).value == "StoppingCalorimeter") {
          name += "_StoppingCalorimeter";
          StoppingCalorimeterSD* aStoppingCalorimeterSD = new StoppingCalorimeterSD(name);
          SDman->AddNewDetector(aStoppingCalorimeterSD);
          ((*iter).first)->SetSensitiveDetector(aStoppingCalorimeterSD);
        }
        else if ((*vit).value == "Tracker") {
          name += "_Tracker";
          TrackerSD* aTrackerSD = new TrackerSD(name);
          SDman->AddNewDetector(aTrackerSD);
          ((*iter).first)->SetSensitiveDetector(aTrackerSD);
        }
        else if ((*vit).value == "PhotonDetector") {
          name += "_PhotonDetector";
          PhotonSD* aPhotonSD = new PhotonSD(name);
          SDman->AddNewDetector(aPhotonSD);
          ((*iter).first)->SetSensitiveDetector(aPhotonSD);
        }
        
        G4cout << "Attaching sensitive Detector: " << (*vit).value
               << " to Volume:  " << ((*iter).first)->GetName()
               << " name: " << name 
               << G4endl;
        
        Volumes -> push_back(name);
      }
      SDman->ListTree();
      //            else if ((*vit).type == "Color") {
      //                    std::cout << "!!!!!!!!!!!! Attaching Color: " << (*vit).value
      //                            << " to Volume:  " << ((*iter).first) << std::endl;
      //            }
    }
    
    
    Solids -> push_back( (((*iter).first)->GetSolid())->GetName() );
    if( (((*iter).first)->GetSolid())->GetEntityType() == "G4Box" )
    {
      SolidsXHalfLength -> push_back( ((G4Box*)(((*iter).first)->GetSolid()))->GetXHalfLength() );
      SolidsYHalfLength -> push_back( ((G4Box*)(((*iter).first)->GetSolid()))->GetYHalfLength() );
      SolidsZHalfLength -> push_back( ((G4Box*)(((*iter).first)->GetSolid()))->GetZHalfLength() );
    }
    else
    {
      SolidsXHalfLength -> push_back( -999. );
      SolidsYHalfLength -> push_back( -999. );
      SolidsZHalfLength -> push_back( -999. );
    }
  }
  
  G4cout << G4endl;
  
  
  G4float var;
  var = parser->GetConstant("numcol");
  numXCells = var+1;
  var = parser->GetConstant("numrow");
  numYCells = var+1;
  var = parser->GetConstant("numlay");
  numZLayers = var+1;  
  
  return World;
}
