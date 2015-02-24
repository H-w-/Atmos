//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B2aDetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B2aDetectorConstruction.cc
/// \brief Implementation of the B2aDetectorConstruction class

#include <cstring>
#include <sstream>
#include <cmath>
 
#include "B2aDetectorConstruction.hh"
#include "B2aDetectorMessenger.hh"
//#include "B2TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4ThreadLocal 
G4GlobalMagFieldMessenger* B2aDetectorConstruction::fMagFieldMessenger = 0;



B2aDetectorConstruction::B2aDetectorConstruction()
:G4VUserDetectorConstruction(), 
 fNbOfChambers(10),
 fLogicTarget(NULL), fLogicChamber(NULL), 
 fTargetMaterial(NULL), fChamberMaterial(NULL), 
 fStepLimit(NULL),
 fCheckOverlaps(true)
{
  fMessenger = new B2aDetectorMessenger(this);

  fLogicChamber = new G4LogicalVolume*[fNbOfChambers];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
B2aDetectorConstruction::~B2aDetectorConstruction()
{
  delete [] fLogicChamber; 
  delete fStepLimit;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4VPhysicalVolume* B2aDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();



  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4double chamber_x = 5*km;
  G4double chamber_y = 5*km;
  G4double chamber_z = 5*km;
  G4double chamberSpacing = 10*km; // from chamber center to center!

/*  G4double target_x = 5*km;
  G4double target_y = 5*km; 
  G4double target_z = 5*km;*/
 
/*  G4double tracker_x = 5*km;  
  G4double tracker_y = ((10+1)*chamberSpacing)/2; //fNbofChmbers isn't working here don't know why
  G4double tracker_z = 5*km;*/

  G4double world_x = 5*km;
  G4double world_y = 1.1 * ((10+1)*chamberSpacing)/2;
  G4double world_z = 5*km;
  G4double firstPosition = -((10+1)*chamberSpacing)/2 + chamberSpacing;

void B2aDetectorConstruction::DefineMaterials()
{
  // Material definition 

  G4NistManager* nistManager = G4NistManager::Instance();

  // Air defined using NIST Manager
  nistManager->FindOrBuildMaterial("G4_AIR");

  G4double density_air = 1.205*mg/cm3;
  G4Material* air  = G4Material::GetMaterial("G4_AIR");
  G4double scale_h = 8.4*km;

  for (G4int i = 0; i < fNbOfChambers; i++){
    G4double height = (i+1)*chamberSpacing; //remember the bottom half of world is -ve but this has to be +ve
    std::stringstream ss;
    ss << "name_" << i;
    G4String name = ss.str();
    new G4Material(name,(exp(-height/scale_h))*density_air,air);
    nistManager->FindOrBuildMaterial(name);
  }

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B2aDetectorConstruction::DefineVolumes()
{
  G4Material* air  = G4Material::GetMaterial("G4_AIR");


  G4Box* worldBox
    = new G4Box("World", world_x, world_y, world_z); //its size
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                          worldBox,   //its solid
                          air,      //its material
                          "World"); //its name
  
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,               // no rotation
                 G4ThreeVector(), // at (0,0,0)
                 worldLV,         // its logical volume
                 "World",         // its name
                 0,               // its mother  volume
                 false,           // no boolean operations
                 0,               // copy number
                 fCheckOverlaps); // checking overlaps 

  // Target
  
 // G4ThreeVector target_pos = G4ThreeVector(0,target_y-world_y,0); //should be at bottom just touch bottom edge 

  /*G4Box* targetBox
    = new G4Box("Target",target_x,target_y,target_z);

  fLogicTarget // the pointer bit is in header file, so it can be accesses form other .cc files
    = new G4LogicalVolume(targetBox, fTargetMaterial,"Target",0,0,0);

  new G4PVPlacement(0,               // no rotation
                    target_pos,  // at (x,y,z)
                    fLogicTarget,    // its logical volume
                    "Target",        // its name
                    worldLV,         // its mother volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps */


/*  G4Box* trackerBox
    = new G4Box("Tracker",tracker_x,tracker_y,tracker_z);
  G4LogicalVolume* trackerLV
    = new G4LogicalVolume(trackerBox, air, "Tracker",0,0,0);  
  new G4PVPlacement(0,               // no rotation
                    G4ThreeVector(), // at (x,y,z)
                    trackerLV,       // its logical volume
                    "Tracker",       // its name
                    worldLV,         // its mother  volume
                    false,           // no boolean operations
                    0,               // copy number
                    fCheckOverlaps); // checking overlaps */

  // Visualization attributes

  G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));

  worldLV      ->SetVisAttributes(boxVisAtt);
 // fLogicTarget ->SetVisAttributes(boxVisAtt);
//  trackerLV    ->SetVisAttributes(boxVisAtt);



  for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {

      G4double Yposition = firstPosition + copyNo * chamberSpacing;
//      G4double rmax =  rmaxFirst + copyNo * rmaxIncr;

      G4Box* chamberBox
        = new G4Box("Chamber_solid", chamber_x, chamber_y, chamber_z);

      fLogicChamber[copyNo] =
              new G4LogicalVolume(chamberBox,fChamberMaterial,"Chamber_LV",0,0,0);

      fLogicChamber[copyNo]->SetVisAttributes(chamberVisAtt);

      new G4PVPlacement(0,                            // no rotation
                        G4ThreeVector(0,Yposition,0), // at (x,y,z)
                        fLogicChamber[copyNo],        // its logical volume
                        "Chamber_PV",                 // its name
                        worldLV,                    // its mother  volume
                        false,                        // no boolean operations
                        copyNo,                       // copy number
                        fCheckOverlaps);              // checking overlaps 

  }

  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  G4double maxStep = 0.5*chamber_x;
  fStepLimit = new G4UserLimits(maxStep);
//  trackerLV->SetUserLimits(fStepLimit);
 
  // Always return the physical world

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2aDetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

/*  G4String trackerChamberSDname = "B2/TrackerChamberSD";
  B2TrackerSD* aTrackerSD = new B2TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  // Setting aTrackerSD to all logical volumes with the same name 
  // of "Chamber_LV".
  SetSensitiveDetector("Chamber_LV", aTrackerSD, true);

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void B2aDetectorConstruction::SetTargetMaterial(G4String materialName)
{
/*  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial = 
              nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fTargetMaterial = pttoMaterial;
        if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
        G4cout 
          << G4endl 
          << "----> The target is made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetTargetMaterial : "
          << materialName << " not found" << G4endl;
     }
  }*/
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetChamberMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fChamberMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fChamberMaterial = pttoMaterial;
        for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {
            if (fLogicChamber[copyNo]) fLogicChamber[copyNo]->
                                               SetMaterial(fChamberMaterial);
        }
        G4cout 
          << G4endl 
          << "----> The chambers are made of " << materialName << G4endl;
     } else {
        G4cout 
          << G4endl 
          << "-->  WARNING from SetChamberMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2aDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
  fCheckOverlaps = checkOverlaps;
}  
