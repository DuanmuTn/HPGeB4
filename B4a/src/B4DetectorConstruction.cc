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
// 
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fHPGePV(nullptr),
   fFingerPV(nullptr),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_Ge");
  nistManager->FindOrBuildMaterial("G4_Al");
  nistManager->FindOrBuildMaterial("G4_Cu");
  
  // Liquid Nitrogen material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  new G4Material("liquidNitrogen", z=7., a= 28.01*g/mole, density= 0.81*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
//  G4int nofLayers = 10;
//  G4double absoThickness = 10.*mm;
//  G4double fingerThickness =  5.*mm;
//  G4double calorSizeXY  = 10.*cm;
	G4double HPGeR = 3.145*cm;
	G4double HPGeh = 5.42*cm;
	G4double foilH = 0.01*mm;
	G4double foilR = 4.*cm;
	G4double shellH = 5*mm;
	G4double shellR = 4.*cm;
	G4double deadgap = 0.9*cm;
	G4double fingerH = 4.14*cm;
	G4double fingerR = 4.35*mm;
	G4double HPGeH = HPGeh-deadgap;
	G4double startAngle = 0.*deg;
	G4double spanAngle = 360.*deg;
	G4double R = 3.*cm;
	G4double foil2H = 0.2*mm;
	G4double foil2R = 3.*mm;
	G4double foil2dis = 11.9*cm;

//  auto layerThickness = absoThickness + fingerThickness;
//  auto calorThickness = nofLayers * layerThickness;
  auto worldSizeXY = 8.*cm;//3. * HPGeR;
  auto worldSizeZ  = 30.*cm;//3. * HPGeH; 
  auto HPGedistance =foilH+shellH+deadgap;  
  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto WallMaterial = G4Material::GetMaterial("G4_Pb");
  auto HPGeMaterial = G4Material::GetMaterial("G4_Ge");
  auto foilMaterial = G4Material::GetMaterial("G4_Al");
  auto shellMaterial = G4Material::GetMaterial("G4_Al");
  auto fingerMaterial = G4Material::GetMaterial("G4_Cu");
  auto SourceMaterial = G4Material::GetMaterial("Galactic");

  
  if ( ! defaultMaterial || ! WallMaterial || ! HPGeMaterial || !  foilMaterial ||! fingerMaterial ||! shellMaterial) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //
  G4Tubs* worldS 
    = new G4Tubs("World",           // its name
                 0.*mm, worldSizeXY, worldSizeZ/2, startAngle, spanAngle); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Calorimeter
/*  //  
  auto calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size
                         
  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
*/  //                                 
  // Layer
  //
/*  auto layerS 
    = new G4Box("Layer",           // its name
                 calorSizeXY/2, calorSizeXY/2, layerThickness/2); // its size
                         
  auto layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

  new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 calorLV,          // its mother
                 kZAxis,           // axis of replication
                 nofLayers,        // number of replica
                 layerThickness);  // witdth of replica
  
*/  //                               
  // Absorber
  //
/*  auto BallS 
    = new G4Box("BAll",            // its name
                 3./2*cm, 3./2*cm, 3./2*cm); // its size
                         
  auto BallLV
    = new G4LogicalVolume(
                 BallS,        // its solid
                 HPGeMaterial, // its material
                 "BALL");          // its name
                                   
  fBallPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,-3./2*cm), // its position
                 BallLV,       // its logical volume                         
                 "BALL",           // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
*/
  auto HPGeS 
    = new G4Tubs("HpGe",            // its name
                 fingerR, HPGeR, HPGeH/2, startAngle, spanAngle); // its size
                         
  auto HPGeLV
    = new G4LogicalVolume(
                 HPGeS,        // its solid
                 HPGeMaterial, // its material
                 "HpGe");          // its name
                                   
  fHPGePV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,HPGedistance + HPGeH/2), // its position
                 HPGeLV,       // its logical volume                         
                 "HpGe",           // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  auto foilS 
    = new G4Tubs("Foil",            // its name
                 0.*mm, foilR, foilH/2, startAngle, spanAngle); // its size
                         
  auto foilLV
    = new G4LogicalVolume(
                 foilS,        // its solid
                 foilMaterial, // its material
                 "Foil");          // its name
                                   
  ffoilPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,foilH/2), // its position
                 foilLV,       // its logical volume                         
                 "Foil",           // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  auto foil2S 
    = new G4Tubs("Foil2",            // its name
                 0.*mm, foil2R, foil2H/2, startAngle, spanAngle); // its size
                         
  auto foil2LV
    = new G4LogicalVolume(
                 foil2S,        // its solid
                 foilMaterial, // its material
                 "Foil2");          // its name
                                   
  ffoil2PV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,(-1)*foil2dis), // its position
                 foil2LV,       // its logical volume                         
                 "Foil2",           // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  auto shellS 
    = new G4Tubs("Shell",            // its name
                 0.*mm, shellR, shellH/2, startAngle, spanAngle); // its size
                         
  auto shellLV
    = new G4LogicalVolume(
                 shellS,        // its solid
                 shellMaterial, // its material
                 "Shell");          // its name
                                   
  fshellPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,foilH + shellH/2), // its position
                 shellLV,       // its logical volume                         
                 "Shell",           // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  auto deadgapS 
    = new G4Tubs("Deadgap",            // its name
                 0.*mm, HPGeR, deadgap/2, startAngle, spanAngle); // its size
                         
  auto deadgapLV
    = new G4LogicalVolume(
                 deadgapS,        // its solid
                 HPGeMaterial, // its material
                 "Deadgap");          // its name
                                   
  fdeadgapPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0.,foilH+shellH+deadgap/2), // its position
                 deadgapLV,       // its logical volume                         
                 "Deadgap",           // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  //                              
  // Finger
  //
  auto fingerS 
    = new G4Tubs("Finger",             // its name
                 0.*mm, fingerR, fingerH/2, startAngle, spanAngle); // its size
                         
  auto fingerLV
    = new G4LogicalVolume(
                 fingerS,             // its solid
                 fingerMaterial,      // its material
                 "Finger");           // its name
                                   
  fFingerPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., (HPGeH-fingerH/2)+HPGedistance), // its position
                 fingerLV,            // its logical volume                         
                 "Finger",            // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  auto GefingerS 
    = new G4Tubs("GeFinger",             // its name
                 0.*mm, fingerR, (HPGeH-fingerH)/2, startAngle, spanAngle); // its size
                         
  auto GefingerLV
    = new G4LogicalVolume(
                 GefingerS,             // its solid
                 HPGeMaterial,      // its material
                 "GeFinger");           // its name
                                   
  fGeFingerPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., (HPGeH-fingerH)/2+HPGedistance), // its position
                 GefingerLV,            // its logical volume                         
                 "GeFinger",            // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
 /* 
   G4double Sx = 6.*cm;
   G4double Sy = 4.*cm;
   G4double Sz = 1.*cm;
   G4double distance = 0.*cm;
   G4double p_pos_z = distance+Sz/2;
   G4ThreeVector pos_Source(0,0, -1.*p_pos_z);
   G4Box* SourceTubsS 
      = new G4Box("Source", Sx/2, Sy/2, Sz/2); 
   G4LogicalVolume* SourceLV
      = new G4LogicalVolume(SourceTubsS, 
	  						foilMaterial,
							"Source");
   fSourcePV
      = new G4PVPlacement(0,                 //no rotation
                          pos_Source,          //at position
			  "Source",            //its name
			  SourceLV,           //its logical volume
			  worldPV,          //its mother phyiscal volume
			  false,             //no boolean operation
			  2);                //copy number
*/
// /*
//Creating the particle source
	G4double p_pos_x=0.*cm;
	G4double p_pos_y=0.*cm;
	G4double p_ang_min=0.*deg;
	G4double p_ang_max=360.*deg;
	G4double source_distance=0.*cm;
	G4double sourceH=0.1*mm;
	G4double sourceRmin=0.*mm;
	G4double sourceRmax=6.*mm;



G4double ang_Sourcex = 0.*deg;
G4double ang_Sourcey = 0.*deg;
G4double ang_Sourcez = 0.*deg;
   G4RotationMatrix* ang_Source=0;
   ang_Source=new G4RotationMatrix();
   ang_Source->rotateX(ang_Sourcex);
   ang_Source->rotateY(ang_Sourcey);
   ang_Source->rotateZ(ang_Sourcez);
	G4double p_pos_z=(sourceH/2+source_distance)*mm;
   G4ThreeVector pos_Source(p_pos_x, p_pos_y, -1.*p_pos_z);
   G4Tubs* SourceTubsS 
      = new G4Tubs("Source", sourceRmin, sourceRmax, sourceH, p_ang_min, p_ang_max); 
     
   G4LogicalVolume* SourceLV
      = new G4LogicalVolume(SourceTubsS, 
	  						SourceMaterial,
							"Source");
   fSourcePV
      = new G4PVPlacement(ang_Source,                 //no rotation
                          pos_Source,          //at position
			  "Source",            //its name
			  SourceLV,           //its logical volume
			  worldPV,          //its mother phyiscal volume
			  false,             //no boolean operation
			  2);                //copy number
   

  //Creating source body(Cu back)
  //
  	G4double Sbodyr=3.*cm;
	G4double SbodyH=2.*mm;

  	G4double Sbodypos=(sourceH+source_distance+SbodyH/2)*mm;
  auto SbodyMaterial = G4Material::GetMaterial("G4_Cu");
   G4Tubs* SbodyTubsS 
      = new G4Tubs("Sbody", 0.*mm, Sbodyr, SbodyH/2, p_ang_min, p_ang_max); 
     
   G4LogicalVolume* SbodyLV
      = new G4LogicalVolume(SbodyTubsS, 
	  						SbodyMaterial,
							"Sbody");
   G4ThreeVector pos_Sbody(0.*mm, 0.*mm, -1.*Sbodypos*mm);
   fSbodyPV
      = new G4PVPlacement(0,                 //no rotation
                          pos_Sbody,          //at position
			  "Sbody",            //its name
			  SbodyLV,           //its logical volume
			  worldPV,          //its mother phyiscal volume
			  false,             //no boolean operation
			  2);                //copy number
  // print parameters
  // */
/*  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << nofLayers << " layers of: [ "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
    << " + "
    << fingerThickness/mm << "mm of " << fingerMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
 */ 
  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  HPGeLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
