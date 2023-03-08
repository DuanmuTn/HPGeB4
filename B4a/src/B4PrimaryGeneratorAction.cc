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
/// \file B4PrimaryGeneratorAction.cc
/// \brief Implementation of the B4PrimaryGeneratorAction class

#include "B4PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::B4PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(nullptr)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  //
  auto particleDefinition 
    = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  G4double Mx=0;
  G4double My=0;
  G4double Mz=0;

  fParticleGun->SetParticleEnergy(0.477*MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume 
  // from G4LogicalVolumeStore
  //
/*  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();  
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  } 
 */ 
  // Set gun position
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;

  G4ParticleDefinition* particle
                    = particleTable->FindParticle(particleName="gamma");
  fParticleGun->SetParticleDefinition(particle);
  
const double pi=3.141592653;
//set random beam direction
  G4double theta=G4UniformRand()*2*pi;
  G4double  phi =G4UniformRand()*2*pi;
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)));
 Mx+=sin(theta)*cos(phi);
 My+=sin(theta)*sin(phi);
 Mz+=cos(theta);
  //set particle random energy
//G4double position = -0.5*(fDetector->GetWorldSizeX());
//
//set source position
  G4double distance=0.*cm;//distance of source and detector(cm)
  G4double radius=6.*mm;
  G4double Cx = 6.*cm;
  G4double Cy = 4.*cm;
  G4double Cz = 1.*cm;
  G4int rand=10*G4UniformRand();
  /*
  //here defines the collected cup
  G4double position = distance+Cz/2+pow(-1,rand%2)*G4UniformRand()*Cz/2;
  G4double pos1=pow(-1,rand%2)*(Cx/2)*G4UniformRand()*mm;
  G4double pos2=pow(-1,rand%2)*(Cy/2)*G4UniformRand()*mm;
 //  */
//  /*
  
  G4double position = distance;//+pow(-1,rand%2)*0.01*G4UniformRand()*cm;
  G4double pos1=pow(-1,rand%2)*radius*G4UniformRand()*mm;
  G4double pos2=pow(-1,rand%2)*radius*G4UniformRand()*mm;
//  */
  fParticleGun->SetParticlePosition(G4ThreeVector(pos1,pos2,-1.*position));

  fParticleGun->GeneratePrimaryVertex(anEvent);
  G4int n=anEvent->GetEventID();
  if(n%100==0){
  G4cout<<"momentum is "<<sin(theta)*cos(phi)<<"  "<<sin(theta)*sin(phi)<<"  "<<cos(theta)<<G4endl;
  G4cout<<"position is "<<"  "<<pos1<<"  "<<pos2<<"  "<<position<<G4endl;
G4cout<<"the sum of momentum is :"<<Mx<<" "<<My<<" "<<Mz<<G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

