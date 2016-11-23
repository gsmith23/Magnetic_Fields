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

// This class file describes how to generate the primary events

// The default particle beam is an ion (F18), at rest, randomly distributed
// within a zone inside a patient and is defined in
// B3PrimaryGeneratorAction::GeneratePrimaries().
// The type of a primary particle can be changed with G4ParticleGun commands
// (see run2.mac).


#include "B3PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// This following is the class contructor and
// set the default properties of the primary beam

B3PrimaryGeneratorAction::B3PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0)
{
 
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  
  // default particle kinematic
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle
    = particleTable->FindParticle("chargedgeantino");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,.45,0.));
  fParticleGun->SetParticleEnergy(1*eV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,-1.,0.));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
B3PrimaryGeneratorAction::~B3PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
  if (particle == G4ChargedGeantino::ChargedGeantino()) {
    //fluorine
    G4int Z = 9, A = 18;
    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.*keV;
    
    G4ParticleDefinition* ion
      = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
  }
  
  // fixed position  
  G4double x0  = 0*cm, y0  = 45.*cm, z0  = 0*cm;
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  
  //Create the vertex of the primary particle passing it to the G4Event pointer
  //
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

