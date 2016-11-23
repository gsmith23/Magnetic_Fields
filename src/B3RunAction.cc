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
// $Id: B3RunAction.cc 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file B3RunAction.cc
/// \brief Implementation of the B3RunAction class

#include "B3RunAction.hh"
#include "B3PrimaryGeneratorAction.hh"
#include "B3Run.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3RunAction::B3RunAction()
 : G4UserRunAction()
{  
  //add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3RunAction::~B3RunAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* B3RunAction::GenerateRun()
{ return new B3Run; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3RunAction::BeginOfRunAction(const G4Run* run)
{ 
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3RunAction::EndOfRunAction(const G4Run* run)
{
  //retrieve the number of events produced in the run
  G4int nofEvents = run->GetNumberOfEvent();

  //do nothing, if no events were processed
  if (nofEvents == 0) return;
  
  // Run conditions
  // This retrieves the UserPrimaryGeneratorAction object: it is retrieved through the 
  // G4RunManager. 
  //
  // Following the B3ActionInitialization, the UserPrimaryGeneratorAction 
  // exists for all workers (-> Build()) but not for the master (-> BuildForMaster()). The 
  // B3RunAction instead exists for the master and for the worker. So, when the function is 
  // executed by the master, no pointer is found for the primary generator and 
  // generatorAction = NULL

  const B3PrimaryGeneratorAction* generatorAction = static_cast<const B3PrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
 

  G4String partName;
  if (generatorAction) 
  {
      // The GetParticleGun() is defined inside the B3PrimaryGeneratorAction.hh file and
      // it returns a pointer to the current concrete implementation of the G4VPrimaryGeneratorAction
      // Note that GetParticleDefinition return a 'const' type in the case of the use of the G4ParticleGun but not
      // in the case of the G4ParticleGeneralSource. Hence the GetParticleGun in the B3PrimaryGeneratorAction.hh
      // must not be defined 'const' in the case of the GeneralParticleSource
      //
    G4ParticleDefinition* particle = generatorAction -> GetParticleGun() -> GetParticleDefinition();
    partName = particle->GetParticleName();
  }  
  
  //results
  //
  const B3Run* b3Run = static_cast<const B3Run*>(run);
  G4int nbGoodEvents = b3Run->GetNbGoodEvents();
  G4double sumDose   = b3Run->GetSumDose();                           
        
  //print
  //
  if (IsMaster())
  {
    G4cout
     << "\n--------------------End of Global Run-----------------------"
     << " \n The run was " << nofEvents << " events ";
  }
  else
  {
    G4cout
     << "\n--------------------End of Local Run------------------------"
     << " \n The run was " << nofEvents << " "<< partName;
  }      
  G4cout
     << "; Nb of 'good' e+ annihilations: " << nbGoodEvents
     << "\n Total dose in patient : " << G4BestUnit(sumDose,"Dose")   
     << "\n------------------------------------------------------------\n"
     << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
