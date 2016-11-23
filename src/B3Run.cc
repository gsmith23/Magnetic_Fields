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
// $Id: B3Run.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file B3Run.cc
/// \brief Implementation of the B3Run class

#include "B3Run.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3Run::B3Run()
 : G4Run(), 
   fCollID_cryst(-1),
   fCollID_patient(-1),
   fPrintModulo(10000),
   fGoodEvents(0),
   fSumDose(0.)
{ }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3Run::~B3Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Called at the end of each event
void B3Run::RecordEvent(const G4Event* event)
{
  // Retrieve the collectionID corresponding to hits in the crystal
  // and in the patient. This is done only at the first event
  if ( fCollID_cryst < 0 ) {
   fCollID_cryst 
     = G4SDManager::GetSDMpointer()->GetCollectionID("crystal/edep");
   //G4cout << " fCollID_cryst: " << fCollID_cryst << G4endl;   
  }

  if ( fCollID_patient < 0 ) {
   fCollID_patient 
     = G4SDManager::GetSDMpointer()->GetCollectionID("patient/dose");
   //G4cout << " fCollID_patient: " << fCollID_patient << G4endl;   
  }

  G4int evtNb = event->GetEventID();
  
  if (evtNb%fPrintModulo == 0) { 
    G4cout << "\n---> end of event: " << evtNb << G4endl;
  }      
  
  //Hits collections
  //  
  // Get all hits-collections available for this events: there should be two 
  // hits-collection, one of hits in the patient and one of hits in the 
  // crystals. They are created in the UserGeometry.
  G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  if(!HCE) return;
               
  //Energy in crystals : identify 'good events'
  //
  const G4double eThreshold = 500*keV;
  G4int nbOfFired = 0;
   
  //ok, let's start the game: retrieve the hits-collection in the crystals.
  //This comes from a Geant4 multiscorer of type "G4PSEnergyDeposit", which scores 
  //energy deposit.
  G4THitsMap<G4double>* evtMap = 
    static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_cryst));
               
  std::map<G4int,G4double*>::iterator itr;
  //loop over the hits of the hit collection: in principle many crystals can 
  //be fired.
  //Notice: we are not interested in what's the energy deposited in one crystal 
  //or in the full array, but only in how many crystals the deposition is above
  //500 keV.

  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    G4double edep = *(itr->second);
    //count how many events have energy deposition above threshold
    if (edep > eThreshold) nbOfFired++;

    //these are the ID's of the detectors fired. 
    ///G4int copyNb  = (itr->first);

    ///G4cout << "\n  cryst" << copyNb << ": " << edep/keV << " keV ";
  }  
  if (nbOfFired == 2) fGoodEvents++;
  
  
  //Dose deposit in patient
  //
  G4double dose = 0.;

  //now retrieve the hits-collection in the patient
  //This comes from a Geant4 multiscorer of type "G4PSDoseDeposit" which 
  //scores total dose
  evtMap = static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_patient));
               
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    //the loop has only one execution, since the first value of the map is the 
    //copy number. Only one patient is available here.
    ///G4int copyNb  = (itr->first);
    dose = *(itr->second);
  }

  //Sum the dose delivered in this event to the cumulative run value
  fSumDose += dose;
  
  G4Run::RecordEvent(event);      
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3Run::Merge(const G4Run* aRun)
{
  const B3Run* localRun = static_cast<const B3Run*>(aRun);

  //do the merge of the information coming from the individual threads.
  fGoodEvents += localRun->fGoodEvents;
  fSumDose    += localRun->fSumDose;

  G4Run::Merge(aRun); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
