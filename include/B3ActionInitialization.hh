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
// $Id: B3ActionInitialization.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file B3ActionInitialization.hh
/// \brief Definition of the B3ActionInitialization class (Mandatory)

#ifndef B3ActionInitialization_h
#define B3ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

/// Action initialization class.
///

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class B3ActionInitialization : public G4VUserActionInitialization
{
public:
  ///Constructor
  B3ActionInitialization();
  ///Destructor
  virtual ~B3ActionInitialization();

  ///Basic interface intherited from the base class G4VUserActionInitialization

  ///Registers User Actions for the master (only Run Action)
  virtual void BuildForMaster() const;

  ///Register User Actions for the workers
  virtual void Build() const;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
