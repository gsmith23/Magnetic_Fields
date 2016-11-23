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
// $Id: B3DetectorConstruction.hh 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file B3DetectorConstruction.hh
/// @brief Definition of the B3DetectorConstruction class (Mandatory)

#ifndef B3DetectorConstruction_h
#define B3DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials (with their physical properties) and detector geometry.
///
/// This example simulates schematically a Positron Emitted Tomography system (PET).
/// A F18 source, marking patient's cancer cells, decays emitting a positron. 
/// The detection of gamma's, emitted back-to-back in the subsequent positron-electron annihilation, 
/// is performed by scintillating crystals, positioned to form a ring, with an appropriate rotation matrix. 
/// Several copies of Ring are placed in the full detector.

class B3DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    /// constructor
    B3DetectorConstruction();
    /// destructor
    virtual ~B3DetectorConstruction();

  public:
    /// Defines the detector geometry and returns a pointer to the physical World Volume
    virtual G4VPhysicalVolume* Construct();
    /// Register some of the detector's volumes as "sensitive"
    virtual void ConstructSDandField();
               
  private:
    /// Defines all the materials the detector is made of.
    void DefineMaterials();

    G4bool  fCheckOverlaps;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

