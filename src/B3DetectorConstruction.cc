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
// $Id: B3DetectorConstruction.cc 71323 2013-06-13 16:54:23Z gcosmo $
//
/// \file B3DetectorConstruction.cc
/// \brief Implementation of the B3DetectorConstruction class

#include "B3DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4GlobalMagFieldMessenger.hh"
///// Task 1c.1
// Add here the include for the Magnetic Field Messenger

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::B3DetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true)
{
  // **Material definition**
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::~B3DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::DefineMaterials()
{

  G4double a; // mass of a mole
  G4double z; // mean number of protons
  G4String name, symbol;

  // Definition of element Oxygen
  a=16.00*g/mole;
  G4Element*  elO = new G4Element(name="Oxygen", symbol="O", z=8., a);
  
  a=28.09*g/mole;
  G4Element*  elSi = new G4Element(name="Silicon", symbol="Si", z=14., a);

  a=174.97*g/mole;
  G4Element*  elLu = new G4Element(name="Lutetium", symbol="Lu", z=71., a);
  
  a=88.91*g/mole;
  G4Element*  elYt = new G4Element(name="Yttrium", symbol="Lu", z=39., a);
  
  G4double density;
  density = 7.4*g/cm3;

  // Declare the material with its density and number of components
  G4Material * LSO;
  LSO = new G4Material("Lu2SiO5",//name
		       density,//density
		       3);//number of elements
  
  LSO->AddElement(elLu,2);
  LSO->AddElement(elSi,1);
  LSO->AddElement(elO,5);

  
  density = 7.1*g/cm3;
  G4Material * LYSO;
  LYSO = new G4Material("LYSO",//name
			density,//density
			4);//number of elements
  
  LYSO->AddElement(elLu,18);
  LYSO->AddElement(elYt,2);
  LYSO->AddElement(elSi,10);
  LYSO->AddElement(elO,50);
  
  

  // Dump the Table of registered materials 
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B3DetectorConstruction::Construct()
{  
  //
  // First of all, fix the general parameters
  //
  // Gamma detector Parameters
  //
  // Crystals sizes
  // Notice that the crystals will be rotated to be placed
  // in the ring. What is called "cryst_dX" corresponds to 
  // to the z direction in the final system. The "cryst_dZ" is 
  // the "thickness" of the crystal in the ring
  //
  ///// Task 1b.1: 
  // Modify the present geometry by slightly changing 
  // cryst_dY  from 6. cm to 10. cm.
  G4double cryst_dX = 10*cm, cryst_dY = 10*cm, cryst_dZ = 3*cm;

  
  // Number of crystals per ring, and number of rings.
  ///// Task 1b.1: 
  // Change the crystals granularity in the detector, 
  // for example decrease the number of crystals per ring 
  // from 32 to 16
  G4int nb_cryst = 16;
  G4int nb_rings = 9;


  // The inner radius of the ring is accomodated such that  
  // the nb_cryst crystals do not overlap. The approach presented 
  // here is pratically equivalent to ring_R1*twopi = nb_cryst*cryst_dY,
  // i.e. the crystals fill entirely the inner circumference.
  // The outer radius is calculated accordingly, by adding the thickness 
  // of the crystals and an extra term 1/cos. The computation might 
  // screw up if there are too-few crystals.
  //
  G4double dPhi = twopi/nb_cryst, half_dPhi = 0.5*dPhi;
  G4double cosdPhi = std::cos(half_dPhi);
  G4double tandPhi = std::tan(half_dPhi);
  // 
  G4double ring_R1 = 0.5*cryst_dY/tandPhi;
  G4double ring_R2 = (ring_R1+cryst_dZ)/cosdPhi;
  //
  // total length of the detector: nb_rings, each having 
  // thickness of cryst_dX along the z-axis (in the final reference 
  // frame).
  G4double detector_dZ = nb_rings*cryst_dX;
  //
  G4NistManager* nist = G4NistManager::Instance();

  // **Retrieve Nist Materials** 
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  //G4Material* cryst_mat   = nist->FindOrBuildMaterial("Lu2SiO5");
  G4Material* cryst_mat   = nist->FindOrBuildMaterial("LYSO");
  //G4Material* cryst_mat   = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");

    
  //     
  // ***** World *****
  //
  G4double world_sizeXY = 1.0*m;
  G4double world_sizeZ  = 1.0*m;

  // Create the world volume as a box. This is big enough to contain            
  // all the detector rings.                                                    
 
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ); //its size
  
  // World Logical Volume definition
    
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        default_mat,         //its material
                        "World");            //its name
                      
  // World Physical Volume Placement at (0,0,0)              
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps 
                 
  //
  // ring
  //
  // define one ring as an "envelope" made of air. This will be filled 
  // by the crystals (dautgher volumes). The logical volume of the ring 
  // (which contains all daughters) will be then placed many times
  //
  G4Tubs* solidRing =
    new G4Tubs("Ring",      //name
	       ring_R1,      //inner radius
	       ring_R2,      //outer radius
	       0.5*cryst_dX, //height
	       0.,           //start angle
	       twopi);       //spanning angle
      
  G4LogicalVolume* logicRing =                         
    new G4LogicalVolume(solidRing,           //its solid
                        default_mat,         //its material
                        "Ring");             //its name
                    
  //     
  // define crystal
  //
  // define one crystal now. It is made a bit smaller than the 
  // real size, to allow for wrapping. The same logical volume is 
  // placed many times in its mother volume, to produce the actual 
  // ring of crystals.
  //
  G4double gap = 0.5*mm;        //a gap for wrapping
  G4double dX = cryst_dX - gap, dY = cryst_dY - gap;
  G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, cryst_dZ/2);
                     
  G4LogicalVolume* logicCryst = 
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name
               
  // place crystals within a ring: loop 
  //
  for (G4int icrys = 0; icrys < nb_cryst ; icrys++) {
    G4double phi = icrys*dPhi;
    // create a rotation matrix... (identity, by defauly)
    G4RotationMatrix rotm  = G4RotationMatrix();
    //... and apply rotations. Notice that all rotations are 
    // referred to the mother volume.
    rotm.rotateY(90*deg); 
    rotm.rotateZ(phi);
    // Calculate position with respect to the reference frame 
    // of the mother volume
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);     
    G4ThreeVector position = (ring_R1+0.5*cryst_dZ)*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);
                      
    // Place the crystal with the appropriate transformation
    new G4PVPlacement(transform,             //rotation,position
                      logicCryst,            //its logical volume
                      "crystal",             //its name
                      logicRing,             //its mother  volume
                      false,                 //no boolean operation
                      icrys,                 //copy number
                      fCheckOverlaps);       // checking overlaps 
  }
                                                      
  //
  // full detector. This is an "envelope" made of air which contains 
  // the rings. It is a big hollow cylinder, which will then host all rings 
  // as daughters (multiple placement). Notice that the cylinder is "hollow" 
  // so the central part will still belong to the mother (world) volume. 
  //
  G4Tubs* solidDetector =
    new G4Tubs("Detector", ring_R1, ring_R2, 0.5*detector_dZ, 0., twopi);
      
  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(solidDetector,       //its solid
                        default_mat,         //its material
                        "Detector");         //its name
                                 
  // 
  // place rings within detector 
  //
  G4double OG = -0.5*(detector_dZ + cryst_dX);
  for (G4int iring = 0; iring < nb_rings ; iring++) {
    OG += cryst_dX;
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,OG), //position
                      logicRing,             //its logical volume
                      "ring",                //its name
                      logicDetector,         //its mother  volume
                      false,                 //no boolean operation
                      iring,                 //copy number
                      fCheckOverlaps);       // checking overlaps 
  }
                       
  //
  // place detector in world, right in the centter
  //                    
//   new G4PVPlacement(0,                       //no rotation
//                     G4ThreeVector(),         //at (0,0,0)
//                     logicDetector,           //its logical volume
//                     "Detector",              //its name
//                     logicWorld,              //its mother  volume
//                     false,                   //no boolean operation
//                     0,                       //copy number
//                     fCheckOverlaps);         // checking overlaps 
  
  //
  // patient
  //
  // create a cylinder to simulate the brain of the patient. The 
  // G4_BRAIN_ICRP NIST material is used. The patient is placed at 
  // the center of the world.
  G4double patient_radius = 8*cm;
  G4double patient_dZ = 10*cm;  
  G4Material* patient_mat = nist->FindOrBuildMaterial("G4_BRAIN_ICRP");
    
  G4Tubs* solidPatient =
    new G4Tubs("Patient", 0., patient_radius, 0.5*patient_dZ, 0., twopi);
      
  G4LogicalVolume* logicPatient =                         
    new G4LogicalVolume(solidPatient,        //its solid
                        patient_mat,         //its material
                        "PatientLV");        //its name
               
  //
  // place patient in world
  //                    
//   new G4PVPlacement(0,                       //no rotation
//                     G4ThreeVector(),         //at (0,0,0)
//                     logicPatient,            //its logical volume
//                     "Patient",               //its name
//                     logicWorld,              //its mother  volume
//                     false,                   //no boolean operation
//                     0,                       //copy number
//                     fCheckOverlaps);         // checking overlaps 

  ///// Task 1b.2
  // Here create the solid for patient's body
  
  G4double body_dZ = 150*cm;
  G4Tubs*  solidBody = new G4Tubs("Body",                // the name
				  0,
				  12*cm,
				  0.5*body_dZ,// in z direction 
				  0.,         // in radians
				  twopi);     // in radians   

                             
  ///// Task 1b.3
  // Here recover the body-tissue material from NIST DB.
  G4Material* body_mat = nist->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP");
  
  // The create the body logical volume and set its visualization attributes
  // following intructions on the exercise text

  G4LogicalVolume* logicBody =                         
    new G4LogicalVolume(solidBody,        //its solid
                        body_mat,         //its material
                        "BodyLV");        //its name
  
  G4Colour red = G4Colour(1.0, 0.0, 0.0);
  
  logicBody->SetVisAttributes(new G4VisAttributes(red));
  
  ///// Task 1b.4 
  // Here place an instance of the body-logical volume in the world volume
//   new G4PVPlacement(0,                        //no rotation
//                     G4ThreeVector(0,0,-0.5*(body_dZ + patient_dZ + 0.002)),
//                     logicBody,                //its logical volume
//                     "Body",                   //its name
//                     logicWorld,               //its mother  volume
//                     false,                    //no boolean operation
//                     0,                        //copy number
//                     fCheckOverlaps);         // checking overlaps 
             
  // Visualization attributes
  // the two "envelopes" that were created to accomodate crystals and 
  // rings are set as invisible: they are made out of air and they are 
  // not real physical objects. They help to code the geometry in Geant4.
  logicRing->SetVisAttributes (G4VisAttributes::Invisible);
  logicDetector->SetVisAttributes (G4VisAttributes::Invisible);    
  
  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl; 

  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::ConstructSDandField()
{
  //
  // Register some of the volumes as "sensitive" and decide the 
  // type of sensitivity that they have
  //
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  
  // declare crystal as a MultiFunctionalDetector scorer
  //  
  // Create a new scorer (G4MultiFunctionalDetector) and set its 
  // "capability" to G4PSEnergyDeposit (will score total energy deposit)
  G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
  // Attach the scorer to the logical volume
  SetSensitiveDetector("CrystalLV",cryst);
  
  // declare patient as a MultiFunctionalDetector scorer
  //  
  // Create a new scorer (G4MultiFunctionalDetector) and set its 
  // "capability" to G4PSDoseDeposit (will score total dose)
  G4MultiFunctionalDetector* patient = new G4MultiFunctionalDetector("patient");
  G4VPrimitiveScorer* primitiv2 = new G4PSDoseDeposit("dose");
  patient->RegisterPrimitive(primitiv2);
  // Attach the scorer to the proper logical volume
  SetSensitiveDetector("PatientLV",patient);

  ///// Task 1c.1
  // Here is the place to set the magnetic field
  
  G4double Bz = 1 * tesla;

  G4GlobalMagFieldMessenger* fMagFieldMessenger = new G4GlobalMagFieldMessenger(G4ThreeVector(0,0,Bz));
  
  fMagFieldMessenger->SetVerboseLevel(1);
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
