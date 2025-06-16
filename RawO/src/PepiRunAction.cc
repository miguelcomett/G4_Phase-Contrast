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
// $Id: PepiRunAction.cc 94031 2015-11-05 11:54:38Z ihrivnac $
//
/// \file PepiRunAction.cc
/// \brief Implementation of the PepiRunAction class

#include "PepiRunAction.hh"
#include "PepiRun.hh"
#include "PepiEventAction.hh"
#include "G4LogicalVolume.hh"
#include "G4AccumulableManager.hh"
#include "G4Run.hh"
#include "G4LogicalVolumeStore.hh"
#include "PepiDetectorConstruction.hh"
#include "PepiPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4THitsMap.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <stdio.h>
#include <vector>
#include <sys/time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PepiRunAction::PepiRunAction(G4bool /*isOnMaster*/)
: G4UserRunAction(),
  fSeed(-1),        // RANLUX seed
  fLuxury(3),       // RANLUX luxury level (3 is default)
  fSaveRndm(-1),
  dose(0)
//  fObjMGDose(0)
{   
  fSDName.push_back(G4String("PixiRadSD"));   
  fSDName.push_back(G4String("IonChamberSD"));
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;
  const G4double picogray  = 1.e-12*gray;

  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PepiRunAction::~PepiRunAction()
{ 
  fSDName.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* PepiRunAction::GenerateRun()
{
  // - Generate new RUN object, which is specially
  // dedicated for MultiFunctionalDetector scheme.
  // Detailed description can be found in PepiRun.hh/cc.
  return new PepiRun(fSDName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiRunAction::BeginOfRunAction(const G4Run* aRun)
{ 

  if (!IsMaster()) //it is a slave, do nothing else
  {
    G4cout << ">>> Run " << std::setw(4) << aRun->GetRunID() << " starts on slave." << G4endl;
    return;
  }

  // - Each time a NEW RUN starts, reseed the random numbers engine
  // with a seed taken from the time in seconds from the Unix Epoch...
  fSeed = time(NULL);
  G4Random::setTheSeed(fSeed,fLuxury);
  // G4Random::setTheSeed(123456789,fLuxury);
  // G4Random::showEngineStatus();
      
  // - Save Rndm status
  if (fSaveRndm > 0) G4Random::saveEngineStatus("beginOfRun.rndm");

  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();
//  PepiEventAction* pepiEventAction = (PepiEventAction*)(G4RunManager::GetRunManager()->GetUserEventAction());
//  pepiEventAction->ResetDose();

  G4cout << ">>> Frame " << std::setw(4) << aRun->GetRunID() + 1 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PepiRunAction::EndOfRunAction(const G4Run* aRun)
{
  if (!IsMaster()) return;

  PepiRun* pepiRun = (PepiRun*)aRun;
  
  // - Get accumulated quantities for this RUN.
  fRunMaps = pepiRun->GetHCofThisRun();

  if (fSaveRndm > 0) G4Random::saveEngineStatus("endOfRun.rndm");
   // Merge accumulables
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
//  G4double edep2 = fEdep2.GetValue();

//  G4double rms = edep2 - edep*edep/nofEvents;
//  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;
  const auto detConstruction = static_cast<const PepiDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detConstruction->GetScoringVolume()->GetMass();
  dose = edep/mass;
//  G4double rmsDose = rms/mass;
//  G4cout
//     << " Cumulated dose per run, in scoring volume : "
//     << G4BestUnit(dose,"Dose")
//     << G4endl
//     << "------------------------------------------------------------"
//     << G4endl
//     << G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PepiRunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}
