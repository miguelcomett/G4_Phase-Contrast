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
// $Id: PepiRunAction.hh 94031 2015-11-05 11:54:38Z ihrivnac $
//
/// \file PepiRunAction.hh
/// \brief Definition of the PepiRunAction class

#ifndef PepiRunAction_h
#define PepiRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4THitsMap.hh"
#include "globals.hh"
#include <vector>
#include "G4Accumulable.hh"

class G4Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PepiRunAction : public G4UserRunAction
{ 
  public:

    PepiRunAction(G4bool isOnMaster);
    virtual ~PepiRunAction();

  public:
    
    virtual G4Run*  GenerateRun();
    virtual void    BeginOfRunAction(const G4Run*);
    virtual void    EndOfRunAction(const G4Run*);
    void AddEdep (G4double edep);
    // - Get all HitsMaps
    std::vector<G4THitsMap<G4double>*> GetHCofThisRun() const {return fRunMaps;};
//    G4double GetObjectDose(G4String);
    G4double GetObjectDose() const {return  dose;};
  private:

    G4long  fSeed;
    G4int   fLuxury;
    G4int   fSaveRndm;

    std::vector<G4THitsMap<G4double>*> fRunMaps;

   G4double dose;
//    G4double fObjMGDose;
    G4Accumulable<G4double> fEdep=0.;
    G4Accumulable<G4double> fEdep2=0.;

    std::vector<G4String> fSDName; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
