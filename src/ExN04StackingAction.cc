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
// $Id: ExN04StackingAction.cc 78055 2013-12-03 08:27:48Z gcosmo $
//
/// \file parallel/ParN04/src/ExN04StackingAction.cc
/// \brief Implementation of the ExN04StackingAction class
//

#include "ExN04StackingAction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

ExN04StackingAction::ExN04StackingAction()
{ 
}

ExN04StackingAction::~ExN04StackingAction()
{
}

G4ClassificationOfNewTrack 
ExN04StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
    const G4ParticleDefinition *particle = aTrack->GetParticleDefinition();
    return ("neutron" == particle->GetParticleName() ||
            "proton" == particle->GetParticleName())
        ? fUrgent : fKill;
}

void ExN04StackingAction::NewStage()
{
}
    
void ExN04StackingAction::PrepareNewEvent()
{ 
}


