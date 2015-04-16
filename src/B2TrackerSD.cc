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
// $Id: B2TrackerSD.cc 87359 2014-12-01 16:04:27Z gcosmo $
//
/// \file B2TrackerSD.cc
/// \brief Implementation of the B2TrackerSD class

#include "B2TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4ProcessType.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerSD::B2TrackerSD(const G4String& name,
                         const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerSD::~B2TrackerSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new B2TrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B2TrackerSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{  
  G4Track *aTrack = aStep->GetTrack();
  G4StepPoint *preStep = aStep->GetPreStepPoint();

  G4String name = aTrack->GetParticleDefinition()->GetParticleName();

  // store only neutrons
  if ("neutron" == name) {
    // only store hits entering Target
    if ("Target" == aTrack->GetNextVolume()->GetName() &&
        NULL != preStep->GetProcessDefinedStep() &&
        fTransportation == preStep->GetProcessDefinedStep()->GetProcessType()) {

      B2TrackerHit* newHit = new B2TrackerHit();

      newHit->SetTrackID            (aTrack->GetTrackID());
      newHit->SetChamberNb          (preStep->GetTouchableHandle()->GetCopyNumber());
      newHit->SetEdep               (aStep->GetTotalEnergyDeposit());
      newHit->SetPos                (aStep->GetPostStepPoint()->GetPosition());
      newHit->SetKineticEnergy      (aTrack->GetKineticEnergy());
      newHit->SetMomentumDirection  (aTrack->GetMomentumDirection());
      newHit->SetMomentum           (aTrack->GetMomentum());
      newHit->SetVelocity           (aTrack->GetVelocity());
      newHit->SetTotalEnergy        (aTrack->GetTotalEnergy());

      fHitsCollection->insert(newHit);
    }
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int nofHits = fHitsCollection->entries();
  G4cout << G4endl
        << "-------->Hits Collection: in this event they are " << nofHits 
        << " hits in the tracker chambers: " << G4endl;

  std::ofstream fw("hits.txt",std::ofstream::app); // open a file
  for ( G4int i=0; i<nofHits; i++ ) {
    fw << ((B2TrackerHit *)(*fHitsCollection)[i])->ToString() << std::endl;
  }
  fw << std::endl;
  fw.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
