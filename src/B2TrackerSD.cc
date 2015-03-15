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


  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();

  //if (edep==0.) return false;


  G4String name = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
  // if a particle is unseen so far -- register
  if (particles.find(name) == particles.end())
    particles[name] = 1;
  // otherwise increase the amount of times we have seen it
  else
    ++particles[name];


  if ("neutron" == name) {
    // only store hits on target
    if ("Target" == aStep->GetTrack()->GetLogicalVolumeAtVertex()->GetName()) {
      B2TrackerHit* newHit = new B2TrackerHit();

      newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
      newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
                                                   ->GetCopyNumber());
      newHit->SetEdep(edep);
      newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());
      newHit ->SetKineticEnergy (aStep->GetTrack()->GetKineticEnergy()); //G4double
      newHit->SetMomentumDirection (aStep->GetTrack()->GetMomentumDirection()); // G4ThreeVector
      newHit->SetMomentum (aStep->GetTrack()->GetMomentum()); // G4ThreeVector
      newHit->SetVelocity (aStep->GetTrack()->GetVelocity()); // G4double
      newHit->SetTotalEnergy(aStep->GetTrack()->GetTotalEnergy()); //G4double
      // G4cout << name
      //     << " " << aStep->GetTrack()->GetParticleDefinition()->GetInstanceID()
      //     << " " << aStep->GetTrack()->GetKineticEnergy() << G4endl;
      fHitsCollection->insert( newHit );
    }
  } else if ("proton" != name) {
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    }




  // secondaries
  // const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep(); // std::vector<const G4Track*>*
  // std::vector<const G4Track*>::const_iterator it = secondaries->begin();
  // for (; it != secondaries->end(); ++it) {
  //   if ((*it)->GetParticleDefinition()->GetParticleName() == "neutron") {
  //     G4cout << "FOR FUCKS SAKE " 
  //         << (*it)->GetKineticEnergy() 
  //         << " " << (*it)->GetParticleDefinition()->GetPDGStable() 
  //         << " " << (*it)->GetParticleDefinition()->IsShortLived() 
  //         << " " << (*it)->GetParticleDefinition()->GetParticleType() 
  //         << G4endl;
  //   }
  // }

  

  //newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int nofHits = fHitsCollection->entries();
  G4cout << G4endl
        << "-------->Hits Collection: in this event they are " << nofHits 
        << " hits in the tracker chambers: " << G4endl;
  std::ofstream fw("hits.txt"); // open a file
  for ( G4int i=0; i<nofHits; i++ ) {
    fw << ((B2TrackerHit *)(*fHitsCollection)[i])->ToString() << std::endl;
  }
  fw.close();
  // iterator of a map (string -> int), used to loop through saved particle counts
  std::map<G4String, int>::iterator it = particles.begin();
  for (; it != particles.end(); ++it) {
    G4cout << it->first << ": " << it->second << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
