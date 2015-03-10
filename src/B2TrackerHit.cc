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
// $Id: B2TrackerHit.cc 69706 2013-05-13 09:12:40Z gcosmo $
//
/// \file B2TrackerHit.cc
/// \brief Implementation of the B2TrackerHit class

#include "B2TrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>
#include <sstream>
#include <cstring>

G4ThreadLocal G4Allocator<B2TrackerHit>* B2TrackerHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerHit::B2TrackerHit()
 : G4VHit(),
   fTrackID(-1),
   fChamberNb(-1),
   fEdep(0.),
   fPos(G4ThreeVector()),
   fKineticEnergy(-1),
   fMomentumDirection(G4ThreeVector()),
   fVelocity(-1),
   fMomentum(G4ThreeVector()),
   fTotalEnergy(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerHit::~B2TrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2TrackerHit::B2TrackerHit(const B2TrackerHit& right)
  : G4VHit()
{
  fTrackID   = right.fTrackID;
  fChamberNb = right.fChamberNb;
  fEdep      = right.fEdep;
  fPos       = right.fPos;
  fMomentumDirection = right.fMomentumDirection;
  fMomentum = right.fMomentum;
  fVelocity = right.fVelocity;
  fTotalEnergy = right.fTotalEnergy;
  fKineticEnergy = right.fKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const B2TrackerHit& B2TrackerHit::operator=(const B2TrackerHit& right)
{
  fTrackID   = right.fTrackID;
  fChamberNb = right.fChamberNb;
  fEdep      = right.fEdep;
  fPos       = right.fPos;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int B2TrackerHit::operator==(const B2TrackerHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2TrackerHit::Print()
{
  G4cout
     << "  trackID: " << fTrackID 
     << " chamberNb: " << fChamberNb
     << "Edep: " << std::setw(7) << G4BestUnit(fEdep,"Energy")
     << " Position: " << std::setw(7) << G4BestUnit( fPos,"Length")
     << " KE: " << G4BestUnit( fKineticEnergy, "Energy")
     << " Momentum:" << fMomentum
     << " MomentumDirection:" << fMomentumDirection
     << " Velocity: " << fVelocity
     << " Total Energy: " << G4BestUnit(fTotalEnergy, "Energy")
     << G4endl;
}

//G4double MD_x = fMomentumDirection.x; 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
std::string B2TrackerHit::ToString() // need to sort out vectors
{
  std::stringstream ss;
  ss
  << fEdep<< ","
  << fPos.x()<< ","
  << fPos.y()<< ","
  << fPos.z()<< ","
  << fKineticEnergy << ","
  << fMomentum.x() << ","
  << fMomentum.y() << ","
  << fMomentum.z() << ","
  << fMomentumDirection.x() << ","
  << fMomentumDirection.y() << ","
  << fMomentumDirection.z() << ","
  << fVelocity << ","
  << fTotalEnergy;
  return ss.str(); 
}