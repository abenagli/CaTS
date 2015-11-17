#include "TrackInformation.hh"

using namespace CLHEP;



G4Allocator<TrackInformation> aTrackInformationAllocator;

TrackInformation::TrackInformation()
{
  particleDefinition = 0;
  particleName = "";
  particleTrackID = 0;
  particlePosition = G4ThreeVector(0.,0.,0.);
  particleMomentum = G4ThreeVector(0.,0.,0.);
  particleEnergy = 0.;
  particleProdTime = 0.;
  particleTime = 0.;
  particleIsEM = false;
  particleIsNeutron = false;
  
  parentDefinition = 0;
  parentName = "";
  parentTrackID = 0;
  parentPosition = G4ThreeVector(0.,0.,0.);
  parentMomentum = G4ThreeVector(0.,0.,0.);
  parentEnergy = 0.;
  parentTime = 0.;
  parentIsEM = false;
  parentIsNeutron = false;
}



TrackInformation::TrackInformation(const G4Track* aTrack)
{
  particleDefinition = aTrack->GetDefinition();
  particleName = aTrack->GetDefinition()->GetParticleName();
  particleTrackID = aTrack->GetTrackID();
  particlePosition = aTrack->GetPosition();
  particleMomentum = aTrack->GetMomentum();
  particleEnergy = aTrack->GetTotalEnergy();
  particleTime = aTrack->GetGlobalTime();
  particleIsEM = false;
  particleIsNeutron = false;
  
  parentDefinition = aTrack->GetDefinition();
  parentName = aTrack->GetDefinition()->GetParticleName();
  parentTrackID = aTrack->GetTrackID();
  parentPosition = aTrack->GetPosition();
  parentMomentum = aTrack->GetMomentum();
  parentEnergy = aTrack->GetTotalEnergy();
  parentTime = aTrack->GetGlobalTime();
  parentIsEM = false;
  parentIsNeutron = false;
}



TrackInformation::TrackInformation(const TrackInformation* aTrackInfo)
{
  particleDefinition = aTrackInfo->particleDefinition;
  particleName = aTrackInfo->particleName;
  particleTrackID = aTrackInfo->particleTrackID;
  particlePosition = aTrackInfo->particlePosition;
  particleMomentum = aTrackInfo->particleMomentum;
  particleEnergy = aTrackInfo->particleEnergy;
  particleProdTime = aTrackInfo->particleProdTime;
  particleTime = aTrackInfo->particleTime;
  particleIsEM = aTrackInfo->particleIsEM;
  particleIsNeutron = aTrackInfo->particleIsNeutron;
  
  parentDefinition = aTrackInfo->parentDefinition;
  parentName = aTrackInfo->parentName;
  parentTrackID = aTrackInfo->parentTrackID;
  parentPosition = aTrackInfo->parentPosition;
  parentMomentum = aTrackInfo->parentMomentum;
  parentEnergy = aTrackInfo->parentEnergy;
  parentTime = aTrackInfo->parentTime;
  parentIsEM = aTrackInfo->parentIsEM;
  parentIsNeutron = aTrackInfo->parentIsNeutron;
}



void TrackInformation::SetParticleInformation(const TrackInformation* aTrackInfo)
{
  particleDefinition = aTrackInfo->particleDefinition;
  particleName = aTrackInfo->particleName;
  particleTrackID = aTrackInfo->particleTrackID;
  particlePosition = aTrackInfo->particlePosition;
  particleMomentum = aTrackInfo->particleMomentum;
  particleEnergy = aTrackInfo->particleEnergy;
  particleTime = aTrackInfo->particleTime;
  if( aTrackInfo->GetParentIsEM() ) particleIsEM = true;
  if( aTrackInfo->GetParentIsNeutron() ) particleIsNeutron = true;
}



void TrackInformation::SetParentInformation(const TrackInformation* aTrackInfo)
{
  parentDefinition = aTrackInfo->particleDefinition;
  parentName = aTrackInfo->particleName;
  parentTrackID = aTrackInfo->particleTrackID;
  parentPosition = aTrackInfo->particlePosition;
  parentMomentum = aTrackInfo->particleMomentum;
  parentEnergy = aTrackInfo->particleEnergy;
  parentTime = aTrackInfo->particleTime;
  parentIsEM = aTrackInfo->particleIsEM;
  parentIsNeutron = aTrackInfo->particleIsNeutron;
}



TrackInformation::~TrackInformation()
{}



void TrackInformation::Print() const
{
  G4cout << ">>>>>> TrackInformation::Print()::track ID " << particleTrackID 
         << " (" << particleName << ")"
         << " at " << particlePosition
         << " with " << particleEnergy/GeV
         << " GeV" << G4endl;
  
  G4cout << ">>>>>> TrackInformation::Print()::parent track ID " << parentTrackID 
         << " (" << parentName << ")"
         << " at " << parentPosition
         << " with " << parentEnergy/GeV
         << " GeV" << G4endl;
}
