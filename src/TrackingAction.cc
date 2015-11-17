#include "TrackingAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4TrackingManager.hh"

using namespace CLHEP;



TrackingAction::TrackingAction()
{}



TrackingAction::~TrackingAction()
{}



void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  //---------------------
  // tracking information
  
  if( aTrack->GetUserInformation() == 0 )
  {
    TrackInformation* aTrackInfo = new TrackInformation(aTrack);
    G4Track* theTrack = (G4Track*)aTrack;
    theTrack->SetUserInformation( aTrackInfo );
  }
  else
  {
    TrackInformation* oldTrackInfo = (TrackInformation*)(aTrack->GetUserInformation());
    
    TrackInformation* aTrackInfo = new TrackInformation(aTrack);
    oldTrackInfo -> SetParticleInformation( aTrackInfo );
    G4Track* theTrack = (G4Track*)aTrack;
    theTrack -> SetUserInformation( oldTrackInfo );
    
    delete aTrackInfo;
  }
  
  
  TrackInformation* aTrackInfo = (TrackInformation*)( aTrack->GetUserInformation() );
  
  if( aTrackInfo->GetParticleName() == "neutron" )
  {
    aTrackInfo -> SetParticleIsNeutron();
  }
}



void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  TrackInformation* aTrackInfo = (TrackInformation*)( aTrack->GetUserInformation() );
  
  //G4cout << "aTrackInfo::GetParticleName: " << aTrackInfo->GetParticleName() << G4endl;
  
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if( secondaries )
  {
    for(unsigned int i = 0; i < secondaries->size(); ++i)
    {
      G4Track* secTrack = (*secondaries)[i];
      
      TrackInformation* newTrackInfo = new TrackInformation((*secondaries)[i]);
      newTrackInfo -> SetParentInformation( aTrackInfo );
      newTrackInfo -> SetParticleProdTimeInformation( secTrack->GetGlobalTime()/picosecond );
      
      // if(newTrackInfo->GetParticleName() == "mu-" || newTrackInfo->GetParticleName() == "mu+")
      // {
      //   G4cout << newTrackInfo->GetParticleName() << " FOUND from " << aTrackInfo->GetParticleName() << " at " << newTrackInfo->GetParticleProdTime()/1000 << G4endl;
      // }
      
      if( (aTrackInfo->GetParticleIsEM()) || 
          (newTrackInfo->GetParticleName() == "gamma" && aTrackInfo->GetParticleName() == "pi0") ||
          (newTrackInfo->GetParticleName() == "e-"    && aTrackInfo->GetParticleName() == "pi0") ||
          (newTrackInfo->GetParticleName() == "e+"    && aTrackInfo->GetParticleName() == "pi0") ||
          (newTrackInfo->GetParticleName() == "gamma" && aTrackInfo->GetParticleName() == "eta") ||
          (newTrackInfo->GetParticleName() == "e-"    && aTrackInfo->GetParticleName() == "eta") ||
          (newTrackInfo->GetParticleName() == "e+"    && aTrackInfo->GetParticleName() == "eta") )
      {
        newTrackInfo -> SetParentIsEM();
        newTrackInfo -> SetParticleIsEM();
      }
      
      if( (aTrackInfo->GetParticleIsNeutron()) || 
          (aTrackInfo->GetParticleName() == "neutron") )
      {
        newTrackInfo -> SetParentIsNeutron();
        newTrackInfo -> SetParticleIsNeutron();
      }
      
      secTrack -> SetUserInformation( newTrackInfo );
    }
  }
}
