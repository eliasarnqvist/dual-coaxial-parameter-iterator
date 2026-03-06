#include "EPhysicsList.hh"

EPhysicsList::EPhysicsList()
{
    // Register the physics lists with G4VModularPhysicsList
    // EM physics
    RegisterPhysics(new G4EmStandardPhysics_option4());
    // Decay physics
    RegisterPhysics(new G4DecayPhysics());
    // Radioactive decay
    RegisterPhysics(new G4RadioactiveDecayPhysics());
    // Hadronic physics
    RegisterPhysics(new G4HadronPhysicsFTFP_BERT());

    // define flags for nuclear gamma de-excitation model
    G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
    deex->SetCorrelatedGamma(true);
}

EPhysicsList::~EPhysicsList()
{
    
}
