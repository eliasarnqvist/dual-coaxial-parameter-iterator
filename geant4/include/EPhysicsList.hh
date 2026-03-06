#ifndef EPHYSICSLIST_HH
#define EPHYSICSLIST_HH

// Use G4VModularPhysicsList to make physics lists
// It is one of three ways to make physics lists
#include "G4VModularPhysicsList.hh"

// Physics lists for describing particles and particle interactions
// Electromagnetic interactions (gammas, electrons, positrons, ...)
#include "G4EmStandardPhysics_option4.hh" // most accurate and recommended for ~MeV and below
// Defines standard list of particles and decay processes
#include "G4DecayPhysics.hh"
// Defines radioactive decay of isotopes
#include "G4RadioactiveDecayPhysics.hh"
// For hadronic interactions
#include "G4HadronPhysicsFTFP_BERT.hh"

// For handeling nuclear de-excitation 
#include "G4NuclideTable.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"

#include "G4SystemOfUnits.hh"

class EPhysicsList : public G4VModularPhysicsList
{
public:
    EPhysicsList();
    ~EPhysicsList();
};

#endif