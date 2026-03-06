#ifndef EPRIMARYGENERATOR_HH
#define EPRIMARYGENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ChargedGeantino.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4GenericMessenger.hh"
#include "G4Geantino.hh"

#include "G4GenericMessenger.hh"

class EPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
    EPrimaryGenerator();
    ~EPrimaryGenerator();

    virtual void GeneratePrimaries(G4Event *);

private:
    G4ParticleGun *fParticleGun;

    G4GenericMessenger *fMessengerSource;
    G4double sourceRadiusR;
    G4bool selectBackground;
    G4bool selectFilterSource;
    
    std::vector<G4double> xx, yy;
};

#endif