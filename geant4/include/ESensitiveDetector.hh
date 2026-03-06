#ifndef ESENSITIVEDETECTOR_HH
#define ESENSITIVEDETECTOR_HH

#include "G4VSensitiveDetector.hh"

#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

class ESensitiveDetector : public G4VSensitiveDetector
{
public:
    ESensitiveDetector(G4String);
    ~ESensitiveDetector();

private:
    G4double fTotalEnergyDeposited_a;
    G4double fTotalEnergyDeposited_b;
    G4double fLastHitTime_ab;

    G4double detectorTimeResolution = 100.0 * ns;
    
    virtual void Initialize(G4HCofThisEvent *) override;
    virtual void EndOfEvent(G4HCofThisEvent *) override;
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory *) override;

    void FlushWindow();
};

#endif