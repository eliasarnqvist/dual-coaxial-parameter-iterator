#ifndef ERUNACTION_HH
#define ERUNACTION_HH

#include "G4UserRunAction.hh"

#include "G4Run.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4GenericMessenger.hh"

class ERunAction : public G4UserRunAction
{
public:
    ERunAction();
    ~ERunAction();

    virtual void BeginOfRunAction(const G4Run *);
    virtual void EndOfRunAction(const G4Run *);

private:
    G4GenericMessenger *fMessenger;
    G4String fileName;
};

#endif