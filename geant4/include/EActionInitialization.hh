#ifndef EACTIONINITIALIZATION_HH
#define EACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

#include "EPrimaryGenerator.hh"
#include "ERunAction.hh"

class EActionInitialization : public G4VUserActionInitialization
{
public:
    EActionInitialization();
    ~EActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

#endif