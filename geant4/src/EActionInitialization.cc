#include "EActionInitialization.hh"

EActionInitialization::EActionInitialization()
{

}

EActionInitialization::~EActionInitialization()
{
    
}

void EActionInitialization::BuildForMaster() const
{
    ERunAction *runAction = new ERunAction();
    SetUserAction(runAction);
}

void EActionInitialization::Build() const
{
    EPrimaryGenerator *generator = new EPrimaryGenerator();
    SetUserAction(generator);

    ERunAction *runAction = new ERunAction();
    SetUserAction(runAction);
}