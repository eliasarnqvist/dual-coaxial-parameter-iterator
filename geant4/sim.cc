// Include standard input output library so we can print to the console
#include <iostream>

// Include the Geant4 run manager and multi-threaded run manager
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
// Include user interface and visualization libraries
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

// Include our own physics list, detector construction and action initialization
#include "EPhysicsList.hh"
#include "EDetectorConstruction.hh"
#include "EActionInitialization.hh"

int main(int argc, char** argv)
{
    // Prepare user interface pointer
    G4UIExecutive *ui = 0;

    // Check if the Geant4 installation supports multi-threading and use it if possible
    #ifdef G4MULTITHREADED
        G4MTRunManager *runManager = new G4MTRunManager;
    #else
        G4RunManager *runManager = new G4RunManager;
    #endif

    // Three things that always need to be set up with run manager
    // Physics list
    runManager->SetUserInitialization(new EPhysicsList());
    // Detector construction
    runManager->SetUserInitialization(new EDetectorConstruction());
    // Action initialization
    runManager->SetUserInitialization(new EActionInitialization());

    // Here, argc is the argument count, so if we run ./sim then arg=1
    // If there are no command line arguments make ui pointer
    if (argc == 1)
    {
        ui = new G4UIExecutive(argc, argv);
    }

    // Initialize the visualization manager
    G4VisManager *visManager = new G4VisExecutive();
    visManager->Initialize();
    // Get the UI manager
    G4UImanager *UImanager = G4UImanager::GetUIpointer();

    // If ui executive is created, execute the visualization macro and start the session
    // This would open the separate window for visualization
    // Otherwise, we run Geant4 in batch mode from the command line with macro files
    if(ui)
    {
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
    }
    else
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }

    // Clean up and delete the visualization manager and run manager, avoid some warnings by using these
    delete visManager;
    delete runManager;

    return 0;
}