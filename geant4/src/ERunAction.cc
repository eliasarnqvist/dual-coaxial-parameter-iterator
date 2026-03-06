#include "ERunAction.hh"

ERunAction::ERunAction()
{
    fMessenger = new G4GenericMessenger(this, "/E_file_settings/", "Settings for the filename");
	fMessenger->DeclareProperty("fileName", fileName, "Name of the output file");
    fileName = "output_temporary.root";

    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    analysisManager->CreateH2("Edet_ab", "Detected energy", 1024, 0, 5. * MeV, 1024, 0, 5. * MeV);
    analysisManager->CreateH1("Edet_a", "Detected energy", 1024, 0, 5. * MeV);
    analysisManager->CreateH1("Edet_b", "Detected energy", 1024, 0, 5. * MeV);

    analysisManager->CreateNtuple("Detectors", "Detectors");
    // analysisManager->CreateNtupleIColumn("iEvent");
    analysisManager->CreateNtupleDColumn("energy_a");
    analysisManager->CreateNtupleDColumn("energy_b");
    analysisManager->FinishNtuple(0);
}

ERunAction::~ERunAction()
{
    
}

void ERunAction::BeginOfRunAction(const G4Run *run)
{
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    // Output root file
    analysisManager->OpenFile(fileName);
}

void ERunAction::EndOfRunAction(const G4Run *run)
{
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    // Write data to our output root file and then close the file
    analysisManager->Write();
    analysisManager->CloseFile();

    G4int runID = run->GetRunID();

    G4cout << "Finishing run " << runID << G4endl;
}
