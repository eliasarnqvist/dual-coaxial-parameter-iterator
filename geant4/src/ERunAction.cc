#include "ERunAction.hh"

ERunAction::ERunAction()
{
    fMessenger = new G4GenericMessenger(this, "/E_file_settings/", "Settings for the filename");
	fMessenger->DeclareProperty("fileName", fileName, "Name of the output file");
    fileName = "output_temporary.root";

    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    analysisManager->CreateH2("E_ab_histo", "Detected energy 2D coincidence histogram", 1024, 0, 5. * MeV, 1024, 0, 5. * MeV);

    analysisManager->CreateH1("E_a_histo", "Detected energy 1D histogram", 1024, 0, 5. * MeV);

    analysisManager->CreateH1("E_b_histo", "Detected energy 1D histogram", 1024, 0, 5. * MeV);

    analysisManager->CreateNtuple("E_ab_list", "Detected energy coincidence list");
    analysisManager->CreateNtupleDColumn("energy_a");
    analysisManager->CreateNtupleDColumn("energy_b");
    analysisManager->FinishNtuple(0);

    analysisManager->CreateNtuple("E_a_list", "Detected energy list");
    analysisManager->CreateNtupleDColumn("energy_a");
    analysisManager->FinishNtuple(1);

    analysisManager->CreateNtuple("E_b_list", "Detected energy list");
    analysisManager->CreateNtupleDColumn("energy_b");
    analysisManager->FinishNtuple(2);
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
