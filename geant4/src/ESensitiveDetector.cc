#include "ESensitiveDetector.hh"

ESensitiveDetector::ESensitiveDetector(G4String name) : G4VSensitiveDetector(name)
{

}

ESensitiveDetector::~ESensitiveDetector()
{
    
}

void ESensitiveDetector::Initialize(G4HCofThisEvent *)
{
    fTotalEnergyDeposited_a = 0.0;
    fTotalEnergyDeposited_b = 0.0;
    fLastHitTime_ab = -1.0;
}

void ESensitiveDetector::FlushWindow()
{
    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
    if (fTotalEnergyDeposited_a > 0) {
        analysisManager->FillH1(0, fTotalEnergyDeposited_a);

        analysisManager->FillNtupleDColumn(1, 0, fTotalEnergyDeposited_a);
        analysisManager->AddNtupleRow(1);
    }
    if (fTotalEnergyDeposited_b > 0) {
        analysisManager->FillH1(1, fTotalEnergyDeposited_b);

        analysisManager->FillNtupleDColumn(2, 0, fTotalEnergyDeposited_b);
        analysisManager->AddNtupleRow(2);
    }
    if (fTotalEnergyDeposited_a > 0 && fTotalEnergyDeposited_b > 0)
    {
        analysisManager->FillH2(0, fTotalEnergyDeposited_a, fTotalEnergyDeposited_b);

        analysisManager->FillNtupleDColumn(0, 0, fTotalEnergyDeposited_a);
        analysisManager->FillNtupleDColumn(0, 1, fTotalEnergyDeposited_b);
        analysisManager->AddNtupleRow(0);
    }

    // reset the deposited energy
    fTotalEnergyDeposited_a = 0;
    fTotalEnergyDeposited_b = 0;
}

void ESensitiveDetector::EndOfEvent(G4HCofThisEvent *)
{
    FlushWindow();
}

G4bool ESensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    G4double fEnergyDeposited = aStep->GetTotalEnergyDeposit();

    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int copyNo = touchable->GetCopyNumber(3);
    // the number in GetCopyNumber is the order of the mother logic volume (0=this one, 1=1st mother, ...)
    // copyNo for detectors used here: a=0, b=1
    // G4cout << "Copy number: " << copyNo << G4endl;

    if (fEnergyDeposited > 0) {
        G4double fHitTime = aStep->GetPreStepPoint()->GetGlobalTime();

        if (copyNo == 0) {
            // detector a was hit
            if (fLastHitTime_ab < 0.0) {
                fLastHitTime_ab = fHitTime;
            } else if ((fHitTime - fLastHitTime_ab) > detectorTimeResolution) {
                FlushWindow();
                fLastHitTime_ab = fHitTime;
            }
            // accumulate deposited energy
            fTotalEnergyDeposited_a += fEnergyDeposited;
        } else if (copyNo == 1) {
            // detector b was hit
            if (fLastHitTime_ab < 0.0) {
                fLastHitTime_ab = fHitTime;
            } else if ((fHitTime - fLastHitTime_ab) > detectorTimeResolution) {
                FlushWindow();
                fLastHitTime_ab = fHitTime;
            }
            // accumulate deposited energy
            fTotalEnergyDeposited_b += fEnergyDeposited;
        }

        return true;
    } else {
        return false;
    }
}


