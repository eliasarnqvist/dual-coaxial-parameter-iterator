#ifndef EDETECTORCONSTRUCTION_HH
#define EDETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4MultiUnion.hh"

#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4SDManager.hh"

#include "ESensitiveDetector.hh"

#include "G4GenericMessenger.hh"

#include "G4Region.hh"
#include "G4ProductionCuts.hh"

class EDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    EDetectorConstruction();
    virtual ~EDetectorConstruction();

    virtual G4VPhysicalVolume *Construct();

private:
    G4LogicalVolume *logicDetector_a;
    G4LogicalVolume *logicDetector_b;
    ESensitiveDetector *sensDet = nullptr;

    G4Material *MatWorld, *MatGe, *MatAl, *MatCu, *MatVac, *MatCFRP;

    G4GenericMessenger *fMessengerDetector, *fMessengerSource;
    G4double detectorDiameterR, detectorLengthR, sourceDistanceR;
    G4bool selectFilterSource, selectNTypeInsteadOfPType;

    virtual void ConstructSDandField();
    virtual G4LogicalVolume *ConstructHPGe(G4LogicalVolume* logicWorld, G4double detectorDistance, G4double detectorAngle, G4double detectorDiameter, G4double detectorLength, G4bool selectNTypeInsteadOfPType, G4int copyNo, G4bool checkOverlaps);
    virtual void ConstructFilterSource(G4LogicalVolume* logicWorld, G4bool checkOverlaps, G4double filterDiameter, G4double filterHeight, G4double beakerThickness);

    G4Region *regionThinDeadLayer = new G4Region("regionThinDeadLayer");
    G4Region *regionThickDeadLayer = new G4Region("regionThickDeadLayer");
    G4Region *regionActiveRegion = new G4Region("regionActiveRegion");
    G4Region *regionCapWindow = new G4Region("regionCapWindow");
};

#endif