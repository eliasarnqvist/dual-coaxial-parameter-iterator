#include "EDetectorConstruction.hh"

EDetectorConstruction::EDetectorConstruction()
{
    fMessengerDetector = new G4GenericMessenger(this, "/E_detector/", "Settings for the detector");
    fMessengerDetector->DeclareProperty("detectorDiameter", detectorDiameterR, "Diameter of the detectors (mm)");
    fMessengerDetector->DeclareProperty("detectorLength", detectorLengthR, "Length of the detectors (mm)");

    fMessengerSource = new G4GenericMessenger(this, "/E_source/", "Settings for the source");
    fMessengerSource->DeclareProperty("selectFilterSource", selectFilterSource, "If true, make a filter source");

    detectorDiameterR = 60.;
    detectorLengthR = 60.;

    selectFilterSource = true;
    selectNTypeInsteadOfPType = false;
}

EDetectorConstruction::~EDetectorConstruction()
{
    
}

G4VPhysicalVolume *EDetectorConstruction::Construct()
{
    G4double detectorDiameter = detectorDiameterR * mm;
    G4double detectorLength = detectorLengthR * mm;

    G4bool checkOverlaps = true;

    // These materials will be used
    // Define materials
    G4NistManager *nist = G4NistManager::Instance();
    MatWorld = nist->FindOrBuildMaterial("G4_AIR");
    MatGe = nist->FindOrBuildMaterial("G4_Ge");
    MatAl = nist->FindOrBuildMaterial("G4_Al");
    MatCu = nist->FindOrBuildMaterial("G4_Cu");
    MatVac = nist->FindOrBuildMaterial("G4_Galactic");
    // Carbon fiber reinforced polymer (CFRP) material
    MatCFRP = new G4Material("CFRP", 1.4 * g/cm3, 2);
    MatCFRP->AddElement(nist->FindOrBuildElement("C"), 0.85);
    MatCFRP->AddElement(nist->FindOrBuildElement("O"), 0.15);

    // World
    G4double xWorld = 1.0 * m;
    G4double yWorld = 1.0 * m;
    G4double zWorld = 1.0 * m;
    G4Box *solidWorld = new G4Box("solidWorld", 0.5 * xWorld, 0.5 * yWorld, 0.5 * zWorld);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, MatWorld, "logicWorld");
    G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, checkOverlaps);

    // G4double testRadius = 104.0 * mm;
    // G4Orb *solidSURE = new G4Orb("solidSURE", testRadius);
    // G4LogicalVolume *logicSURE = new G4LogicalVolume(solidSURE, MatWorld, "logicSURE");
    // G4VPhysicalVolume *physSURE = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicSURE, "physSURE", logicWorld, false, 0, checkOverlaps);

    if (selectFilterSource)
    {
        // These dimensions need to be changed in EPrimaryGenerator as well!
        G4double filterDiameter = 31.0 * 2 * mm;
        G4double filterHeight = 13.8 * mm;
        G4double beakerThickness = 1.0 * mm;
        ConstructFilterSource(logicWorld, checkOverlaps, filterDiameter, filterHeight, beakerThickness);

        G4double detectorFrontPosition = 0.5 * filterHeight + beakerThickness;
        logicDetector_a = ConstructHPGe(logicWorld, detectorFrontPosition, 0. * deg, detectorDiameter, detectorLength, 0, checkOverlaps);
        logicDetector_b = ConstructHPGe(logicWorld, detectorFrontPosition, 180. * deg, detectorDiameter, detectorLength, 1, checkOverlaps);
    }
    else
    {
        logicDetector_a = ConstructHPGe(logicWorld, 0. * mm, 0. * deg, detectorDiameter, detectorLength, 0, checkOverlaps);
        logicDetector_b = ConstructHPGe(logicWorld, 0. * mm, 180. * deg, detectorDiameter, detectorLength, 1, checkOverlaps);
        // logicDetector_b = ConstructHPGe(logicWorld, detectorDistance, detectorAngle, detectorDiameter, detectorLength, capWallThickness, 1, checkOverlaps);
    }

    // Set smaller range cuts for regions with small geometries, as the default is 1 mm
    G4ProductionCuts *cutsThinDeadLayer = new G4ProductionCuts();
    cutsThinDeadLayer->SetProductionCut(0.2 * um); // note the unit
    regionThinDeadLayer->SetProductionCuts(cutsThinDeadLayer);

    G4ProductionCuts *cutsThickDeadLayer = new G4ProductionCuts();
    cutsThickDeadLayer->SetProductionCut(0.5 * mm);
    regionThickDeadLayer->SetProductionCuts(cutsThickDeadLayer);

    G4ProductionCuts *cutsActiveRegion = new G4ProductionCuts();
    cutsActiveRegion->SetProductionCut(1 * mm);
    regionActiveRegion->SetProductionCuts(cutsActiveRegion);

    G4ProductionCuts *cutsCapWindow = new G4ProductionCuts();
    cutsCapWindow->SetProductionCut(0.5 * mm);
    regionCapWindow->SetProductionCuts(cutsCapWindow);

    return physWorld;
}

void EDetectorConstruction::ConstructSDandField()
{
    ESensitiveDetector *sensDet = new ESensitiveDetector("ESensitiveDetector");
    logicDetector_a->SetSensitiveDetector(sensDet);
    logicDetector_b->SetSensitiveDetector(sensDet);
    G4SDManager::GetSDMpointer()->AddNewDetector(sensDet);
}

G4LogicalVolume *EDetectorConstruction::ConstructHPGe(G4LogicalVolume* logicWorld, G4double detectorDistance, G4double detectorAngle, G4double detectorDiameter, G4double detectorLength, G4int copyNo, G4bool checkOverlaps)
{
    G4double capWallThickness = 1.5 * mm;
    G4double vacFrontSpace = 4.0 * mm;
    G4double detectorOuterDeadLayer = 0.3 * um; // note the unit
    G4double detectorInnerDeadLayer = 0.7 * mm;

    G4double holderRingThickness = 2.0 * mm;
    G4double holderThickness = 1.0 * mm;
    G4double holderBottomSpace = 5.0 * mm;
    G4double holderBottomThickness = 5.0 * mm;
    G4double holderLength = detectorLength + holderBottomThickness + holderBottomSpace;
    G4double holderOuterDiameter = detectorDiameter + 2 * holderThickness;
    G4double holderRingOuterDiameter = detectorDiameter + 2 * holderRingThickness;
    G4double holderRingLength = 10.0 * mm;
    G4double holderRingSpacing = 10.0 * mm;

    G4double detectorBevelRadius = 4.0 * mm;
    G4double detectorHoleDiameter = 10.0 * mm;
    G4double detectorHoleDepth = 0.75 * detectorLength;
    G4double detectorInternalBevelRadius = detectorBevelRadius - detectorOuterDeadLayer;
    G4double detectorInternalHoleDiameter = detectorHoleDiameter + 2 * detectorInnerDeadLayer;
    G4double detectorInternalHoleDepth = detectorHoleDepth + detectorInnerDeadLayer;
    G4double detectorInternalDiameter = detectorDiameter - 2 * detectorOuterDeadLayer;
    G4double detectorInternalLength = detectorLength - detectorOuterDeadLayer;

    G4double coldfingerDiameter = 15.0 * mm;
    G4double coldfingerLength = 10.0 * mm;

    G4double capBottomThickness = 5.0 * mm;
    G4double capOuterLength = capWallThickness + vacFrontSpace + holderLength + coldfingerLength + capBottomThickness;
    G4double capSideSpace = 3.0 * mm;
    G4double capOuterDiameter = detectorDiameter + 2 * holderRingThickness + 2 * capSideSpace + 2 * capWallThickness;
    
    G4double vacDiameter = capOuterDiameter - 2 * capWallThickness;
    G4double vacLength = capOuterLength - capWallThickness - capBottomThickness;

    G4double windowThickness = 0.6 * mm;
    G4double windowDiameter = 0.8 * capOuterDiameter;
    G4double windowVacThickness = capWallThickness - windowThickness;
    
    G4double sliceAngle = 180. * deg;

    // Cap
    G4Tubs *solidCap = new G4Tubs("solidCap", 0. * mm, 0.5 * capOuterDiameter, 0.5 * capOuterLength, 0. * deg, sliceAngle);
    G4LogicalVolume *logicCap = new G4LogicalVolume(solidCap, MatAl, "logicCap");
    G4Rotate3D myRoration(detectorAngle, G4ThreeVector(0., 1., 0.));
    G4Translate3D myTranslation(G4ThreeVector(0., 0., 0.5 * (capOuterLength) + detectorDistance));
    G4Transform3D myTransform = (myRoration) * (myTranslation);
    G4VPhysicalVolume *physCap = new G4PVPlacement(myTransform, logicCap, "physCap", logicWorld, false, copyNo, checkOverlaps);
    
    // Window
    G4Tubs *solidWindow = new G4Tubs("solidWindow", 0. * mm, 0.5 * windowDiameter, 0.5 * windowThickness, 0. * deg, sliceAngle);
    G4LogicalVolume *logicWindow = new G4LogicalVolume(solidWindow, MatCFRP, "logicWindow");
    G4VPhysicalVolume *physWindow = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5 * (windowThickness - capOuterLength)), logicWindow, "physWindow", logicCap, false, 0, checkOverlaps);

    G4Tubs *solidWindowVacuum = new G4Tubs("solidWindowVacuum", 0. * mm, 0.5 * windowDiameter, 0.5 * windowVacThickness, 0. * deg, sliceAngle);
    G4LogicalVolume *logicWindowVacuum = new G4LogicalVolume(solidWindowVacuum, MatVac, "logicWindowVacuum");
    G4VPhysicalVolume *physWindowVacuum = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5 * (windowVacThickness - capOuterLength) + windowThickness), logicWindowVacuum, "physWindowVacuum", logicCap, false, 0, checkOverlaps);

    // Vacuum inside cap
    G4Tubs *solidVacuum = new G4Tubs("solidVacuum", 0. * mm, 0.5 * vacDiameter, 0.5 * vacLength, 0. * deg, sliceAngle);
    G4LogicalVolume *logicVacuum = new G4LogicalVolume(solidVacuum, MatVac, "logicVacuum");
    G4VPhysicalVolume *physVacuum = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5 * (capWallThickness - capBottomThickness)), logicVacuum, "physVacuum", logicCap, false, 0, checkOverlaps);

    // Crystal (including the dead layers)
    G4Tubs *solidCrystal1 = new G4Tubs("solidCrystal1", 0. * mm, 0.5 * detectorDiameter - detectorBevelRadius, 0.5 * detectorBevelRadius, 0. * deg, sliceAngle);
    G4Tubs *solidCrystal2 = new G4Tubs("solidCrystal2", 0. * mm, 0.5 * detectorDiameter, 0.5 * (detectorLength - detectorHoleDepth - detectorBevelRadius), 0. * deg, sliceAngle);
    G4Tubs *solidCrystal3 = new G4Tubs("solidCrystal3", 0.5 * detectorHoleDiameter, 0.5 * detectorDiameter, 0.5 * detectorHoleDepth, 0. * deg, sliceAngle);
    G4Torus *solidCrystal4 = new G4Torus("solidCrystal4", 0. * mm, detectorBevelRadius, 0.5 * detectorDiameter - detectorBevelRadius, 0. * deg, sliceAngle);
    G4MultiUnion* munionCrystal = new G4MultiUnion("munionCrystal");
    munionCrystal->AddNode(*solidCrystal1, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., vacFrontSpace + 0.5 * (detectorBevelRadius))));
    munionCrystal->AddNode(*solidCrystal2, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., vacFrontSpace + detectorBevelRadius + 0.5 * (detectorLength - detectorHoleDepth - detectorBevelRadius))));
    munionCrystal->AddNode(*solidCrystal3, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., vacFrontSpace + detectorLength - 0.5 * detectorHoleDepth)));
    munionCrystal->AddNode(*solidCrystal4, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., vacFrontSpace + detectorBevelRadius)));
    munionCrystal->Voxelize();
    G4LogicalVolume *logicCrystal = new G4LogicalVolume(munionCrystal, MatGe, "logicCrystal");
    G4VPhysicalVolume *physCrystal = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5 * (- vacLength)), logicCrystal, "physCrystal", logicVacuum, false, 0, checkOverlaps);

    // Crystal (only the sensitive part)
    G4Tubs *solidDetector1 = new G4Tubs("solidDetector1", 0. * mm, 0.5 * detectorInternalDiameter - detectorInternalBevelRadius, 0.5 * detectorInternalBevelRadius, 0. * deg, sliceAngle);
    G4Tubs *solidDetector2 = new G4Tubs("solidDetector2", 0. * mm, 0.5 * detectorInternalDiameter, 0.5 * (detectorInternalLength - detectorInternalHoleDepth - detectorInternalBevelRadius), 0. * deg, sliceAngle);
    G4Tubs *solidDetector3 = new G4Tubs("solidDetector3", 0.5 * detectorInternalHoleDiameter, 0.5 * detectorInternalDiameter, 0.5 * detectorInternalHoleDepth, 0. * deg, sliceAngle);
    G4Torus *solidDetector4 = new G4Torus("solidDetector4", 0. * mm, detectorInternalBevelRadius, 0.5 * detectorInternalDiameter - detectorInternalBevelRadius, 0. * deg, sliceAngle);
    G4MultiUnion* munionDetector = new G4MultiUnion("munionDetector");
    munionDetector->AddNode(*solidDetector1, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., detectorOuterDeadLayer + 0.5 * (detectorInternalBevelRadius))));
    munionDetector->AddNode(*solidDetector2, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., detectorOuterDeadLayer + detectorInternalBevelRadius + 0.5 * (detectorInternalLength - detectorInternalHoleDepth - detectorInternalBevelRadius))));
    munionDetector->AddNode(*solidDetector3, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., detectorOuterDeadLayer + detectorInternalLength - 0.5 * detectorInternalHoleDepth)));
    munionDetector->AddNode(*solidDetector4, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., detectorOuterDeadLayer + detectorInternalBevelRadius)));
    munionDetector->Voxelize();
    G4LogicalVolume *logicDetector = new G4LogicalVolume(munionDetector, MatGe, "logicDetector");
    G4VPhysicalVolume *physDetector = new G4PVPlacement(0, G4ThreeVector(0., 0., vacFrontSpace), logicDetector, "physDetector", logicCrystal, false, 0, checkOverlaps);

    // Crystal (only inner dead layer)
    G4Tubs *solidCrystalInner1 = new G4Tubs("solidCrystalInner1", 0. * mm, 0.5 * detectorInternalHoleDiameter, 0.5 * detectorInnerDeadLayer, 0. * deg, sliceAngle);
    G4Tubs *solidCrystalInner2 = new G4Tubs("solidCrystalInner2", 0.5 * detectorHoleDiameter, 0.5 * detectorInternalHoleDiameter, 0.5 * detectorHoleDepth, 0. * deg, sliceAngle);
    G4MultiUnion* munionCrystalInner = new G4MultiUnion("munionCrystalInner");
    munionCrystalInner->AddNode(*solidCrystalInner1, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., 0.5 * detectorInnerDeadLayer)));
    munionCrystalInner->AddNode(*solidCrystalInner2, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., detectorInnerDeadLayer + 0.5 * detectorHoleDepth)));
    munionCrystalInner->Voxelize();
    G4LogicalVolume *logicCrystalInner = new G4LogicalVolume(munionCrystalInner, MatGe, "logicCrystalInner");
    G4VPhysicalVolume *physCrystalInner = new G4PVPlacement(0, G4ThreeVector(0., 0., vacFrontSpace + detectorLength - detectorHoleDepth - detectorInnerDeadLayer), logicCrystalInner, "physCrystalInner", logicCrystal, false, 0, checkOverlaps);

    // Holder
    G4Tubs *solidHolderSides = new G4Tubs("solidHolderSides", 0.5 * detectorDiameter, 0.5 * holderOuterDiameter, 0.5 * holderLength, 0. * deg, sliceAngle);
    G4Tubs *solidHolderTopRing = new G4Tubs("solidHolderTopRing", 0.5 * detectorDiameter, 0.5 * holderRingOuterDiameter, 0.5 * holderRingLength, 0. * deg, sliceAngle);
    G4Tubs *solidHolderBottomRing = new G4Tubs("solidHolderBottomRing", 0.5 * detectorDiameter, 0.5 * holderRingOuterDiameter, 0.5 * holderRingLength, 0. * deg, sliceAngle);
    G4Tubs *solidHolderBottom = new G4Tubs("solidHolderBottom", 0. * mm, 0.5 * detectorDiameter, 0.5 * holderBottomThickness, 0. * deg, sliceAngle);
    G4MultiUnion* munionHolder = new G4MultiUnion("munionHolder");
    munionHolder->AddNode(*solidHolderSides, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., vacFrontSpace + 0.5 * holderLength)));
    munionHolder->AddNode(*solidHolderTopRing, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., vacFrontSpace + 0.5 * holderRingLength)));
    munionHolder->AddNode(*solidHolderBottomRing, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., vacFrontSpace + holderRingSpacing + 1.5 * holderRingLength)));
    munionHolder->AddNode(*solidHolderBottom, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., vacFrontSpace + holderLength - 0.5 * holderBottomThickness)));
    munionHolder->Voxelize();
    G4LogicalVolume *logicHolder = new G4LogicalVolume(munionHolder, MatAl, "logicHolder");
    G4VPhysicalVolume *physHolder = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5 * (- vacLength)), logicHolder, "physHolder", logicVacuum, false, 0, checkOverlaps);

    // Cold finger
    G4Tubs *solidColdfinger = new G4Tubs("solidColdfinger", 0. * mm, 0.5 * coldfingerDiameter, 0.5 * coldfingerLength, 0. * deg, sliceAngle);
    G4LogicalVolume *logicColdfinger = new G4LogicalVolume(solidColdfinger, MatCu, "logicColdfinger");
    G4VPhysicalVolume *physColdfinger = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5 * (vacLength - coldfingerLength)), logicColdfinger, "physColdfinger", logicVacuum, false, 0, checkOverlaps);

    // Set appropriate range cuts by assigning logic volumes to regions
    logicDetector->SetRegion(regionActiveRegion);
    regionActiveRegion->AddRootLogicalVolume(logicDetector);

    logicCrystal->SetRegion(regionThinDeadLayer);
    regionThinDeadLayer->AddRootLogicalVolume(logicCrystal);

    // logicWindow->SetRegion(regionThickDeadLayer);
    // regionThickDeadLayer->AddRootLogicalVolume(logicWindow);

    logicWindow->SetRegion(regionCapWindow);
    regionCapWindow->AddRootLogicalVolume(logicWindow);


    // Show pretty colors in the visualization
    // G4VisAttributes *capVisAtt = new G4VisAttributes(G4Color(1.0, 0.0, 0.0, 0.5));
    // capVisAtt->SetForceSolid(true);
    // logicCap->SetVisAttributes(capVisAtt);
    G4VisAttributes *detVisAtt = new G4VisAttributes(G4Color(1.0, 1.0, 0.0, 0.5));
    detVisAtt->SetForceSolid(true);
    logicDetector->SetVisAttributes(detVisAtt);
    // G4VisAttributes *holVisAtt = new G4VisAttributes(G4Color(0.0, 1.0, 1.0, 0.5));
    // holVisAtt->SetForceSolid(true);
    // logicHolder->SetVisAttributes(holVisAtt);

    return logicDetector;
}

void EDetectorConstruction::ConstructFilterSource(G4LogicalVolume* logicWorld, G4bool checkOverlaps, G4double filterDiameter, G4double filterHeight, G4double beakerThickness)
{
    G4NistManager *nist = G4NistManager::Instance();
    G4Material* MatBeaker = nist->FindOrBuildMaterial("G4_POLYSTYRENE");

    // For the filter material
    G4Element *eSi = new G4Element("Silicon", "Si", 14., 28.085*g/mole);
    G4Element *eNa = new G4Element("Sodium", "Na", 11., 22.990*g/mole);
    G4Element *eAl = new G4Element("Aluminium", "Al", 13., 26.982*g/mole);
    G4Element *eBa = new G4Element("Barium", "Ba", 56., 137.33*g/mole);
    G4Element *eCa = new G4Element("Calcium", "Ca", 20., 40.078*g/mole);
    G4Element *eB = new G4Element("Boron", "B", 5., 10.81*g/mole);
    G4Element *eZn = new G4Element("Zink", "Zn", 30., 65.38*g/mole);
    G4Element *eK = new G4Element("Potasium", "K", 19., 39.098*g/mole);
    G4Element *eO = new G4Element("Oxygen", "O", 8., 15.999*g/mole);
    G4Element *eC = new G4Element("Carbon", "C", 6., 12.011*g/mole);
    // Density of pressed filter (total mass = 40 g)
    G4Material *MatFilter = new G4Material("FOI_filter", 0.96*g/cm3, 10);
    MatFilter->AddElement(eSi, 0.282524);
    MatFilter->AddElement(eNa, 0.060861);
    MatFilter->AddElement(eAl, 0.032884);
    MatFilter->AddElement(eBa, 0.034297);
    MatFilter->AddElement(eCa, 0.024823);
    MatFilter->AddElement(eB , 0.030021);
    MatFilter->AddElement(eZn, 0.024166);
    MatFilter->AddElement(eK , 0.023044);
    MatFilter->AddElement(eO , 0.452465);
    MatFilter->AddElement(eC , 0.034916);

    // G4double filterDiameter = 31.0 * 2 * mm;
    // G4double filterHeight = 13.8 * mm;

    // G4double beakerThickness = 1.0 * mm;
    G4double beakerDiameter = 2 * beakerThickness + filterDiameter;
    G4double beakerHeight = 2 * beakerThickness + filterHeight;

    G4Tubs *solidBeaker = new G4Tubs("solidBeaker", 0., 0.5 * beakerDiameter, 0.5 * beakerHeight, 0., 360.*degree);
    G4LogicalVolume *logicBeaker = new G4LogicalVolume(solidBeaker, MatBeaker, "logicBeaker");
    G4VPhysicalVolume *physBeaker = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicBeaker, "physBeaker", logicWorld, false, 0, checkOverlaps);

    G4Tubs *solidFilter = new G4Tubs("solidFilter", 0., 0.5 * filterDiameter, 0.5 * filterHeight, 0., 360.*degree);
    G4LogicalVolume *logicFilter = new G4LogicalVolume(solidFilter, MatFilter, "logicFilter");
    G4VPhysicalVolume *physFilter = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicFilter, "physFilter", logicBeaker, false, 0, checkOverlaps);

    // Show pretty colors in the visualization
    G4VisAttributes *filterVisAtt = new G4VisAttributes(G4Color(1.0, 0.0, 0.0, 0.5));
    filterVisAtt->SetForceSolid(true);
    logicFilter->SetVisAttributes(filterVisAtt);
}


