#include "EPrimaryGenerator.hh"

EPrimaryGenerator::EPrimaryGenerator()
{
    fMessengerSource = new G4GenericMessenger(this, "/E_source/", "Settings for the source");
    fMessengerSource->DeclareProperty("selectBackground", selectBackground, "If true, simulate a background source");
    fMessengerSource->DeclareProperty("sourceRadius", sourceRadiusR, "Radius of the source (mm)");
    fMessengerSource->DeclareProperty("selectFilterSource", selectFilterSource, "If true, make a filter source");

    selectBackground = false;
    selectFilterSource = true;
    sourceRadiusR = 104.;

    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);

    fParticleGun->SetParticleEnergy(0 * eV);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
}

EPrimaryGenerator::~EPrimaryGenerator()
{
    delete fParticleGun;
}

void EPrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
    // If the particle has not bee assigned yet, assign it as either a gamma ray or an ion
    if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino())
    {
        if (selectBackground)
        {
            // Read gamma ray flux CDF (should work both in the build directory and outside!)
            std::ifstream datafile;
            datafile.open("../resources/flux_cdf.dat");
            if (!datafile.is_open()) {
                G4cout << "Could not open ../resources/flux_cdf.dat, trying resources/flux_cdf.dat" << G4endl;
                datafile.open("resources/flux_cdf.dat");
            }
            G4double x, y;
            if (datafile.is_open()) {
                while (datafile >> x >> y) {
                    xx.push_back(x);
                    yy.push_back(y);
                }
                datafile.close();
                G4cout << "Read " << xx.size() << " entries from ../resources/flux_cdf.dat" << G4endl;
            } else {
                G4cerr << "Error: Could not open either ../resources/flux_cdf.dat or resources/flux_cdf.dat" << G4endl;
            }
    
            G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
            G4ParticleDefinition *gammaRay = particleTable->FindParticle("gamma");

            fParticleGun->SetParticleDefinition(gammaRay);
        }
        else
        {
            G4int Z = 57, A = 140;
            G4double excitEnergy = 0. * keV;
            G4IonTable *ionTable = G4IonTable::GetIonTable();
            G4ParticleDefinition *ion = ionTable->GetIon(Z, A, excitEnergy);

            fParticleGun->SetParticleDefinition(ion);
            fParticleGun->SetParticleCharge(0. * eplus);
        }
    }

    // If we simulate a gamma ray, the energy needs to follow the CDF and the position needs to be random
    if (selectBackground)
    {
        // Set the energy of the gamma ray
        G4double CFD_sample = G4UniformRand();
        G4int max_index = xx.size() - 1;
        G4int j = 0;
        while (j < max_index && yy[j] < CFD_sample) {
            j++;
        }
        // G4cout << "CDF sample " << CFD_sample << " out of 1" << G4endl;
        // G4cout << "Gamma energy " << xx[j] << " keV" << G4endl;
        // Spread the emitted gamma-rays equally over the whole bin
        G4double gammaSpread = 4.78515625 * (G4UniformRand() - 0.5);
        G4double gammaEnergy = (xx[j] + gammaSpread) * keV;
        // G4cout << "Gamma energy " << gammaEnergy << " MeV" << G4endl;
        fParticleGun->SetParticleEnergy(gammaEnergy);

        G4double sourceRadius = sourceRadiusR * mm;

        // Random position
        G4double phi = 2  * pi * G4UniformRand();
        G4double theta = std::acos(1 - 2 * G4UniformRand());

        G4double rhat_x = std::sin(theta) * std::cos(phi);
        G4double rhat_y = std::sin(theta) * std::sin(phi);
        G4double rhat_z = std::cos(theta);

        G4double thetahat_x = std::cos(theta) * std::cos(phi);
        G4double thetahat_y = std::cos(theta) * std::sin(phi);
        G4double thetahat_z = - std::sin(theta);

        G4double phihat_x = - std::sin(phi);
        G4double phihat_y = std::cos(phi);
        G4double phihat_z = 0;

        G4double disk_phi = G4UniformRand() * 2 * pi;
        G4double disk_r = std::sqrt(G4UniformRand() * sourceRadius * sourceRadius);

        G4double x = sourceRadius * rhat_x + disk_r * (std::cos(disk_phi) * thetahat_x + std::sin(disk_phi) * phihat_x);
        G4double y = sourceRadius * rhat_y + disk_r * (std::cos(disk_phi) * thetahat_y + std::sin(disk_phi) * phihat_y);
        G4double z = sourceRadius * rhat_z + disk_r * (std::cos(disk_phi) * thetahat_z + std::sin(disk_phi) * phihat_z);

        // New random position but pointing towards detector
        fParticleGun->SetParticlePosition(G4ThreeVector(x,y,z));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-rhat_x,-rhat_y,-rhat_z));
    }
    else if (selectFilterSource)
    {
        // New random position inside filter source
        G4double filterDiameter = 31.0 * 2 * mm;
        G4double filterHeight = 13.8 * mm;
        // Randomized position (corrected for uniform sampling in cylinder)
        G4double rho = G4UniformRand() * (filterDiameter/2) * (filterDiameter/2);
        G4double phi = G4UniformRand() * 2 * pi;

        G4double x = std::cos(phi) * std::sqrt(rho);
        G4double y = std::sin(phi) * std::sqrt(rho);
        G4double z = filterHeight * (G4UniformRand() - 0.5);

        fParticleGun->SetParticlePosition(G4ThreeVector(x,y,z));
    }

    // Create vertex
    fParticleGun->GeneratePrimaryVertex(anEvent);
}


