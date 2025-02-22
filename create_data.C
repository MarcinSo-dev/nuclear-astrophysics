#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>

void create_data() {
    // Create and open the ROOT file
    TFile *file = new TFile("beam_intensity.root", "RECREATE");

    // Create directories
    TDirectory *dataDir = file->mkdir("data");
    dataDir->cd();

    // Create a TTree for constants
    TTree *constantsTree = new TTree("constants", "Experimental Constants");

    // Define and store only necessary constants
    double S_n = 6.512;   // Threshold energy (MeV)
    double N_A = 6.023e23;  // Avogadro's number (atoms/mol)
    double A = 197;   // Atomic mass (gold)
    double rho = 19.283;  // Density (g/cm^3)
    double x = 25.4e-4;   // Sample thickness (cm)
    double t_IRR = 59 * 60;  // Irradiation time (s)
    double t_delay = 19 * 60;  // Delay time (s)
    double t_measurement = 12 * 3600;  // Measurement time (s)
    double half_life = 6.1669 * 24 * 3600;  // Half-life (s)

    // Store only these in the ROOT file
    constantsTree->Branch("S_n", &S_n, "S_n/D");
    constantsTree->Branch("N_A", &N_A, "N_A/D");
    constantsTree->Branch("A", &A, "A/D");
    constantsTree->Branch("rho", &rho, "rho/D");
    constantsTree->Branch("x", &x, "x/D");
    constantsTree->Branch("t_IRR", &t_IRR, "t_IRR/D");
    constantsTree->Branch("t_delay", &t_delay, "t_delay/D");
    constantsTree->Branch("t_measurement", &t_measurement, "t_measurement/D");
    constantsTree->Branch("half_life", &half_life, "half_life/D");

    // Fill tree with a single entry
    constantsTree->Fill();
    constantsTree->Write();

    // Create a TTree for cross-section data
    TTree *crossSectionTree = new TTree("cross_section", "Cross-Section Data");

    // Define variables
    double E, sigma, Dsigma;

    // Set branches
    crossSectionTree->Branch("E", &E, "E/D");
    crossSectionTree->Branch("sigma", &sigma, "sigma/D");
    crossSectionTree->Branch("Dsigma", &Dsigma, "Dsigma/D");

    // Load cross-section data from file
    std::ifstream infile("cross_section.txt");
    if (!infile) {
        std::cerr << "Error: Cannot open cross_section.txt!" << std::endl;
        return;
    }

    while (infile >> E >> sigma >> Dsigma) {
        crossSectionTree->Fill();
    }
    infile.close();

    // Write cross-section data to ROOT file
    crossSectionTree->Write();

    // Cleanup
    file->Close();
    delete file;
}