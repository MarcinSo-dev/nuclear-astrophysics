#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

void import_data_from_txt() {
    // Create and open the ROOT file
    TFile *file = new TFile("beam_intensity.root", "RECREATE");

    // Create directories
    TDirectory *dataDir = file->mkdir("data");
    dataDir->cd();

    // Create a TTree for constants
    TTree *constantsTree = new TTree("constants", "Experimental Constants");

    // Define variables
    std::string name;
    double value;

    // Set branches
    constantsTree->Branch("name", &name);
    constantsTree->Branch("value", &value, "value/D");

    // Load constants from file
    std::ifstream constFile("constants.txt");
    if (!constFile) {
        std::cerr << "Error: Cannot open constants.txt!" << std::endl;
        return;
    }

    std::string line;
    while (std::getline(constFile, line)) {
        std::istringstream iss(line);
        iss >> name >> value; // Read name and value
        constantsTree->Fill();
    }
    constFile.close();

    // Write constants data to ROOT file
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
    
    // Skip the first row
    std::string dummyLine;
    std::getline(infile, dummyLine);

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