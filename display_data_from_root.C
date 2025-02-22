#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <iostream>

void display_data_from_root() {
    // Open the ROOT file
    TFile *file = new TFile("beam_intensity.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open beam_intensity.root!" << std::endl;
        return;
    }

    // Access the "data" directory
    TDirectory *dataDir = (TDirectory*)file->Get("data");
    if (!dataDir) {
        std::cerr << "Error: Cannot find directory 'data'!" << std::endl;
        file->Close();
        return;
    }

    // Retrieve the constants tree
    TTree *constantsTree = (TTree*)dataDir->Get("constants");
    if (constantsTree) {
        std::cout << "\nExperimental Constants:" << std::endl;
        
        std::string *name = nullptr;
        double value;
        constantsTree->SetBranchAddress("name", &name);
        constantsTree->SetBranchAddress("value", &value);

        for (Long64_t i = 0; i < constantsTree->GetEntries(); i++) {
            constantsTree->GetEntry(i);
            std::cout << *name << ": " << value << std::endl;
        }
    } else {
        std::cerr << "Error: Cannot find tree 'constants'!" << std::endl;
    }

    // Retrieve the cross-section tree
    TTree *crossSectionTree = (TTree*)dataDir->Get("cross_section");
    if (crossSectionTree) {
        std::cout << "\nCross-Section Data:" << std::endl;
        
        double E, sigma, Dsigma;
        crossSectionTree->SetBranchAddress("E", &E);
        crossSectionTree->SetBranchAddress("sigma", &sigma);
        crossSectionTree->SetBranchAddress("Dsigma", &Dsigma);

        for (Long64_t i = 0; i < crossSectionTree->GetEntries(); i++) {
            crossSectionTree->GetEntry(i);
            std::cout << "E = " << E << " MeV, sigma = " << sigma << " barn, Dsigma = " << Dsigma << " barn" << std::endl;
        }
    } else {
        std::cerr << "Error: Cannot find tree 'cross_section'!" << std::endl;
    }

    // Close file
    file->Close();
    delete file;
}
