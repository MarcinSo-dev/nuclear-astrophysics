#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <iostream>
#include <iomanip>

void display_data() {
    // Open the ROOT file
    TFile *file = new TFile("beam_intensity.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open beam_intensity.root!" << std::endl;
        return;
    }

    // Load Constants Tree
    TTree *constantsTree = (TTree*)file->Get("data/constants");
    if (!constantsTree) {
        std::cerr << "Error: Cannot find 'constants' tree!" << std::endl;
        file->Close();
        return;
    }

    std::cout << "\n📌 Stored Constants in ROOT File:\n";

    // Iterate over branches (constants)
    TObjArray *branches = constantsTree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); i++) {
        TBranch *branch = (TBranch*)branches->At(i);
        if (branch) {
            double value;
            constantsTree->SetBranchAddress(branch->GetName(), &value);
            constantsTree->GetEntry(0);
            std::cout << branch->GetName() << " = " << value << std::endl;
        }
    }

    // Load Cross-Section Tree
    TTree *crossSectionTree = (TTree*)file->Get("data/cross_section");
    if (!crossSectionTree) {
        std::cerr << "Error: Cannot find 'cross_section' tree!" << std::endl;
        file->Close();
        return;
    }

    std::cout << "\n📌 Stored Cross-Section Data in ROOT File:\n";
    std::cout << "E (MeV)  |  Sigma (barn)  |  Dsigma (barn)\n";
    std::cout << "-------------------------------------------\n";

    // Set branches
    double E, sigma, Dsigma;
    crossSectionTree->SetBranchAddress("E", &E);
    crossSectionTree->SetBranchAddress("sigma", &sigma);
    crossSectionTree->SetBranchAddress("Dsigma", &Dsigma);

    // Print all rows
    int nEntries = crossSectionTree->GetEntries();
    for (int i = 0; i < nEntries; i++) {
        crossSectionTree->GetEntry(i);
        std::cout << std::fixed << std::setprecision(3)
                  << E << "  |  " << sigma << "  |  " << Dsigma << std::endl;
    }

    // Cleanup
    file->Close();
    delete file;
}