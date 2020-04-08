#include <TCanvas.h>
#include <TH1I.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"

using namespace std;

void printDelt_RE4R2CH20(){

    TTree *tree = new TTree("treedata","Data from csv file"); 

    // Uses the Read1File method: the first entry is the csv file name, the second entry contains the step voltage
    // applied in the RPC, the current of the first month and the current of the last month. The last entry is the
    // csv format (comma).
    tree->ReadFile("Chamber_RE-4_R2_CH20.csv","VMON/D,IMONJ/D,IMONF/D",','); 
    
    double IMONJ, IMONF, VMON;  

        
    // Setting the variables on tree branches.
    tree->SetBranchAddress("VMON", &VMON);
    tree->SetBranchAddress("IMONJ", &IMONJ);
    tree->SetBranchAddress("IMONF", &IMONF);

    // Number of tree entries, and arrays for voltage and current variation.
    const auto n = tree->GetEntries();
    double V[n];
    double deltaI[n];
    
    // Loop for take tree data.
    for (unsigned i = 0; i !=n; ++i){
        tree->GetEntry(i);
        V[i] = VMON;
        deltaI[i] = IMONJ-IMONF;
        
    
    }
    // Canvas and graph.
    TCanvas *canv = new TCanvas("Plot", "Current variation analysis", 700, 600);
    TGraph *gr = new TGraph(n, V, deltaI);
    gr ->Draw("ap");

    // Marker style (cross).
    gr->SetMarkerStyle(47);
    // Marker size (1 - default).
    gr->SetMarkerSize(1.2);
    // Marker color (kpurple=890).
    gr->SetMarkerColor(890);
    // Title and axis.
    gr->SetTitle("#DeltaI for RE-4_R2_CH20 - Days: 20.01.2020 and 17.02.2020;  V (V); #DeltaI (#muA)");

    // Linear functions for plotting.
    TF1 *pol1 = new TF1("pol1","pol1");
    // Keeping the fit results.
    TFitResultPtr fit = gr->Fit(pol1,"S M");
    // Fit parameters.
    gStyle->SetOptFit(1111);
    TPaveStats *st = (TPaveStats*) gr->FindObject("stats");
    st->SetX1NDC(0.1); //new x start position
    st->SetX2NDC(0.5); //new x end position
    fit->Print();
    gr->Draw("AP");
    canv->SaveAs("Current variation for RE-4_R2_CH20 .pdf");
    

}

void printDelt_RE4R2CH15(){

    TTree *tree = new TTree("treedata","Data from csv file"); 

    // Uses the Read1File method: the first entry is the csv file name, the second entry contains the step voltage
    // applied in the RPC, the current of the first month and the current of the last month. The last entry is the
    // csv format (comma).
    tree->ReadFile("Chamber_RE-4_R2_CH15.csv","VMON/D,IMONJ/D,IMONF/D",','); 
    
    double IMONJ, IMONF, VMON;  
    
    // Setting the variables on tree branches.
    tree->SetBranchAddress("VMON", &VMON);
    tree->SetBranchAddress("IMONJ", &IMONJ);
    tree->SetBranchAddress("IMONF", &IMONF);

    // Number of tree entries, and arrays for voltage and current variation.
    const auto n = tree->GetEntries();
    double V[n];
    double deltaI[n];
    

    for (unsigned i = 0; i !=n; ++i){
        tree->GetEntry(i);
        V[i] = VMON;
        deltaI[i] = IMONJ-IMONF;
        
    
    }
    // Criação do canvas.
    TCanvas *canv = new TCanvas("Plot", "Current variation analysis", 700, 600);
    TGraph *gr = new TGraph(n, V, deltaI);
    gr ->Draw("ap");

    // Marker style (cross).
    gr->SetMarkerStyle(47);
    // Marker size (1 - default).
    gr->SetMarkerSize(1.2);
    // Marker color (kpurple=890).
    gr->SetMarkerColor(890);
    // Title.
    gr->SetTitle("#DeltaI for RE-4_R2_CH15 - Days: 20.01.2020 and 17.02.2020;  V (V); #DeltaI (#muA)");

    // Linear functions for plotting.
    TF1 *pol1 = new TF1("pol1","pol1");
    // Keeping the fit results.
    TFitResultPtr fit = gr->Fit(pol1,"S M");
    // Fit parameters.
    gStyle->SetOptFit(1111);
    TPaveStats *st = (TPaveStats*) gr->FindObject("stats");
    st->SetX1NDC(0.1); //new x start position
    st->SetX2NDC(0.5); //new x end position
    fit->Print();
    gr->Draw("AP");
    canv->SaveAs("Current variation for RE-4_R2_CH15 .pdf");
    

}

void printDelt_RE4R2CH14(){

    TTree *tree = new TTree("treedata","Data from csv file"); 

    // Uses the Read1File method: the first entry is the csv file name, the second entry contains the step voltage
    // applied in the RPC, the current of the first month and the current of the last month. The last entry is the
    // csv format (comma).
    tree->ReadFile("Chamber_RE-4_R2_CH14.csv","VMON/D,IMONJ/D,IMONF/D",','); 
    
    double IMONJ, IMONF, VMON;  
    
    // Setting the variables on tree branches.
    tree->SetBranchAddress("VMON", &VMON);
    tree->SetBranchAddress("IMONJ", &IMONJ);
    tree->SetBranchAddress("IMONF", &IMONF);

    // Number of tree entries, and arrays for voltage and current variation.
    const auto n = tree->GetEntries();
    double V[n];
    double deltaI[n];
    

    for (unsigned i = 0; i !=n; ++i){
        tree->GetEntry(i);
        V[i] = VMON;
        deltaI[i] = IMONJ-IMONF;
        
    
    }
    // Criação do canvas.
    TCanvas *canv = new TCanvas("Plot", "Current variation analysis", 700, 600);
    TGraph *gr = new TGraph(n, V, deltaI);
    gr ->Draw("ap");

    // Marker style (cross).
    gr->SetMarkerStyle(47);
    // Marker size (1 - default).
    gr->SetMarkerSize(1.2);
    // Marker color (kpurple=890).
    gr->SetMarkerColor(890);
    // Title.
    gr->SetTitle("#DeltaI for RE-4_R2_CH14 - Days: 20.01.2020 and 17.02.2020;  V (V); #DeltaI (#muA)");

    // Linear functions for plotting.
    TF1 *pol1 = new TF1("pol1","pol1");
    // Keeping the fit results.
    TFitResultPtr fit = gr->Fit(pol1,"S M");
    // Fit parameters.
    gStyle->SetOptFit(1111);
    TPaveStats *st = (TPaveStats*) gr->FindObject("stats");
    st->SetX1NDC(0.1); //new x start position
    st->SetX2NDC(0.5); //new x end position
    fit->Print();
    gr->Draw("AP");
    canv->SaveAs("Current variation for RE-4_R2_CH14 .pdf");
    

}

void printDelt_RE4R2CH13(){

    TTree *tree = new TTree("treedata","Data from csv file"); 

    // Uses the Read1File method: the first entry is the csv file name, the second entry contains the step voltage
    // applied in the RPC, the current of the first month and the current of the last month. The last entry is the
    // csv format (comma).
    tree->ReadFile("Chamber_RE-4_R2_CH13.csv","VMON/D,IMONJ/D,IMONF/D",','); 
    
    double IMONJ, IMONF, VMON;  
    
    // Setting the variables on tree branches.
    tree->SetBranchAddress("VMON", &VMON);
    tree->SetBranchAddress("IMONJ", &IMONJ);
    tree->SetBranchAddress("IMONF", &IMONF);

    // Number of tree entries, and arrays for voltage and current variation.
    const auto n = tree->GetEntries();
    double V[n];
    double deltaI[n];
    

    for (unsigned i = 0; i !=n; ++i){
        tree->GetEntry(i);
        V[i] = VMON;
        deltaI[i] = IMONJ-IMONF;
        
    
    }
    // Criação do canvas.
    TCanvas *canv = new TCanvas("Plot", "Current variation analysis", 700, 600);
    TGraph *gr = new TGraph(n, V, deltaI);
    gr ->Draw("ap");

    // Marker style (cross).
    gr->SetMarkerStyle(47);
    // Marker size (1 - default).
    gr->SetMarkerSize(1.2);
    // Marker color (kpurple=890).
    gr->SetMarkerColor(890);
    // Title.
    gr->SetTitle("#DeltaI for RE-4_R2_CH13 - Days: 20.01.2020 and 17.02.2020;  V (V); #DeltaI (#muA)");

    // Linear functions for plotting.
    TF1 *pol1 = new TF1("pol1","pol1");
    // Keeping the fit results.
    TFitResultPtr fit = gr->Fit(pol1,"S M");
    // Fit parameters.
    gStyle->SetOptFit(1111);
    TPaveStats *st = (TPaveStats*) gr->FindObject("stats");
    st->SetX1NDC(0.1); //new x start position
    st->SetX2NDC(0.5); //new x end position
    fit->Print();
    gr->Draw("AP");
    canv->SaveAs("Current variation for RE-4_R2_CH13 .pdf");
    

}

void printDelt_RE4R2CH12(){

    TTree *tree = new TTree("treedata","Data from csv file"); 

    // Uses the Read1File method: the first entry is the csv file name, the second entry contains the step voltage
    // applied in the RPC, the current of the first month and the current of the last month. The last entry is the
    // csv format (comma).
    tree->ReadFile("Chamber_RE-4_R2_CH12.csv","VMON/D,IMONJ/D,IMONF/D",','); 
    
    double IMONJ, IMONF, VMON;  
    
    // Setting the variables on tree branches.
    tree->SetBranchAddress("VMON", &VMON);
    tree->SetBranchAddress("IMONJ", &IMONJ);
    tree->SetBranchAddress("IMONF", &IMONF);

    // Number of tree entries, and arrays for voltage and current variation.
    const auto n = tree->GetEntries();
    double V[n];
    double deltaI[n];
    

    for (unsigned i = 0; i !=n; ++i){
        tree->GetEntry(i);
        V[i] = VMON;
        deltaI[i] = IMONJ-IMONF;
        
    
    }
    // Criação do canvas.
    TCanvas *canv = new TCanvas("Plot", "Current variation analysis", 700, 600);
    TGraph *gr = new TGraph(n, V, deltaI);
    gr ->Draw("ap");

    // Marker style (cross).
    gr->SetMarkerStyle(47);
    // Marker size (1 - default).
    gr->SetMarkerSize(1.2);
    // Marker color (kpurple=890).
    gr->SetMarkerColor(890);
    // Title.
    gr->SetTitle("#DeltaI for RE-4_R2_CH12 - Days: 20.01.2020 and 17.02.2020;  V (V); #DeltaI (#muA)");

    // Linear functions for plotting.
    TF1 *pol1 = new TF1("pol1","pol1");
    // Keeping the fit results.
    TFitResultPtr fit = gr->Fit(pol1,"S M");
    // Fit parameters.
    gStyle->SetOptFit(1111);
    TPaveStats *st = (TPaveStats*) gr->FindObject("stats");
    st->SetX1NDC(0.1); //new x start position
    st->SetX2NDC(0.5); //new x end position
    fit->Print();
    gr->Draw("AP");
    canv->SaveAs("Current variation for RE-4_R2_CH12.pdf");
    

}

void printDelt_RE4R2CH10(){

    TTree *tree = new TTree("treedata","Data from csv file"); 

    // Uses the Read1File method: the first entry is the csv file name, the second entry contains the step voltage
    // applied in the RPC, the current of the first month and the current of the last month. The last entry is the
    // csv format (comma).
    tree->ReadFile("Chamber_RE-4_R2_CH10.csv","VMON/D,IMONJ/D,IMONF/D",','); 
    
    double IMONJ, IMONF, VMON;  
    
    // Setting the variables on tree branches.
    tree->SetBranchAddress("VMON", &VMON);
    tree->SetBranchAddress("IMONJ", &IMONJ);
    tree->SetBranchAddress("IMONF", &IMONF);

    // Number of tree entries, and arrays for voltage and current variation.
    const auto n = tree->GetEntries();
    double V[n];
    double deltaI[n];
    

    for (unsigned i = 0; i !=n; ++i){
        tree->GetEntry(i);
        V[i] = VMON;
        deltaI[i] = IMONJ-IMONF;
        
    
    }
    // Criação do canvas.
    TCanvas *canv = new TCanvas("Plot", "Current variation analysis", 700, 600);
    TGraph *gr = new TGraph(n, V, deltaI);
    gr ->Draw("ap");

    // Marker style (cross).
    gr->SetMarkerStyle(47);
    // Marker size (1 - default).
    gr->SetMarkerSize(1.2);
    // Marker color (kpurple=890).
    gr->SetMarkerColor(890);
    // Title.
    gr->SetTitle("#DeltaI for RE-4_R2_CH10 - Days: 20.01.2020 and 17.02.2020;  V (V); #DeltaI (#muA)");

    // Linear functions for plotting.
    TF1 *pol1 = new TF1("pol1","pol1");
    // Keeping the fit results.
    TFitResultPtr fit = gr->Fit(pol1,"S M");
    // Fit parameters.
    gStyle->SetOptFit(1111);
    TPaveStats *st = (TPaveStats*) gr->FindObject("stats");
    st->SetX1NDC(0.1); //new x start position
    st->SetX2NDC(0.5); //new x end position
    fit->Print();
    gr->Draw("AP");
    canv->SaveAs("Current variation for RE-4_R2_CH10.pdf");
    

}