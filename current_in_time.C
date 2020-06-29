#include "TCanvas.h"
#include "TTree.h"
#include "TGraph.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TF1.h"
#include "TStyle.h"
#include "date.h"

using namespace std;


// Function to plot read data from csv and generate graphs.
void dataRead(string chamber, int const nPoints, string dt[], double volt){ //TString chamber, int numDate, TString datemeas
   
    // Tree for read date from a csc file.
    TTree *tree = new TTree("Data for scan","Data from csv file"); 
    // Method to read the file .csv with two variables.
   
    tree->ReadFile("Current_Scan_AllData_2019_03_13_REAL.csv","CHAMBER_NAME_CMSSW/C:d1/C:d2/C:d3/C:vmon/D:imon/D",',');
    // Variables to put on the tree branch.
    double current, voltage;
    char chamb[20], date[20];
    // Setting the branch address with those variables.
    tree->SetBranchAddress("CHAMBER_NAME_CMSSW", &chamb);
    tree->SetBranchAddress("d1", &date);
    tree->SetBranchAddress("imon", &current);
    tree->SetBranchAddress("vmon", &voltage);
    // Creation of the arrays current and date for plot.
    const auto n = tree->GetEntries();
    double curr[nPoints];
    double tt[nPoints];
    // take the tree values and puts on the arrays.
    
    //std::cout << dt[] << std::endl;
            

    // Loop for fill time array.
    tt[0] = 0;
    
    std::string auxYear1, auxYear2;
    
    std::string auxMonth1, auxMonth2;
    
    std::string auxDay1, auxDay2;

    int year1, year2, month1, month2, day1, day2;

    // day, month, year.
    Date dt1 = {1, 2, 2000}; 
    Date dt2 = {1, 2, 2000}; 
    for (unsigned it = 1; it != nPoints; ++it ){
        // Conversion from string format "yyyy.mm.dd" to interger format yyyy, mm, dd
        auxYear1 = dt[it-1];
        auxYear1.resize(4);
        year1 = std::stoi(auxYear1);
        
        auxYear2 = dt[it];
        if (auxYear2 == "") continue;
        auxYear2.resize(4);
        year2 = std::stoi(auxYear2);
               
        auxMonth1 = dt[it-1];
        auxMonth1 = auxMonth1.substr(5,2);
        month1 = std::stoi(auxMonth1);

        auxMonth2 = dt[it];
        if (auxMonth2 == "") continue;        
        auxMonth2 = auxMonth2.substr(5,2);
        month2 = std::stoi(auxMonth2);
        
        auxDay1 = dt[it-1];
        auxDay1 = auxDay1.substr(8,2);
        day1 = std::stoi(auxDay1);

        auxDay2 = dt[it];
        auxDay2 = auxDay2.substr(8,2);
        day2 = std::stoi(auxDay2);

        // day, month, year.
        dt1 = {day1, month1, year1}; 
        dt2 = {day2, month2, year2};

        tt[it] = (getDifference(dt1, dt2) + tt[it-1]);
        
        //cout << tt[it]<< endl;
    }
        
   
    
    std::string auxDate = "";
    for (unsigned j = 0; j!=nPoints; ++j){
        //cout << dt[j] << endl;
        for (unsigned i = 0; i != n; ++i){
            tree->GetEntry(i);
            
            if(chamb == chamber){
                
                if(voltage == volt){ 
                    //std::cout << voltage << std::endl;
                    //std::cout << chamb << std::endl;
                    auxDate = date;
                    auxDate.resize(10);
                    
                       
                        
                        if (auxDate == dt[j]){
                            curr[j] = current;
                        }                   
                }
            }     
        }
    }
    
    cout << nPoints << endl;
    for (unsigned q = 0; q != nPoints; ++q){
        cout << "tempo: " << tt[q] << " e corrente: " << curr[q] << endl;
    }
    
    // Beautiful canvas
    TCanvas *c1 = new TCanvas("c1","Exponential behavior of current",200,10,700,500);
   
    // Graph with the current and date. 
    TGraph *gr1 = new TGraph(nPoints, tt, curr);
    gr1->Draw("AP");

    // Marker style (cross).
    gr1->SetMarkerStyle(47);
    // Marker size (1 - default).
    gr1->SetMarkerSize(1.2);
    // Marker color (kpurple=890).
    gr1->SetMarkerColor(890);
    // Title and axis.
    TString a = "Current for ";
    TString b = " at ";
    int kVolt = volt/1000;
    TString c = std::to_string(kVolt);
    TString d = " kV";
    gr1->SetTitle(a + chamber + b + c + d + ";  t (s); I (#muA) ");
        
    // const + exponential function for fitting.
    TF1 *sum = new TF1("sum", "[0]+[1]*TMath::Exp([2]*x)" );
    // Guess the parameters (If the fit goes bad you can do it on the fit panel and then improve this parameters)
    sum->SetParameters(-4, 10, -8e-9);
    // Keeping the fit results.
    TFitResultPtr fit = gr1->Fit(sum,"S M");
    // Fit parameters on canvas.
    gStyle->SetOptFit(1111);
    TPaveStats *st = (TPaveStats*) gr1->FindObject("stats");
    
    st->SetX1NDC(0.9); //new x start position
    st->SetX2NDC(0.5); //new x end position
    gr1->Draw("AP");
    TString cur = "Current in time for ";
    c1->SaveAs(cur + chamber + ".pdf");  
   // c1->Close();  
    
}

// Function to plot the RPC current in time.
/* void plotGraph(){
    dataRead("Chamber_RE+4_R2_CH04");
    dataRead("Chamber_RE+4_R2_CH05");
    dataRead("Chamber_RE+4_R2_CH06");
    dataRead("Chamber_RE+4_R2_CH07");
    dataRead("Chamber_RE+4_R2_CH08");
    dataRead("Chamber_RE+4_R2_CH09");
    

} */




int main(){

    int const nPoints = 3;
    std::string dats[nPoints] = {"2019.08.28","2019.09.02","2019.09.06"};
    
    
    dataRead("RE+4_R2_CH08", nPoints, dats, 6000);
    return 0;
}