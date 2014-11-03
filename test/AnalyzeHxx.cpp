#include <iostream>
#include <utility>
#include <vector>
#include <math.h>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include "koptions.h"
#include "histogram_manager.h"
#include "cutflow_tool.h"

using namespace std;
using namespace kutil;

void usage(){
   cout << "usage:  AnalyzeHxx  [OPTIONS] <input_file> <output_root> <output_dir>\n";
   cout << "\n";
   exit(0);
}


double sensitivity(TH1F * hsig, TH1F * hbkg, double sigtot){
  int n = hsig->GetNbinsX() + 2;
  double best_sens = 0.0;
  int ibest = 0;
  for (int i=0; i<n; i++){
    double s = 0;
    double b = 0;
    for (int j=i; j<n; j++){
      //s += hsig->GetBinContent(i);
      s += sigtot * hsig->GetBinContent(i);
      b += hbkg->GetBinContent(i);
    }    
    double sens = 0;
    if (b > 0.0) sens = s / sqrt(b);
    if (sens > best_sens) { 
       best_sens = sens;
       ibest = i;
    }
  }
  cout << "Best MET cut at:  " << hsig->GetBinLowEdge(ibest) << "\n";
  return best_sens;
}


int main(int argc, char *argv[])
{
   TRandom rng;
   rng.SetSeed(2014);
   
   std::string opt, infile, outroot, outdir;
   koptions options(argc, argv);
   
   //check for the --help option:
   if ( options.find("--help") ) { usage(); }

   //check for unrecognized options (beginning with -- or -)
   while(options.unrecognized_option(opt)){
      cout << "WARNING: unrecognized option:" << opt << "\n";
      usage();
   } 

   if (options.size() != 3){
      usage();
   }

   options.shift_argument(infile);
   options.shift_argument(outroot);
   options.shift_argument(outdir);

   cout << "INFO: reading analysis tree file:    " << infile  << "\n";
   cout << "INFO: writing histogram root file:   " << outroot << "\n";
   cout << "INFO: writing results to directory:  " << outdir  << "\n";

   auto_write aw;

   cutflow_tool cutflow;
   histogram_manager h0mll(new TH1F("h0mll","",60,60.0,120.0));

     h0mll.add_sample(9000, "_DAS");
     h0mll.add_sample(9001, "_PU");
     h0mll.add_sample(9002, "_NO_PU");
     
     cutflow.add_sample_name(9000, "DAS");
     cutflow.add_sample_name(9001, "PU");
     cutflow.add_sample_name(9002, "NO_PU");

     h0mll.add_sample(20, "_hxx_1GeV");
     h0mll.add_sample(21, "_hxx_10GeV");
     h0mll.add_sample(22, "_hxx_100GeV");
     h0mll.add_sample(23, "_hxx_500GeV");
     h0mll.add_sample(24, "_hxx_1000GeV");
     
     cutflow.add_sample_name(20, "hxx_1GeV"); 
     cutflow.add_sample_name(21, "hxx_10GeV"); 
     cutflow.add_sample_name(22, "hxx_100GeV"); 
     cutflow.add_sample_name(23, "hxx_500GeV"); 
     cutflow.add_sample_name(24, "hxx_1000GeV"); 

   
   h0mll.add_auto_write(aw);

   TH1F hfit_bkg("sample_1","",  100,0.0,300.0);
   TH1F hfit_sig1("signal1","",  100,0.0,300.0);
   TH1F hfit_sig10("signal10","",  100,0.0,300.0);
   TH1F hfit_sig100("signal100","",  100,0.0,300.0);
   TH1F hfit_sig500("signal500","",  100,0.0,300.0);
   TH1F hfit_sig1000("signal1000","",  100,0.0,300.0);

   int    nsig = 0;
   double wsig = 0.0;
   hfit_bkg.Sumw2();
   hfit_sig1.Sumw2();
   hfit_sig10.Sumw2();
   hfit_sig100.Sumw2();
   hfit_sig500.Sumw2();
   hfit_sig1000.Sumw2();
   
   //control plots (no pre-selection)
   histogram_manager hcnjet  (new TH1F("hcnjet",  "", 10,  0.0, 10.0),     h0mll, aw);
   histogram_manager hcmjjon (new TH1F("hcmjjon", "", 75,  0.0, 150.0),    h0mll, aw);
   histogram_manager hcmjjoff(new TH1F("hcmjjoff","", 100, 0.0, 500.0),    h0mll, aw);

   //test new variables
   histogram_manager htestx  (new TH1F("htestx",  "", 25,  0.0, 3.2),      h0mll, aw);
   histogram_manager htesty  (new TH1F("htesty",  "", 25,  0.0, 200.0),    h0mll, aw);
 
   cout << "INFO: opening file: " << infile << "\n";

   TFile * file = new TFile(infile.c_str());
   TTree * data = (TTree*) file->Get("SelectedTree");

   if (! data) {
      cout << "ERROR:  could not open tree.\n";
      exit(0);
   }
   
   long int numberOfEntries = data->GetEntries();

   data->Print();
  
   int count = 0;
   int nupdate = numberOfEntries / 20;
   if (nupdate < 1) nupdate=1;
   cout << "PROGRESS:  ";

   int prescale  = 0; 
   int iprescale = 0;
   // Loop over all events
   for(Int_t entry = 0; entry < numberOfEntries; ++entry)
   {
      count++;
      if (count % nupdate == 0) { cout << "."; cout.flush(); }

      if (iprescale < prescale){
         iprescale++;
         continue;
      } else {
         iprescale = 0;
      }
      
      data->GetEntry(entry);

      cutflow.increment(0, 1, 1);      
      //cutflow.increment(0, data.sample, data.weight);      

      // fill limit setting histograms:
      //if (data.sample < 20)  hfit_bkg.Fill(data.nopu_met, data.weight);
      //if (data.sample == 20) hfit_sig1.Fill(data.nopu_met, data.weight);
   }
   cout << "\n";

   cout << "SUMMARY:  read " << count << " of " << data->GetEntries() << " events from analysis tree.\n";
   
   TFile * foutroot = new TFile(outroot.c_str(), "RECREATE");
   foutroot->cd();
   aw.WriteAll();
   foutroot->Close();

   cout << "SUMMARY:  done writing files.\n";

   cout << "Cutflow:  Stage 0 \n";
   cutflow.print(0);

   char name[100];
   TFile * f = NULL;
   TH1F * h = NULL;

   sprintf(name, "%s/mchi1.root", outdir.c_str());
   f = new TFile(name, "RECREATE");
   f->cd();
   h = (TH1F *) hfit_sig1.Clone("signal");
   h->Write();
   hfit_bkg.Write();
   f->Close();

   sprintf(name, "%s/mchi10.root", outdir.c_str());
   f = new TFile(name, "RECREATE");
   f->cd();
   h = (TH1F *) hfit_sig10.Clone("signal");
   h->Write();
   hfit_bkg.Write();
   f->Close();

   sprintf(name, "%s/mchi100.root", outdir.c_str());
   f = new TFile(name, "RECREATE");
   f->cd();
   h = (TH1F *) hfit_sig100.Clone("signal");
   h->Write();
   hfit_bkg.Write();
   f->Close();

   sprintf(name, "%s/mchi500.root", outdir.c_str());
   f = new TFile(name, "RECREATE");
   f->cd();
   h = (TH1F *) hfit_sig500.Clone("signal");
   h->Write();
   hfit_bkg.Write();
   f->Close();

   sprintf(name, "%s/mchi1000.root", outdir.c_str());
   f = new TFile(name, "RECREATE");
   f->cd();
   h = (TH1F *) hfit_sig1000.Clone("signal");
   h->Write();
   hfit_bkg.Write();
   f->Close();

}


