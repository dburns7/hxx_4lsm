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
#include "TNtuple.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include "koptions.h"
#include "histogram_manager.h"
#include "cutflow_tool.h"

using namespace std;
using namespace kutil;

void usage(){
   cout << "usage:  AnalyzeHxx  [OPTIONS] <input_file> <output_root> <output_dir> <run_id>\n";
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
   
   std::string opt, infile, outroot, outdir, runname, runid;
   koptions options(argc, argv);
   
   //check for the --help option:
   if ( options.find("--help") ) { usage(); }

   //check for unrecognized options (beginning with -- or -)
   while(options.unrecognized_option(opt)){
      cout << "WARNING: unrecognized option:" << opt << "\n";
      usage();
   } 

   if (options.size() != 5){
      usage();
   }

   options.shift_argument(infile);
   options.shift_argument(outroot);
   options.shift_argument(outdir);
   options.shift_argument(runname);
   options.shift_argument(runid);

   cout << "INFO: Reading ntuple:                " << infile  << "\n";
   cout << "INFO: Using run name:                " << runname << "\n";
   cout << "INFO: Using run ID:                  " << runid   << "\n";
   cout << "INFO: Writing histogram root file:   " << outroot << "\n";
   cout << "INFO: Writing results to directory:  " << outdir  << "\n"; 

   auto_write aw;

   cutflow_tool cutflow;
   histogram_manager h0mll(new TH1F("h0mll","",60,60.0,120.0));
   
   char title[100];
   sprintf(title, "_%s", runname.c_str());
   int id = atoi(runid.c_str()); 

   cout << "TEST:  " << title << " " << id << endl;

   h0mll.add_sample(id, title);
   cutflow.add_sample_name(id, runname.c_str());
   
   /*h0mll.add_sample(0, "_data");
   cutflow.add_sample_name(0, "DATA");
   h0mll.add_sample(11113, "_2e2mu");
   cutflow.add_sample_name(11113, "2E2MU");
   h0mll.add_sample(11115, "_2e2tau");
   cutflow.add_sample_name(11115, "2E2TAU");
   h0mll.add_sample(11315, "_2mu2tau");
   cutflow.add_sample_name(11315, "2MU2TAU");
   h0mll.add_sample(11111, "_4e");
   cutflow.add_sample_name(11111, "4E");
   h0mll.add_sample(11313, "_4mu");
   cutflow.add_sample_name(11313, "4MU");
   h0mll.add_sample(11515, "_4tau");
   cutflow.add_sample_name(11515, "4TAU");
   */
   
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
   
   //Book histograms
   histogram_manager histZZMass     (new TH1F("histZZMass", "", 100, 70.0, 800.0), h0mll, aw);
   histogram_manager histPFMET      (new TH1F("histPFMET", "", 10, 0.0, 30.0), h0mll, aw);
   histogram_manager histPFMETctrl  (new TH1F("histPFMETctrl", "", 30, 0.0, 300.0), h0mll, aw);

   cout << "INFO: opening file: " << infile << "\n";
  
   TFile* file = NULL;
   TTree* data = NULL;
   if(id < 100){
     file = new TFile(infile.c_str());
     data = (TTree*) file->Get("SelectedTree");
   }
   else{
     file = new TFile(infile.c_str());
     TDirectoryFile * dir = (TDirectoryFile*) file->Get("demo");
     data = (TTree*) dir->Get("hxxtree");
   }

   if (! data) {
      cout << "ERROR:  could not open tree.\n";
      exit(0);
   }
   
   long int numberOfEntries = data->GetEntries();

   //data->Print();
   //data->ls();
  
   int count = 0;
   int nupdate = numberOfEntries / 20;
   if (nupdate < 1) nupdate=1;
   cout << "PROGRESS:  ";

   float ZZMass;
   short genProcessId;
   float PFMET;
   data->SetBranchAddress("ZZMass", &ZZMass);
   data->SetBranchAddress("genProcessId", &genProcessId);
   data->SetBranchAddress("PFMET", &PFMET);

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
      
      //if(PFMET > 75) continue; 

      if(entry < 20) cout << "TEST:  " << id << " " << ZZMass << endl;
      histZZMass .Fill(id, ZZMass, 0.000617832); 
      histPFMET  .Fill(id, PFMET, 0.000617832);
      
      if(ZZMass > 95 && ZZMass < 155) continue;

      histPFMETctrl.Fill(id, PFMET, 0.000617832);
      
      //cutflow.increment(0, 1, 1);      
      //cutflow.increment(0, data.sample, data.weight);      

      // fill limit setting histograms:
      //if (data.sample < 20)  hfit_bkg.Fill(data.nopu_met, data.weight);
      //if (data.sample == 20) hfit_sig1.Fill(data.nopu_met, data.weight);
   }
   cout << "\n";

   cout << "SUMMARY:  read " << count << " of " << data->GetEntries() << " events from analysis tree.\n";
   
     
   TFile * foutroot = new TFile(outroot.c_str(), "UPDATE");
   foutroot->cd();
   aw.WriteAll();
   foutroot->Close();

   cout << "SUMMARY:  done writing files.\n";

   //cout << "Cutflow:  Stage 0 \n";
   //cutflow.print(0);

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


