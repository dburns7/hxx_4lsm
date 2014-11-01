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
#include "hxx_tree.h"
#include "histogram_manager.h"
#include "cutflow_tool.h"

using namespace std;
using namespace kutil;

void usage(){
   cout << "usage:  AnalyzeHxx  [OPTIONS] <input_file> <output_root> <output_dir>\n";
   cout << "\n";
   cout << "  --met_smear=<x>   : extra amount to smear MET.\n";
   cout << "  --num_smear=<n>   : number of times to smear MET.\n";
   cout << "  --fake_rate=<x>   : apply simple (whole event) fake rate <x> to W+jets.\n";
   cout << "  --8tev            : 8 tev settings (default 14 tev).\n";
   cout << "  --pub_plots       : generate publication plots.\n";
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
   hxx_tree data;
   
   std::string opt, infile, outroot, outdir;
   koptions options(argc, argv);
   
   //check for the --help option:
   if ( options.find("--help") ) { usage(); }
   //options.set("--seed=", seed);   
   //if (seed > 0) rng.SetSeed(seed);
   double met_smear = 30.0;
   options.set("--met_smear=", met_smear);   
   int    num_smear = 1;
   options.set("--num_smear=", num_smear);   
   double fake_rate = 0.0;
   options.set("--fake_rate=", fake_rate);   
   int mode_8tev = options.find("--8tev");
   int pub_plots = options.find("--pub_plots");



   double lumi = 300;
   double met_cut = 150.0;

   if (mode_8tev){
     lumi = 20.0;
     met_cut = 75;
     cout << "INFO:  settings for sqrt(s) = 8 TeV\n";
   } else {
     cout << "INFO:  settings for sqrt(s) = 14 TeV\n";
   }
   cout << "INFO:  lumi:     " << lumi << "\n";
   cout << "INFO:  MET cut:  " << met_cut << "\n";

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
   cout << "INFO: MET smearing amount is " << met_smear << "\n";
   cout << "INFO: MET smearing number is " << num_smear << "\n";
   cout << "INFO: W+jets fake rate is " << fake_rate << "\n";

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
/*   TH1F hfit_sig21("signal21","",  100,0.0,300.0);
   TH1F hfit_sig22("signal22","",  100,0.0,300.0);
   TH1F hfit_sig23("signal23","",  100,0.0,300.0);
   TH1F hfit_sig24("signal24","",  100,0.0,300.0);
*/
   int    nsig = 0;
   double wsig = 0.0;
   hfit_bkg.Sumw2();
   hfit_sig1.Sumw2();
/*   hfit_sig21.Sumw2();
   hfit_sig22.Sumw2();
   hfit_sig23.Sumw2();
   hfit_sig24.Sumw2();
*/
/*
   
   //control plots (no pre-selection)
   histogram_manager hcnjet  (new TH1F("hcnjet",  "", 10,  0.0, 10.0),     h0mll, aw);
   histogram_manager hcmjjon (new TH1F("hcmjjon", "", 75,  0.0, 150.0),    h0mll, aw);
   histogram_manager hcmjjoff(new TH1F("hcmjjoff","", 100, 0.0, 500.0),    h0mll, aw);

   //test new variables
   histogram_manager htestx  (new TH1F("htestx",  "", 25,  0.0, 3.2),      h0mll, aw);
   histogram_manager htesty  (new TH1F("htesty",  "", 25,  0.0, 200.0),    h0mll, aw);
*/

   //histograms at stage 0
   histogram_manager h0ID       (new TH1F("h0ID",     "", 25,  0.0, 25.0), h0mll, aw);
   
   histogram_manager h0ept      (new TH1F("h0ept",    "", 100,  0.0, 250.0), h0mll, aw);
   histogram_manager h0eeta     (new TH1F("h0eeta",   "", 60,   -5.0, 5.0),  h0mll, aw);
   histogram_manager h0ephi     (new TH1F("h0ephi",   "", 60,   -5.0, 5.0),  h0mll, aw);
   histogram_manager h0ech      (new TH1F("h0ech",    "", 60,   -5.0, 5.0),  h0mll, aw);
   histogram_manager h0ecnt     (new TH1F("h0ecnt",   "", 10,   0.0, 10.0),  h0mll, aw);
   
   histogram_manager h0mpt      (new TH1F("h0mpt",    "", 100,  0.0, 250.0), h0mll, aw);
   histogram_manager h0meta     (new TH1F("h0meta",   "", 60,   -5.0, 5.0),  h0mll, aw);
   histogram_manager h0mphi     (new TH1F("h0mphi",   "", 60,   -5.0, 5.0),  h0mll, aw);
   histogram_manager h0mch      (new TH1F("h0mch",    "", 60,   -5.0, 5.0),  h0mll, aw);
   histogram_manager h0mcnt     (new TH1F("h0mcnt",   "", 10,   0.0, 10.0),  h0mll, aw);
   
   histogram_manager h0l1pt      (new TH1F("h0l1pt",    "", 200,  0.0, 250.0), h0mll, aw);
   histogram_manager h0l2pt      (new TH1F("h0l2pt",    "", 200,  0.0, 250.0), h0mll, aw);
   histogram_manager h0l3pt      (new TH1F("h0l3pt",    "", 200,  0.0, 250.0), h0mll, aw);
   histogram_manager h0l4pt      (new TH1F("h0l4pt",    "", 200,  0.0, 250.0), h0mll, aw);
   
   histogram_manager h0gpt      (new TH1F("h0gpt",    "", 100,  0.0, 250.0), h0mll, aw);
   histogram_manager h0geta     (new TH1F("h0geta",   "", 60,   -5.0, 5.0),  h0mll, aw);
   histogram_manager h0gphi     (new TH1F("h0gphi",   "", 60,   -5.0, 5.0),  h0mll, aw);
   histogram_manager h0gch      (new TH1F("h0gch",    "", 60,   -5.0, 5.0),  h0mll, aw);
   histogram_manager h0gcnt     (new TH1F("h0gcnt",   "", 10,   0.0, 10.0),  h0mll, aw);
   
   histogram_manager h0mee      (new TH1F("h0mee",    "", 100,  0.0, 250.0), h0mll, aw);
   histogram_manager h0mmm      (new TH1F("h0mmm",    "", 100,  0.0, 250.0), h0mll, aw);
   histogram_manager h0mgg      (new TH1F("h0mgg",    "", 100,  0.0, 250.0), h0mll, aw);
   histogram_manager h0mll_l    (new TH1F("h0mll_l",  "", 100,  0.0, 250.0), h0mll, aw);
   histogram_manager h0mll_sl   (new TH1F("h0mll_sl", "", 100,  0.0, 250.0), h0mll, aw);
   histogram_manager h0mllll    (new TH1F("h0mllll",  "", 100,  0.0, 250.0), h0mll, aw);
   
   histogram_manager h0met_nopu (new TH1F("h0met_nopu", "", 100,  0.0, 250.0), h0mll, aw);
   histogram_manager h0met_pu   (histprop(new TH1F("h0met_pu",   "", 100,  0.0, 250.0)), h0mll, aw);
   histogram_manager h0met_nopu_phi (new TH1F("h0met_nopu_phi", "", 60,   -5.0, 5.0),  h0mll, aw);
   histogram_manager h0met_pu_phi   (new TH1F("h0met_pu_phi",   "", 60,   -5.0, 5.0),  h0mll, aw);
   
   //histograms after analysis cuts
   histogram_manager h1met_nopu (new TH1F("h1met_nopu", "", 100,  0.0, 250.0), h0mll, aw);
 
   cout << "INFO: opening file: " << infile << "\n";

   TFile * file = new TFile(infile.c_str());
   TDirectoryFile * dir = (TDirectoryFile*) file->Get("demo");
   //TDirectoryFile	demo;1	demo
   TTree * tree = (TTree*) dir->Get("hxxtree");

   if (! tree) {
      cout << "ERROR:  could not open tree.\n";
      exit(0);
   }

   data.ReadTree(tree);
   long int numberOfEntries = tree->GetEntries();

   //tree->Print();
  
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
      
      tree->GetEntry(entry);

      h0ID.Fill(data.sample, data.sample);

      cutflow.increment(0, data.sample, data.weight);      
      h0mll.Fill(data.sample, data.mll, data.weight);
      
     //Fill single lepton hists at stage 0 
     for(int i=0; i<data.nelec; i++){
        //cout << data.nelec << " " << data.elec_pt->size() << " " << data.sample << endl;
        h0ept  .Fill(data.sample, data.elec_pt->at(i),     data.weight);
        h0eeta .Fill(data.sample, data.elec_eta->at(i),    data.weight);
        h0ephi .Fill(data.sample, data.elec_phi->at(i),    data.weight);
        h0ech  .Fill(data.sample, data.elec_ch->at(i),     data.weight);
      }
      h0ecnt   .Fill(data.sample, data.nelec,    data.weight);
      for(int i=0; i<data.nmuon; i++){
        h0mpt  .Fill(data.sample, data.muon_pt->at(i),     data.weight);
        h0meta .Fill(data.sample, data.muon_eta->at(i),    data.weight);
        h0mphi .Fill(data.sample, data.muon_phi->at(i),    data.weight);
        h0mch  .Fill(data.sample, data.muon_ch->at(i),     data.weight);
      }
      h0mcnt   .Fill(data.sample, data.nmuon,    data.weight);

     //Fill single photon hists at stage 0 
     for(int i=0; i<data.nphot; i++){
        h0gpt  .Fill(data.sample, data.phot_pt->at(i),     data.weight);
        h0geta .Fill(data.sample, data.phot_eta->at(i),    data.weight);
        h0gphi .Fill(data.sample, data.phot_phi->at(i),    data.weight);
      }
      h0gcnt   .Fill(data.sample, data.nphot,    data.weight);
      
      //Fill multiparticle hists at stage 0
      double mll_l, mll_sl, mllll;
      if(data.nelec == 4 && data.nmuon == 0){
        double mee1, mee2, mee_l, mee_sl, meeee, ept1, ept2, ept3, ept4;
        kinematic_vars_4e(data, mee1, mee2, mee_l, mee_sl, meeee, ept1, ept2, ept3, ept4);
        h0l1pt   .Fill(data.sample, ept1,   data.weight);
        h0l2pt   .Fill(data.sample, ept2,   data.weight);
        h0l3pt   .Fill(data.sample, ept3,   data.weight);
        h0l4pt   .Fill(data.sample, ept4,   data.weight);
        h0mee    .Fill(data.sample, mee1,   data.weight);
        h0mee    .Fill(data.sample, mee2,   data.weight);
        h0mll_l  .Fill(data.sample, mee_l,  data.weight);
        h0mll_sl .Fill(data.sample, mee_sl, data.weight);
        h0mllll  .Fill(data.sample, meeee,  data.weight);
        mll_l  = mee_l;
        mll_sl = mee_sl;
        mllll  = meeee;
      } 
      if(data.nelec == 0 && data.nmuon == 4){
        double mmm1, mmm2, mmm_l, mmm_sl, mmmmm;
        kinematic_vars_4m(data, mmm1, mmm2, mmm_l, mmm_sl, mmmmm);
        h0mmm    .Fill(data.sample, mmm1,   data.weight);
        h0mmm    .Fill(data.sample, mmm2,   data.weight);
        h0mll_l  .Fill(data.sample, mmm_l,  data.weight);
        h0mll_sl .Fill(data.sample, mmm_sl, data.weight);
        h0mllll  .Fill(data.sample, mmmmm,  data.weight);
        mll_l  = mmm_l;
        mll_sl = mmm_sl;
        mllll  = mmmmm;
      }
      if(data.nelec == 2 && data.nmuon == 2){
        double mee, mmm, m_l, m_sl, meemm;
        kinematic_vars_2e2m(data, mee, mmm, m_l, m_sl, meemm);
        h0mee    .Fill(data.sample, mee,    data.weight);
        h0mmm    .Fill(data.sample, mmm,    data.weight);
        h0mll_l  .Fill(data.sample, m_l,    data.weight);
        h0mll_sl .Fill(data.sample, m_sl,   data.weight);
        h0mllll  .Fill(data.sample, meemm,  data.weight);
        mll_l  = m_l;
        mll_sl = m_sl;
        mllll  = meemm;
      }
      double mgg;
      if(data.nphot == 2){
        kinematic_vars_2g(data, mgg);
	h0mgg    .Fill(data.sample, mgg,    data.weight);
      }

      //Fill other hists at stage 0
      h0met_nopu     .Fill(data.sample, data.nopu_met, data.weight);
      h0met_nopu_phi .Fill(data.sample, data.nopu_met_phi, data.weight);


      cutflow.increment(1, data.sample, data.weight);
      
      // Leading dilepton mass cut:
      if(mll_l < 50.0 || mll_l > 106.0) continue;      
      cutflow.increment(2, data.sample, data.weight);
      
      // Subleading dilepton mass cut:
      if(mll_sl < 12.0 || mll_sl > 115.0) continue;
      cutflow.increment(3, data.sample, data.weight);
    
      // Four lepton mass cut:
      if(mllll < 105.0 || mllll > 145.0) continue;
      cutflow.increment(4, data.sample, data.weight);

      // Met cut:
      if(data.nopu_met < 75.0) continue;
      cutflow.increment(5, data.sample, data.weight);

    
      //Fill histograms after analysis cuts
      h1met_nopu. Fill(data.sample, data.nopu_met, data.weight);

      // fill limit setting histograms:
      if (data.sample < 20)  hfit_bkg.Fill(data.nopu_met, data.weight);
      if (data.sample == 20) hfit_sig1.Fill(data.nopu_met, data.weight);

   }
   cout << "\n";

   cout << "SUMMARY:  read " << count << " of " << tree->GetEntries() << " events from analysis tree.\n";
   
   TFile * foutroot = new TFile(outroot.c_str(), "RECREATE");
   foutroot->cd();
   aw.WriteAll();
   foutroot->Close();

   //cout << "SUMMARY:  done writing files.\n";


////////////////////// 4l //////////////////////////
   cout << "Cutflow:  Stage 0 (llll preselection)\n";
   cutflow.print(0);

   cout << "Cutflow:  Stage 1 (after single lepton cuts)\n";
   cutflow.print(1);

   cout << "Cutflow:  Stage 2 (after mll_l cut) \n";
   cutflow.print(2);

   cout << "Cutflow:  Stage 3 (after mll_sl cut) \n";
   cutflow.print(3);
    
   cout << "Cutflow:  Stage 4 (after mllll cut) \n";
   cutflow.print(4);

   cout << "Cutflow:  Stage 5 (after met cut) \n";
   cutflow.print(5);

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

   /*sprintf(name, "%s/mchi10.root", outdir.c_str());
   f = new TFile(name, "RECREATE");
   f->cd();
   h = (TH1F *) hfit_sig21.Clone("signal");
   h->Write();
   hfit_bkg.Write();
   f->Close();

   sprintf(name, "%s/mchi100.root", outdir.c_str());
   f = new TFile(name, "RECREATE");
   f->cd();
   h = (TH1F *) hfit_sig22.Clone("signal");
   h->Write();
   hfit_bkg.Write();
   f->Close();

   sprintf(name, "%s/mchi500.root", outdir.c_str());
   f = new TFile(name, "RECREATE");
   f->cd();
   h = (TH1F *) hfit_sig23.Clone("signal");
   h->Write();
   hfit_bkg.Write();
   f->Close();

   sprintf(name, "%s/mchi1000.root", outdir.c_str());
   f = new TFile(name, "RECREATE");
   f->cd();
   h = (TH1F *) hfit_sig24.Clone("signal");
   h->Write();
   hfit_bkg.Write();
   f->Close();
*/


}


