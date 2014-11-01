#include <iostream>
#include <stdlib.h>

#include "TH1F.h"

#include "histogram_manager.h"

using namespace std;



histogram_manager::histogram_manager(){ init(); }
   
histogram_manager::histogram_manager(TH1F * hbase ){ 
   init(); set_base_histogram(hbase); 
}
   
histogram_manager::histogram_manager(TH1F * hbase, const histogram_manager & hm){ 
   init(); set_base_histogram(hbase); copy_samples_from(hm); 
}

histogram_manager::histogram_manager(TH1F * hbase, const histogram_manager & hm, auto_write & aw){ 
   init(); set_base_histogram(hbase); copy_samples_from(hm); add_auto_write(aw);
}

// add to an automatic write service:
void histogram_manager::add_auto_write(auto_write & aw){ aw.AutoWrite(*this); }


void histogram_manager::add_sample(int sample, const std::string & tag){
   string hname(hbase->GetName());
   hname += tag;
   
   if (! hbase){cout << "ERROR: no base histogram defined!\n"; exit(0); }

   this->sample.push_back (sample);
   this->tag.push_back    (tag);
   //cout << "INFO:  adding histogram named " << hname << " to managed collection.\n";
   TH1F * temp = (TH1F*) hbase->Clone(hname.c_str());
   temp->Sumw2();
   vhist.push_back((TH1F*) temp);

}

void histogram_manager::copy_samples_from(const histogram_manager & hm){
   int n = hm.sample.size();
   for (int i=0; i<n; i++){
      add_sample(hm.sample[i], hm.tag[i]);
   }
}

int histogram_manager::index(int sample){
   int n = this->sample.size();
   for (int i=0; i<n; i++){
      if (sample == this->sample[i]){
         return i;
      }
   }
   if (! hbase) { cout << "ERROR: no base histogram defined!\n"; exit(0); }
   cout << "ERROR:  undefined sample id " << sample << " requested for histogram manager " << hbase->GetName() << "\n";
   exit(0);
   return 0;
}

void histogram_manager::Fill(int sample, double x){
   vhist[index(sample)]->Fill(x);
}

void histogram_manager::Fill(int sample, double x, double w){
   vhist[index(sample)]->Fill(x, w);
}

void histogram_manager::WriteAll(){
   hbase->Write();
   int n = sample.size();
   for (int i=0; i<n; i++){
      vhist[i]->Write();
   }
}


auto_write::auto_write(){
   write_list.clear();
}

// Automatically write to the file list:
void auto_write::AutoWrite(histogram_manager & hm){ 
   write_list.push_back(&hm); 
}

void auto_write::WriteAll(){
   int n = write_list.size();
   for (int i=0; i<n; i++){
      write_list[i]->WriteAll();
   }
}


