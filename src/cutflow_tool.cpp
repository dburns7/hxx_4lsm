#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "cutflow_tool.h"
#include "math.h"

using namespace std;

cutflow_tool::cutflow_tool(){
   evt_counts.clear();
}

void cutflow_tool::increment(int stage, int sample, double weight){
  if (stage >= evt_counts.size()) {
    evt_counts.resize(stage+1);
    wsq_counts.resize(stage+1);
  }
   std::vector<double> & x  = evt_counts[stage]; 
   std::vector<double> & x2 = wsq_counts[stage]; 
   if (sample >= x.size()){
     x.resize(sample+1,0);
     x2.resize(sample+1,0);
   }
   x[sample]  += weight;
   x2[sample] += weight*weight;

}

double cutflow_tool::count(int stage, int sample){
   if (stage >= evt_counts.size()) return 0;
   std::vector<double> & x = evt_counts[stage]; 
   if (sample >= x.size()) return 0;
   return x[sample];
}

void cutflow_tool::print(int stage){
   if (stage >= evt_counts.size()) return;   
   std::vector<double> & x  = evt_counts[stage]; 
   std::vector<double> & x2 = wsq_counts[stage]; 

   for (int i=0; i<x.size(); i++){
      double evts  = x[i];
      double devts = sqrt(x2[i]); 
      if (evts > 0.0) {
         if ((sample_names.size() > i) && (sample_names[i] != ""))
            cout << setw(10) << sample_names[i];
         else 
            cout << "sample " << setw(3) << i;
         cout << "    events:  " << evts << " +/-" << devts << "\n";
      }
   }

   double bkg = 0;
   for (int i=0; i<x.size(); i++){
     if (i < 20) bkg+= x[i];
   }
   cout << "Total Background:  " << bkg << "\n";

   if (bkg <= 0.0) return;
   
   for (int i=0; i<x.size(); i++){
     if (i >= 20){
       if ((sample_names.size() > i) && (sample_names[i] != ""))
	 cout << setw(10) << sample_names[i];
       else 
	 cout << "sample " << setw(3) << i;
       cout << " S / sqrt(B) = " << x[i] / sqrt(bkg) << "\n";
     }
   }



}

void cutflow_tool::add_sample_name(int sample, const std::string & sample_name){
   std::string empty = "";
   if (sample >= sample_names.size()) sample_names.resize(sample+1, empty);
   sample_names[sample] = sample_name;
}
