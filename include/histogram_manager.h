#ifndef HISTOGRAM_MANAGER_H
#define HISTOGRAM_MANAGER_H

#include <vector>
#include <string>

class TH1F;
class auto_write;

class histogram_manager {
public:

   histogram_manager();
   histogram_manager(TH1F * hbase );
   histogram_manager(TH1F * hbase, const histogram_manager & hm);
   histogram_manager(TH1F * hbase, const histogram_manager & hm, auto_write & aw);

   // add to an automatic write service:
   void add_auto_write(auto_write & aw);

   // set the base histogram, for which each sample will have its own copy:
   void set_base_histogram(TH1F * hbase){ this->hbase = hbase; }

   // add a sample and corresponding sample name
   void add_sample         (int sample, const std::string & tag);

   // copy the samples and names from an already initialized histogram_manager
   void copy_samples_from  (const histogram_manager & hm);

   // find the index for this sample
   // (issues cout message and exits on no corresponding sample!)
   int index(int sample);

   // Fill histogram with this sample id:
   void Fill(int sample, double x);
   void Fill(int sample, double x, double w);

   void WriteAll();

   // direct access to everything for your enjoyment:
   TH1F * hbase;
   std::vector<int>          sample;
   std::vector<std::string>  tag;
   std::vector<TH1F *>       vhist;

private:
   void init(){ hbase = NULL; }
};


// Automatically write to the file list:
class auto_write {
public:
   auto_write();
   void AutoWrite(histogram_manager & hm);
   void WriteAll();
private:
   std::vector<histogram_manager *> write_list;
};



#endif
