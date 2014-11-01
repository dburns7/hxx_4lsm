#ifndef CUTFLOW_TOOL_H
#define CUTFLOW_TOOL_H

#include <vector>
#include <string>

class cutflow_tool {
public:
   cutflow_tool();

   void increment(int stage, int sample, double weight);
   double count(int stage, int sample);
   void print(int stage);

   void add_sample_name(int sample, const std::string & sample_name); // optional
private:
   std::vector< std::string > sample_names;
   std::vector< std::vector<double> > evt_counts;
   std::vector< std::vector<double> > wsq_counts;
};

#endif
