#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;
  
int main(){
  std::string line;
  std::string name;
  int id;
  
  ifstream infile;
  infile.open("config.dat");
  
  while(getline(infile, line)){
    cout << "INFO: Line to parse: " << line << endl;
    stringstream ss(line);
    ss >> name >> id;
    cout << "INFO: Name:  " << name << endl;
    cout << "INFO: Id:    " << id << endl;
  }
  
  return 0;
}
