#include <fstream>
#include <iostream>
#include "minimizeND.h"
#include <cmath>
#include "TFile.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TH2S.h"
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace std;

class peak
{
 public:
  //~peak(){delete []y; delete []ys;}

  double y[500];
  double ys[500];
};


class fit:public minimizeND
{
 public:
  int Npeaks;
  fit(int,int,int,double*,bool*);
  ~fit();
  double functND(double*);
  double getValue(int,double*);
  void printValue(int,double*);
  double getBackground(int,double*);
  void print(double*);
  void plot(double*);
  void read2new(string,string);
  void read3new(string,string,string);
  void readeff();
  void effcorr(TH1I* h);

  int istart;
  int istop;

  string filename1;
  string filename2;
  string filename3;

  float xlow;
  float xhigh;
  int low;
  int high;

  int nn;

  float scale;
  int n;
 private:
  double eff[400];
  double x[400];
  double sx[400];
  double yexp[400];
  double syexp[400];


  double * para2;
  bool * useit;
  int NN;

  int count;
  double ymax;

  int solution;

  peak * Peak;
  int intermediate_spin;

};
