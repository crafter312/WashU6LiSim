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

  void readeff();
  void effcorr(TH1I* h);
  void read1new(string fn1);

  int istart;
  int istop;

  string filename1;

  float xlow;
  float xhigh;
  int low;
  int high;

  int nn;

  float scale;
  int n;
 private:
  double eff[500];
  double x[500];
  double sx[500];
  double yexp[500];
  double syexp[500];


  double * para2;
  bool * useit;
  int NN;

  int count;
  double ymax;

  int solution;

  peak * Peak;
  int intermediate_spin;

};
