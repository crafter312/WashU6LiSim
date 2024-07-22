#include "fit.h"
#include <chrono>


int main(int argc, char*argv[])
{

  //const int intermediate_spin = 2;
  const int Npeaks = 2;
  int N = 0;
  int NN = 0;
  double one;
  bool buse;
  bool useit[30];
  double para2[30];
  ifstream fin;
  fin.open("para.dat");
  if (!fin.is_open()) 
  {
    cout << "could not open file para.dat" << endl;
    abort();
  }

  for (;;)
  {
    fin >> one >> buse;
    if (fin.eof()) break;
    if (fin.bad()) break;
    if (NN < Npeaks) para2[NN]=sqrt(one);
    else para2[NN] = one;

    //second variable in para.dat indicates whether we use the line for fit or not
    if (buse)
    {
      N++;
      useit[NN] = true;
    }
    else useit[NN] = false;
    NN++;
  } 

  fin.close();
  fin.clear();

  
  double * para = new double [N];

  //xi is a n x n matrix
  double ** xi;
  xi = new double* [N];
  for (int i=0;i<N;i++) xi[i] = new double [N];

  for (int i=0;i<N;i++)
  {
    for (int j=0;j<N;j++)
    {
      if (i==j) xi[i][j] = 1.;
      else xi[i][j] = 0.;
    }
  }

  double const ftol = .000001;

  int ii =0;
  for (int i=0;i<NN;i++)
  {
    if(useit[i])
    {
      para[ii] = para2[i];
      cout << "i " << i << "para[" << ii << "] " << para[ii] << endl;
      ii++;
    }
  }

  if (N!=ii) abort();

  fit Fit(N,NN,Npeaks,para2,useit);

  string filename1 = "/home/Li7star/li7sim/rootout/sim_p6He_11295_184.root";
  //string filename2 = "/home/Li7star/li7sim/rootout/sim_p6He_11600_410.root";
  string filename2 = "/home/Li7star/li7sim/rootout/sim_p6He_11680_290.root";
  //string filename3 = "/home/Li7star/li7sim/rootout/sim_line_10800_130_line.root";

  Fit.read2new(filename1, filename2);
  
  double chisq_min = Fit.powell(para,xi,ftol);


  ii = 0;
  for (int  i = 0;i<NN;i++)
  {
    if (useit[i])
    {  
      para2[i] = para[ii];
      ii++;
    }
  } 

  cout << "chisq_min = " << chisq_min << endl;
  for (int i=0;i<NN;i++) 
  {
    if (i < Npeaks) cout << i << " " <<  pow(para2[i],2) << endl;
    else cout << i << " " << para2[i] << endl; 
  }

  auto end = std::chrono::system_clock::now();
  time_t end_time = std::chrono::system_clock::to_time_t(end);

  std::ofstream chilog;
  chilog.open("record.txt", std::ios_base::app); // append instead of overwrite
  chilog << "computed at " << std::ctime(&end_time);
  chilog << "    chi2=" << chisq_min << endl;
  chilog << "        f1:" << Fit.filename1 << endl;
  chilog << "        f2:" << Fit.filename2 << endl;
  chilog.close();


  ofstream fileOut;
  fileOut.open("para2.dat");

  for (int i=0;i<NN;i++)
  {
    if (i < Npeaks) fileOut << pow(para2[i],2) << " " << useit[i] << endl;
    else fileOut << para2[i] << " " << useit[i] << endl;
  }
  fileOut.close();

  Fit.plot(para2);

  //Fit.printValue(92,para);  
}
