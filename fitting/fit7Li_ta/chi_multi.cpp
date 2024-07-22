#include "fit.h"
#include <chrono>


int main(int argc, char*argv[])
{

  //const int intermediate_spin = 2;
  const int Npeaks = 1;
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


  const int n = 324;//64;


string files1[n] = {

"sim_ta_4654_050.root", "sim_ta_4654_052.root", "sim_ta_4654_054.root", "sim_ta_4654_056.root", "sim_ta_4654_058.root", "sim_ta_4654_060.root", "sim_ta_4654_062.root", "sim_ta_4654_064.root", "sim_ta_4654_066.root", "sim_ta_4654_068.root", "sim_ta_4654_070.root", "sim_ta_4654_072.root", "sim_ta_4654_074.root", "sim_ta_4654_076.root", "sim_ta_4654_078.root", "sim_ta_4654_080.root", "sim_ta_4654_082.root", "sim_ta_4654_084.root", "sim_ta_4654_086.root", "sim_ta_4654_088.root", "sim_ta_4654_090.root", "sim_ta_4654_92.root", "sim_ta_4654_94.root", "sim_ta_4654_96.root", "sim_ta_4654_98.root", "sim_ta_4654_100.root", "sim_ta_4654_102.root", "sim_ta_4654_104.root", "sim_ta_4654_106.root", "sim_ta_4654_108.root", "sim_ta_4654_110.root", "sim_ta_4654_112.root", "sim_ta_4654_114.root", "sim_ta_4654_116.root", "sim_ta_4654_118.root", "sim_ta_4654_120.root",

"sim_ta_4652_050.root", "sim_ta_4652_052.root", "sim_ta_4652_054.root", "sim_ta_4652_056.root", "sim_ta_4652_058.root", "sim_ta_4652_060.root", "sim_ta_4652_062.root", "sim_ta_4652_064.root", "sim_ta_4652_066.root", "sim_ta_4652_068.root", "sim_ta_4652_070.root", "sim_ta_4652_072.root", "sim_ta_4652_074.root", "sim_ta_4652_076.root", "sim_ta_4652_078.root", "sim_ta_4652_080.root", "sim_ta_4652_082.root", "sim_ta_4652_084.root", "sim_ta_4652_086.root", "sim_ta_4652_088.root", "sim_ta_4652_090.root", "sim_ta_4652_92.root", "sim_ta_4652_94.root", "sim_ta_4652_96.root", "sim_ta_4652_98.root", "sim_ta_4652_100.root", "sim_ta_4652_102.root", "sim_ta_4652_104.root", "sim_ta_4652_106.root", "sim_ta_4652_108.root", "sim_ta_4652_110.root", "sim_ta_4652_112.root", "sim_ta_4652_114.root", "sim_ta_4652_116.root", "sim_ta_4652_118.root", "sim_ta_4652_120.root",

"sim_ta_4650_050.root", "sim_ta_4650_052.root", "sim_ta_4650_054.root", "sim_ta_4650_056.root", "sim_ta_4650_058.root", "sim_ta_4650_060.root", "sim_ta_4650_062.root", "sim_ta_4650_064.root", "sim_ta_4650_066.root", "sim_ta_4650_068.root", "sim_ta_4650_070.root", "sim_ta_4650_072.root", "sim_ta_4650_074.root", "sim_ta_4650_076.root", "sim_ta_4650_078.root", "sim_ta_4650_080.root", "sim_ta_4650_082.root", "sim_ta_4650_084.root", "sim_ta_4650_086.root", "sim_ta_4650_088.root", "sim_ta_4650_090.root", "sim_ta_4650_92.root", "sim_ta_4650_94.root", "sim_ta_4650_96.root", "sim_ta_4650_98.root", "sim_ta_4650_100.root", "sim_ta_4650_102.root", "sim_ta_4650_104.root", "sim_ta_4650_106.root", "sim_ta_4650_108.root", "sim_ta_4650_110.root", "sim_ta_4650_112.root", "sim_ta_4650_114.root", "sim_ta_4650_116.root", "sim_ta_4650_118.root", "sim_ta_4650_120.root",

"sim_ta_4648_050.root", "sim_ta_4648_052.root", "sim_ta_4648_054.root", "sim_ta_4648_056.root", "sim_ta_4648_058.root", "sim_ta_4648_060.root", "sim_ta_4648_062.root", "sim_ta_4648_064.root", "sim_ta_4648_066.root", "sim_ta_4648_068.root", "sim_ta_4648_070.root", "sim_ta_4648_072.root", "sim_ta_4648_074.root", "sim_ta_4648_076.root", "sim_ta_4648_078.root", "sim_ta_4648_080.root", "sim_ta_4648_082.root", "sim_ta_4648_084.root", "sim_ta_4648_086.root", "sim_ta_4648_088.root", "sim_ta_4648_090.root", "sim_ta_4648_92.root", "sim_ta_4648_94.root", "sim_ta_4648_96.root", "sim_ta_4648_98.root", "sim_ta_4648_100.root", "sim_ta_4648_102.root", "sim_ta_4648_104.root", "sim_ta_4648_106.root", "sim_ta_4648_108.root", "sim_ta_4648_110.root", "sim_ta_4648_112.root", "sim_ta_4648_114.root", "sim_ta_4648_116.root", "sim_ta_4648_118.root", "sim_ta_4648_120.root",

"sim_ta_4646_050.root", "sim_ta_4646_052.root", "sim_ta_4646_054.root", "sim_ta_4646_056.root", "sim_ta_4646_058.root", "sim_ta_4646_060.root", "sim_ta_4646_062.root", "sim_ta_4646_064.root", "sim_ta_4646_066.root", "sim_ta_4646_068.root", "sim_ta_4646_070.root", "sim_ta_4646_072.root", "sim_ta_4646_074.root", "sim_ta_4646_076.root", "sim_ta_4646_078.root", "sim_ta_4646_080.root", "sim_ta_4646_082.root", "sim_ta_4646_084.root", "sim_ta_4646_086.root", "sim_ta_4646_088.root", "sim_ta_4646_090.root", "sim_ta_4646_092.root", "sim_ta_4646_094.root", "sim_ta_4646_096.root", "sim_ta_4646_098.root", "sim_ta_4646_100.root", "sim_ta_4646_102.root", "sim_ta_4646_104.root", "sim_ta_4646_106.root", "sim_ta_4646_108.root", "sim_ta_4646_110.root", "sim_ta_4646_112.root", "sim_ta_4646_114.root", "sim_ta_4646_116.root", "sim_ta_4646_118.root", "sim_ta_4646_120.root",

"sim_ta_4644_050.root", "sim_ta_4644_052.root", "sim_ta_4644_054.root", "sim_ta_4644_056.root", "sim_ta_4644_058.root", "sim_ta_4644_060.root", "sim_ta_4644_062.root", "sim_ta_4644_064.root", "sim_ta_4644_066.root", "sim_ta_4644_068.root", "sim_ta_4644_070.root", "sim_ta_4644_072.root", "sim_ta_4644_074.root", "sim_ta_4644_076.root", "sim_ta_4644_078.root", "sim_ta_4644_080.root", "sim_ta_4644_082.root", "sim_ta_4644_084.root", "sim_ta_4644_086.root", "sim_ta_4644_088.root", "sim_ta_4644_090.root", "sim_ta_4644_092.root", "sim_ta_4644_094.root", "sim_ta_4644_096.root", "sim_ta_4644_098.root", "sim_ta_4644_100.root", "sim_ta_4644_102.root", "sim_ta_4644_104.root", "sim_ta_4644_106.root", "sim_ta_4644_108.root", "sim_ta_4644_110.root", "sim_ta_4644_112.root", "sim_ta_4644_114.root", "sim_ta_4644_116.root", "sim_ta_4644_118.root", "sim_ta_4644_120.root",

"sim_ta_4642_050.root", "sim_ta_4642_052.root", "sim_ta_4642_054.root", "sim_ta_4642_056.root", "sim_ta_4642_058.root", "sim_ta_4642_060.root", "sim_ta_4642_062.root", "sim_ta_4642_064.root", "sim_ta_4642_066.root", "sim_ta_4642_068.root", "sim_ta_4642_070.root", "sim_ta_4642_072.root", "sim_ta_4642_074.root", "sim_ta_4642_076.root", "sim_ta_4642_078.root", "sim_ta_4642_080.root", "sim_ta_4642_082.root", "sim_ta_4642_084.root", "sim_ta_4642_086.root", "sim_ta_4642_088.root", "sim_ta_4642_090.root", "sim_ta_4642_092.root", "sim_ta_4642_094.root", "sim_ta_4642_096.root", "sim_ta_4642_098.root", "sim_ta_4642_100.root", "sim_ta_4642_102.root", "sim_ta_4642_104.root", "sim_ta_4642_106.root", "sim_ta_4642_108.root", "sim_ta_4642_110.root", "sim_ta_4642_112.root", "sim_ta_4642_114.root", "sim_ta_4642_116.root", "sim_ta_4642_118.root", "sim_ta_4642_120.root",

"sim_ta_4640_50.root", "sim_ta_4640_52.root", "sim_ta_4640_54.root", "sim_ta_4640_56.root", "sim_ta_4640_58.root", "sim_ta_4640_60.root", "sim_ta_4640_62.root", "sim_ta_4640_64.root", "sim_ta_4640_66.root", "sim_ta_4640_68.root", "sim_ta_4640_70.root", "sim_ta_4640_72.root", "sim_ta_4640_74.root", "sim_ta_4640_76.root", "sim_ta_4640_78.root", "sim_ta_4640_80.root", "sim_ta_4640_82.root", "sim_ta_4640_84.root", "sim_ta_4640_86.root", "sim_ta_4640_88.root", "sim_ta_4640_90.root", "sim_ta_4640_92.root", "sim_ta_4640_94.root", "sim_ta_4640_96.root", "sim_ta_4640_98.root", "sim_ta_4640_100.root", "sim_ta_4640_102.root", "sim_ta_4640_104.root", "sim_ta_4640_106.root", "sim_ta_4640_108.root", "sim_ta_4640_110.root", "sim_ta_4640_112.root", "sim_ta_4640_114.root", "sim_ta_4640_116.root", "sim_ta_4640_118.root", "sim_ta_4640_120.root",

"sim_ta_4638_50.root", "sim_ta_4638_52.root", "sim_ta_4638_54.root", "sim_ta_4638_56.root", "sim_ta_4638_58.root", "sim_ta_4638_60.root", "sim_ta_4638_62.root", "sim_ta_4638_64.root", "sim_ta_4638_66.root", "sim_ta_4638_68.root", "sim_ta_4638_70.root", "sim_ta_4638_72.root", "sim_ta_4638_74.root", "sim_ta_4638_76.root", "sim_ta_4638_78.root", "sim_ta_4638_80.root", "sim_ta_4638_82.root", "sim_ta_4638_84.root", "sim_ta_4638_86.root", "sim_ta_4638_88.root", "sim_ta_4638_90.root", "sim_ta_4638_92.root", "sim_ta_4638_94.root", "sim_ta_4638_96.root", "sim_ta_4638_98.root", "sim_ta_4638_100.root", "sim_ta_4638_102.root", "sim_ta_4638_104.root", "sim_ta_4638_106.root", "sim_ta_4638_108.root", "sim_ta_4638_110.root", "sim_ta_4638_112.root", "sim_ta_4638_114.root", "sim_ta_4638_116.root", "sim_ta_4638_118.root", "sim_ta_4638_120.root"


};




  std::ofstream chilog;
  chilog.open("record_grid.txt", std::ios_base::app); // append instead of overwrite  

  fit Fit(N,NN,Npeaks,para2,useit);

  double chisq_min;
  bool isStop = false;
  for (int j=0; j<n; j++)
  {
  cout << "before write" << endl;
    string filename1 = "/home/Li7star/li7sim/rootout/" + files1[j];

    Fit.read1new(filename1);
    cout << Fit.ND << endl;
    
    chisq_min = Fit.powell(para,xi,ftol);
    cout << "chisq_min = " << chisq_min << endl;

    chilog << chisq_min << "\t" << "f1:" << files1[j] << endl;
    cout << "after write" << endl;
  }

  chilog.close();

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

  //auto end = std::chrono::system_clock::now();
  //time_t end_time = std::chrono::system_clock::to_time_t(end);




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
