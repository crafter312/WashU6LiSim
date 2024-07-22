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


  const int n = 1;
  const int m = 435;

string files1[n] = {
 "sim_p6He_11295_184.root"};

//"sim_p6He_11295_154_wide1.root", "sim_p6He_11295_156_wide1.root", "sim_p6He_11295_158_wide1.root", "sim_p6He_11295_160_wide1.root", "sim_p6He_11295_162_wide1.root", "sim_p6He_11295_164_wide1.root", "sim_p6He_11295_166_wide1.root", "sim_p6He_11295_168_wide1.root", "sim_p6He_11295_170_wide1.root", "sim_p6He_11295_172_wide1.root", "sim_p6He_11295_174_wide1.root", "sim_p6He_11295_176_wide1.root", "sim_p6He_11295_178_wide1.root", "sim_p6He_11295_180_wide1.root", "sim_p6He_11295_182_wide1.root", "sim_p6He_11295_184_wide1.root", "sim_p6He_11295_186_wide1.root", "sim_p6He_11295_188_wide1.root", "sim_p6He_11295_190_wide1.root", "sim_p6He_11295_192_wide1.root", "sim_p6He_11295_194_wide1.root", "sim_p6He_11295_196_wide1.root", "sim_p6He_11295_198_wide1.root", "sim_p6He_11295_200_wide1.root", "sim_p6He_11295_202_wide1.root", "sim_p6He_11295_204_wide1.root", "sim_p6He_11295_206_wide1.root", "sim_p6He_11295_208_wide1.root", "sim_p6He_11295_210_wide1.root", "sim_p6He_11295_212_wide1.root", "sim_p6He_11295_214_wide1.root", "sim_p6He_11295_216_wide1.root",

//"sim_p6He_11300_154_reg2.root", "sim_p6He_11300_156_reg2.root", "sim_p6He_11300_158_reg2.root", "sim_p6He_11300_160_reg2.root", "sim_p6He_11300_162_reg2.root", "sim_p6He_11300_164_reg2.root", "sim_p6He_11300_166_reg2.root", "sim_p6He_11300_168_reg2.root", "sim_p6He_11300_170_reg2.root", "sim_p6He_11300_172_reg2.root", "sim_p6He_11300_174_reg2.root", "sim_p6He_11300_176_reg2.root", "sim_p6He_11300_178_reg2.root", "sim_p6He_11300_180_reg2.root", "sim_p6He_11300_182_reg2.root", "sim_p6He_11300_184_reg2.root", "sim_p6He_11300_186_reg2.root", "sim_p6He_11300_188_reg2.root", "sim_p6He_11300_190_reg2.root", "sim_p6He_11300_192_reg2.root", "sim_p6He_11300_194_reg2.root", "sim_p6He_11300_196_reg2.root", "sim_p6He_11300_198_reg2.root", "sim_p6He_11300_200_reg2.root", "sim_p6He_11300_202_reg2.root", "sim_p6He_11300_204_reg2.root", "sim_p6He_11300_206_reg2.root", "sim_p6He_11300_208_reg2.root", "sim_p6He_11300_210_reg2.root", "sim_p6He_11300_212_reg2.root", "sim_p6He_11300_214_reg2.root", "sim_p6He_11300_216_reg2.root",

//"sim_p6He_11300_154_reg3.root", "sim_p6He_11300_156_reg3.root", "sim_p6He_11300_158_reg3.root", "sim_p6He_11300_160_reg3.root", "sim_p6He_11300_162_reg3.root", "sim_p6He_11300_164_reg3.root", "sim_p6He_11300_166_reg3.root", "sim_p6He_11300_168_reg3.root", "sim_p6He_11300_170_reg3.root", "sim_p6He_11300_172_reg3.root", "sim_p6He_11300_174_reg3.root", "sim_p6He_11300_176_reg3.root", "sim_p6He_11300_178_reg3.root", "sim_p6He_11300_180_reg3.root", "sim_p6He_11300_182_reg3.root", "sim_p6He_11300_184_reg3.root", "sim_p6He_11300_186_reg3.root", "sim_p6He_11300_188_reg3.root", "sim_p6He_11300_190_reg3.root", "sim_p6He_11300_192_reg3.root", "sim_p6He_11300_194_reg3.root", "sim_p6He_11300_196_reg3.root", "sim_p6He_11300_198_reg3.root", "sim_p6He_11300_200_reg3.root", "sim_p6He_11300_202_reg3.root", "sim_p6He_11300_204_reg3.root", "sim_p6He_11300_206_reg3.root", "sim_p6He_11300_208_reg3.root", "sim_p6He_11300_210_reg3.root", "sim_p6He_11300_212_reg3.root", "sim_p6He_11300_214_reg3.root", "sim_p6He_11300_216_reg3.root",

//"sim_p6He_11305_154_reg2.root", "sim_p6He_11305_156_reg2.root", "sim_p6He_11305_158_reg2.root", "sim_p6He_11305_160_reg2.root", "sim_p6He_11305_162_reg2.root", "sim_p6He_11305_164_reg2.root", "sim_p6He_11305_166_reg2.root", "sim_p6He_11305_168_reg2.root", "sim_p6He_11305_170_reg2.root", "sim_p6He_11305_172_reg2.root", "sim_p6He_11305_174_reg2.root", "sim_p6He_11305_176_reg2.root", "sim_p6He_11305_178_reg2.root", "sim_p6He_11305_180_reg2.root", "sim_p6He_11305_182_reg2.root", "sim_p6He_11305_184_reg2.root", "sim_p6He_11305_186_reg2.root", "sim_p6He_11305_188_reg2.root", "sim_p6He_11305_190_reg2.root", "sim_p6He_11305_192_reg2.root", "sim_p6He_11305_194_reg2.root", "sim_p6He_11305_196_reg2.root", "sim_p6He_11305_198_reg2.root", "sim_p6He_11305_200_reg2.root", "sim_p6He_11305_202_reg2.root", "sim_p6He_11305_204_reg2.root", "sim_p6He_11305_206_reg2.root", "sim_p6He_11305_208_reg2.root", "sim_p6He_11305_210_reg2.root", "sim_p6He_11305_212_reg2.root", "sim_p6He_11305_214_reg2.root", "sim_p6He_11305_216_reg2.root",

//"sim_p6He_11305_154_reg3.root", "sim_p6He_11305_156_reg3.root", "sim_p6He_11305_158_reg3.root", "sim_p6He_11305_160_reg3.root", "sim_p6He_11305_162_reg3.root", "sim_p6He_11305_164_reg3.root", "sim_p6He_11305_166_reg3.root", "sim_p6He_11305_168_reg3.root", "sim_p6He_11305_170_reg3.root", "sim_p6He_11305_172_reg3.root", "sim_p6He_11305_174_reg3.root", "sim_p6He_11305_176_reg3.root", "sim_p6He_11305_178_reg3.root", "sim_p6He_11305_180_reg3.root", "sim_p6He_11305_182_reg3.root", "sim_p6He_11305_184_reg3.root", "sim_p6He_11305_186_reg3.root", "sim_p6He_11305_188_reg3.root", "sim_p6He_11305_190_reg3.root", "sim_p6He_11305_192_reg3.root", "sim_p6He_11305_194_reg3.root", "sim_p6He_11305_196_reg3.root", "sim_p6He_11305_198_reg3.root", "sim_p6He_11305_200_reg3.root", "sim_p6He_11305_202_reg3.root", "sim_p6He_11305_204_reg3.root", "sim_p6He_11305_206_reg3.root", "sim_p6He_11305_208_reg3.root", "sim_p6He_11305_210_reg3.root", "sim_p6He_11305_212_reg3.root", "sim_p6He_11305_214_reg3.root", "sim_p6He_11305_216_reg3.root"
//};













//"sim_p6He_11300_154_reg3.root", "sim_p6He_11300_156_reg3.root", "sim_p6He_11300_158_reg3.root", "sim_p6He_11300_160_reg3.root", "sim_p6He_11300_162_reg3.root", "sim_p6He_11300_164_reg3.root", "sim_p6He_11300_166_reg3.root", "sim_p6He_11300_168_reg3.root", "sim_p6He_11300_170_reg3.root", "sim_p6He_11300_172_reg3.root", "sim_p6He_11300_174_reg3.root", "sim_p6He_11300_176_reg3.root", "sim_p6He_11300_178_reg3.root", "sim_p6He_11300_180_reg3.root", "sim_p6He_11300_182_reg3.root", "sim_p6He_11300_184_reg3.root", "sim_p6He_11300_186_reg3.root", "sim_p6He_11300_188_reg3.root", "sim_p6He_11300_190_reg3.root", "sim_p6He_11300_192_reg3.root", "sim_p6He_11300_194_reg3.root", "sim_p6He_11300_196_reg3.root", "sim_p6He_11300_198_reg3.root", "sim_p6He_11300_200_reg3.root", "sim_p6He_11300_202_reg3.root", "sim_p6He_11300_204_reg3.root", "sim_p6He_11300_206_reg3.root", "sim_p6He_11300_208_reg3.root", "sim_p6He_11300_210_reg3.root", "sim_p6He_11300_212_reg3.root", "sim_p6He_11300_214_reg3.root", "sim_p6He_11300_216_reg3.root"





string files2[m] = {//"sim_p6He_11660_340.root"};

"sim_p6He_11560_220.root", "sim_p6He_11560_230.root", "sim_p6He_11560_240.root", "sim_p6He_11560_250.root", "sim_p6He_11560_260.root", "sim_p6He_11560_270.root", "sim_p6He_11560_280.root", "sim_p6He_11560_290.root", "sim_p6He_11560_300.root", "sim_p6He_11560_310.root", "sim_p6He_11560_320.root", "sim_p6He_11560_330.root", "sim_p6He_11560_340.root", "sim_p6He_11560_350.root", "sim_p6He_11560_360.root", "sim_p6He_11560_370.root", "sim_p6He_11560_380.root", "sim_p6He_11560_390.root", "sim_p6He_11560_400.root", "sim_p6He_11560_410.root", "sim_p6He_11560_420.root", "sim_p6He_11560_430.root", "sim_p6He_11560_440.root", "sim_p6He_11560_450.root", "sim_p6He_11560_460.root", "sim_p6He_11560_470.root", "sim_p6He_11560_480.root", "sim_p6He_11560_490.root", "sim_p6He_11560_500.root",

"sim_p6He_11570_220.root", "sim_p6He_11570_230.root", "sim_p6He_11570_240.root", "sim_p6He_11570_250.root", "sim_p6He_11570_260.root", "sim_p6He_11570_270.root", "sim_p6He_11570_280.root", "sim_p6He_11570_290.root", "sim_p6He_11570_300.root", "sim_p6He_11570_310.root", "sim_p6He_11570_320.root", "sim_p6He_11570_330.root", "sim_p6He_11570_340.root", "sim_p6He_11570_350.root", "sim_p6He_11570_360.root", "sim_p6He_11570_370.root", "sim_p6He_11570_380.root", "sim_p6He_11570_390.root", "sim_p6He_11570_400.root", "sim_p6He_11570_410.root", "sim_p6He_11570_420.root", "sim_p6He_11570_430.root", "sim_p6He_11570_440.root", "sim_p6He_11570_450.root", "sim_p6He_11570_460.root", "sim_p6He_11570_470.root", "sim_p6He_11570_480.root", "sim_p6He_11570_490.root", "sim_p6He_11570_500.root",

"sim_p6He_11580_220.root", "sim_p6He_11580_230.root", "sim_p6He_11580_240.root", "sim_p6He_11580_250.root", "sim_p6He_11580_260.root", "sim_p6He_11580_270.root", "sim_p6He_11580_280.root", "sim_p6He_11580_290.root", "sim_p6He_11580_300.root", "sim_p6He_11580_310.root", "sim_p6He_11580_320.root", "sim_p6He_11580_330.root", "sim_p6He_11580_340.root", "sim_p6He_11580_350.root", "sim_p6He_11580_360.root", "sim_p6He_11580_370.root", "sim_p6He_11580_380.root", "sim_p6He_11580_390.root", "sim_p6He_11580_400.root", "sim_p6He_11580_410.root", "sim_p6He_11580_420.root", "sim_p6He_11580_430.root", "sim_p6He_11580_440.root", "sim_p6He_11580_450.root", "sim_p6He_11580_460.root", "sim_p6He_11580_470.root", "sim_p6He_11580_480.root", "sim_p6He_11580_490.root", "sim_p6He_11580_500.root",

"sim_p6He_11590_220.root", "sim_p6He_11590_230.root", "sim_p6He_11590_240.root", "sim_p6He_11590_250.root", "sim_p6He_11590_260.root", "sim_p6He_11590_270.root", "sim_p6He_11590_280.root", "sim_p6He_11590_290.root", "sim_p6He_11590_300.root", "sim_p6He_11590_310.root", "sim_p6He_11590_320.root", "sim_p6He_11590_330.root", "sim_p6He_11590_340.root", "sim_p6He_11590_350.root", "sim_p6He_11590_360.root", "sim_p6He_11590_370.root", "sim_p6He_11590_380.root", "sim_p6He_11590_390.root", "sim_p6He_11590_400.root", "sim_p6He_11590_410.root", "sim_p6He_11590_420.root", "sim_p6He_11590_430.root", "sim_p6He_11590_440.root", "sim_p6He_11590_450.root", "sim_p6He_11590_460.root", "sim_p6He_11590_470.root", "sim_p6He_11590_480.root", "sim_p6He_11590_490.root", "sim_p6He_11590_500.root",

"sim_p6He_11600_220.root", "sim_p6He_11600_230.root", "sim_p6He_11600_240.root", "sim_p6He_11600_250.root", "sim_p6He_11600_260.root", "sim_p6He_11600_270.root", "sim_p6He_11600_280.root", "sim_p6He_11600_290.root", "sim_p6He_11600_300.root", "sim_p6He_11600_310.root", "sim_p6He_11600_320.root", "sim_p6He_11600_330.root", "sim_p6He_11600_340.root", "sim_p6He_11600_350.root", "sim_p6He_11600_360.root", "sim_p6He_11600_370.root", "sim_p6He_11600_380.root", "sim_p6He_11600_390.root", "sim_p6He_11600_400.root", "sim_p6He_11600_410.root", "sim_p6He_11600_420.root", "sim_p6He_11600_430.root", "sim_p6He_11600_440.root", "sim_p6He_11600_450.root", "sim_p6He_11600_460.root", "sim_p6He_11600_470.root", "sim_p6He_11600_480.root", "sim_p6He_11600_490.root", "sim_p6He_11600_500.root",

"sim_p6He_11610_220.root", "sim_p6He_11610_230.root", "sim_p6He_11610_240.root", "sim_p6He_11610_250.root", "sim_p6He_11610_260.root", "sim_p6He_11610_270.root", "sim_p6He_11610_280.root", "sim_p6He_11610_290.root", "sim_p6He_11610_300.root", "sim_p6He_11610_310.root", "sim_p6He_11610_320.root", "sim_p6He_11610_330.root", "sim_p6He_11610_340.root", "sim_p6He_11610_350.root", "sim_p6He_11610_360.root", "sim_p6He_11610_370.root", "sim_p6He_11610_380.root", "sim_p6He_11610_390.root", "sim_p6He_11610_400.root", "sim_p6He_11610_410.root", "sim_p6He_11610_420.root", "sim_p6He_11610_430.root", "sim_p6He_11610_440.root", "sim_p6He_11610_450.root", "sim_p6He_11610_460.root", "sim_p6He_11610_470.root", "sim_p6He_11610_480.root", "sim_p6He_11610_490.root", "sim_p6He_11610_500.root",

"sim_p6He_11620_220.root", "sim_p6He_11620_230.root", "sim_p6He_11620_240.root", "sim_p6He_11620_250.root", "sim_p6He_11620_260.root", "sim_p6He_11620_270.root", "sim_p6He_11620_280.root", "sim_p6He_11620_290.root", "sim_p6He_11620_300.root", "sim_p6He_11620_310.root", "sim_p6He_11620_320.root", "sim_p6He_11620_330.root", "sim_p6He_11620_340.root", "sim_p6He_11620_350.root", "sim_p6He_11620_360.root", "sim_p6He_11620_370.root", "sim_p6He_11620_380.root", "sim_p6He_11620_390.root", "sim_p6He_11620_400.root", "sim_p6He_11620_410.root", "sim_p6He_11620_420.root", "sim_p6He_11620_430.root", "sim_p6He_11620_440.root", "sim_p6He_11620_450.root", "sim_p6He_11620_460.root", "sim_p6He_11620_470.root", "sim_p6He_11620_480.root", "sim_p6He_11620_490.root", "sim_p6He_11620_500.root",

"sim_p6He_11630_220.root", "sim_p6He_11630_230.root", "sim_p6He_11630_240.root", "sim_p6He_11630_250.root", "sim_p6He_11630_260.root", "sim_p6He_11630_270.root", "sim_p6He_11630_280.root", "sim_p6He_11630_290.root", "sim_p6He_11630_300.root", "sim_p6He_11630_310.root", "sim_p6He_11630_320.root", "sim_p6He_11630_330.root", "sim_p6He_11630_340.root", "sim_p6He_11630_350.root", "sim_p6He_11630_360.root", "sim_p6He_11630_370.root", "sim_p6He_11630_380.root", "sim_p6He_11630_390.root", "sim_p6He_11630_400.root", "sim_p6He_11630_410.root", "sim_p6He_11630_420.root", "sim_p6He_11630_430.root", "sim_p6He_11630_440.root", "sim_p6He_11630_450.root", "sim_p6He_11630_460.root", "sim_p6He_11630_470.root", "sim_p6He_11630_480.root", "sim_p6He_11630_490.root", "sim_p6He_11630_500.root",

"sim_p6He_11640_220.root", "sim_p6He_11640_230.root", "sim_p6He_11640_240.root", "sim_p6He_11640_250.root", "sim_p6He_11640_260.root", "sim_p6He_11640_270.root", "sim_p6He_11640_280.root", "sim_p6He_11640_290.root", "sim_p6He_11640_300.root", "sim_p6He_11640_310.root", "sim_p6He_11640_320.root", "sim_p6He_11640_330.root", "sim_p6He_11640_340.root", "sim_p6He_11640_350.root", "sim_p6He_11640_360.root", "sim_p6He_11640_370.root", "sim_p6He_11640_380.root", "sim_p6He_11640_390.root", "sim_p6He_11640_400.root", "sim_p6He_11640_410.root", "sim_p6He_11640_420.root", "sim_p6He_11640_430.root", "sim_p6He_11640_440.root", "sim_p6He_11640_450.root", "sim_p6He_11640_460.root", "sim_p6He_11640_470.root", "sim_p6He_11640_480.root", "sim_p6He_11640_490.root", "sim_p6He_11640_500.root",

"sim_p6He_11650_220.root", "sim_p6He_11650_230.root", "sim_p6He_11650_240.root", "sim_p6He_11650_250.root", "sim_p6He_11650_260.root", "sim_p6He_11650_270.root", "sim_p6He_11650_280.root", "sim_p6He_11650_290.root", "sim_p6He_11650_300.root", "sim_p6He_11650_310.root", "sim_p6He_11650_320.root", "sim_p6He_11650_330.root", "sim_p6He_11650_340.root", "sim_p6He_11650_350.root", "sim_p6He_11650_360.root", "sim_p6He_11650_370.root", "sim_p6He_11650_380.root", "sim_p6He_11650_390.root", "sim_p6He_11650_400.root", "sim_p6He_11650_410.root", "sim_p6He_11650_420.root", "sim_p6He_11650_430.root", "sim_p6He_11650_440.root", "sim_p6He_11650_450.root", "sim_p6He_11650_460.root", "sim_p6He_11650_470.root", "sim_p6He_11650_480.root", "sim_p6He_11650_490.root", "sim_p6He_11650_500.root",

"sim_p6He_11660_220.root", "sim_p6He_11660_230.root", "sim_p6He_11660_240.root", "sim_p6He_11660_250.root", "sim_p6He_11660_260.root", "sim_p6He_11660_270.root", "sim_p6He_11660_280.root", "sim_p6He_11660_290.root", "sim_p6He_11660_300.root", "sim_p6He_11660_310.root", "sim_p6He_11660_320.root", "sim_p6He_11660_330.root", "sim_p6He_11660_340.root", "sim_p6He_11660_350.root", "sim_p6He_11660_360.root", "sim_p6He_11660_370.root", "sim_p6He_11660_380.root", "sim_p6He_11660_390.root", "sim_p6He_11660_400.root", "sim_p6He_11660_410.root", "sim_p6He_11660_420.root", "sim_p6He_11660_430.root", "sim_p6He_11660_440.root", "sim_p6He_11660_450.root", "sim_p6He_11660_460.root", "sim_p6He_11660_470.root", "sim_p6He_11660_480.root", "sim_p6He_11660_490.root", "sim_p6He_11660_500.root",

"sim_p6He_11670_220.root", "sim_p6He_11670_230.root", "sim_p6He_11670_240.root", "sim_p6He_11670_250.root", "sim_p6He_11670_260.root", "sim_p6He_11670_270.root", "sim_p6He_11670_280.root", "sim_p6He_11670_290.root", "sim_p6He_11670_300.root", "sim_p6He_11670_310.root", "sim_p6He_11670_320.root", "sim_p6He_11670_330.root", "sim_p6He_11670_340.root", "sim_p6He_11670_350.root", "sim_p6He_11670_360.root", "sim_p6He_11670_370.root", "sim_p6He_11670_380.root", "sim_p6He_11670_390.root", "sim_p6He_11670_400.root", "sim_p6He_11670_410.root", "sim_p6He_11670_420.root", "sim_p6He_11670_430.root", "sim_p6He_11670_440.root", "sim_p6He_11670_450.root", "sim_p6He_11670_460.root", "sim_p6He_11670_470.root", "sim_p6He_11670_480.root", "sim_p6He_11670_490.root", "sim_p6He_11670_500.root",

"sim_p6He_11680_220.root", "sim_p6He_11680_230.root", "sim_p6He_11680_240.root", "sim_p6He_11680_250.root", "sim_p6He_11680_260.root", "sim_p6He_11680_270.root", "sim_p6He_11680_280.root", "sim_p6He_11680_290.root", "sim_p6He_11680_300.root", "sim_p6He_11680_310.root", "sim_p6He_11680_320.root", "sim_p6He_11680_330.root", "sim_p6He_11680_340.root", "sim_p6He_11680_350.root", "sim_p6He_11680_360.root", "sim_p6He_11680_370.root", "sim_p6He_11680_380.root", "sim_p6He_11680_390.root", "sim_p6He_11680_400.root", "sim_p6He_11680_410.root", "sim_p6He_11680_420.root", "sim_p6He_11680_430.root", "sim_p6He_11680_440.root", "sim_p6He_11680_450.root", "sim_p6He_11680_460.root", "sim_p6He_11680_470.root", "sim_p6He_11680_480.root", "sim_p6He_11680_490.root", "sim_p6He_11680_500.root",

"sim_p6He_11690_220.root", "sim_p6He_11690_230.root", "sim_p6He_11690_240.root", "sim_p6He_11690_250.root", "sim_p6He_11690_260.root", "sim_p6He_11690_270.root", "sim_p6He_11690_280.root", "sim_p6He_11690_290.root", "sim_p6He_11690_300.root", "sim_p6He_11690_310.root", "sim_p6He_11690_320.root", "sim_p6He_11690_330.root", "sim_p6He_11690_340.root", "sim_p6He_11690_350.root", "sim_p6He_11690_360.root", "sim_p6He_11690_370.root", "sim_p6He_11690_380.root", "sim_p6He_11690_390.root", "sim_p6He_11690_400.root", "sim_p6He_11690_410.root", "sim_p6He_11690_420.root", "sim_p6He_11690_430.root", "sim_p6He_11690_440.root", "sim_p6He_11690_450.root", "sim_p6He_11690_460.root", "sim_p6He_11690_470.root", "sim_p6He_11690_480.root", "sim_p6He_11690_490.root", "sim_p6He_11690_500.root",

"sim_p6He_11700_220.root", "sim_p6He_11700_230.root", "sim_p6He_11700_240.root", "sim_p6He_11700_250.root", "sim_p6He_11700_260.root", "sim_p6He_11700_270.root", "sim_p6He_11700_280.root", "sim_p6He_11700_290.root", "sim_p6He_11700_300.root", "sim_p6He_11700_310.root", "sim_p6He_11700_320.root", "sim_p6He_11700_330.root", "sim_p6He_11700_340.root", "sim_p6He_11700_350.root", "sim_p6He_11700_360.root", "sim_p6He_11700_370.root", "sim_p6He_11700_380.root", "sim_p6He_11700_390.root", "sim_p6He_11700_400.root", "sim_p6He_11700_410.root", "sim_p6He_11700_420.root", "sim_p6He_11700_430.root", "sim_p6He_11700_440.root", "sim_p6He_11700_450.root", "sim_p6He_11700_460.root", "sim_p6He_11700_470.root", "sim_p6He_11700_480.root", "sim_p6He_11700_490.root", "sim_p6He_11700_500.root"

};


//"sim_p6He_11625_220.root", "sim_p6He_11625_230.root", "sim_p6He_11625_240.root", "sim_p6He_11625_250.root", "sim_p6He_11625_260.root", "sim_p6He_11625_270.root", "sim_p6He_11625_280.root", "sim_p6He_11625_290.root", "sim_p6He_11625_300.root", "sim_p6He_11625_310.root", "sim_p6He_11625_320.root", "sim_p6He_11625_330.root", "sim_p6He_11625_340.root", "sim_p6He_11625_350.root", "sim_p6He_11625_360.root", "sim_p6He_11625_370.root", "sim_p6He_11625_380.root", "sim_p6He_11625_390.root", "sim_p6He_11625_400.root", "sim_p6He_11625_410.root", "sim_p6He_11625_420.root", "sim_p6He_11625_430.root", "sim_p6He_11625_440.root", "sim_p6He_11625_450.root", "sim_p6He_11625_460.root", "sim_p6He_11625_470.root", "sim_p6He_11625_480.root", "sim_p6He_11625_490.root", "sim_p6He_11625_500.root",

//"sim_p6He_11635_220.root", "sim_p6He_11635_230.root", "sim_p6He_11635_240.root", "sim_p6He_11635_250.root", "sim_p6He_11635_260.root", "sim_p6He_11635_270.root", "sim_p6He_11635_280.root", "sim_p6He_11635_290.root", "sim_p6He_11635_300.root", "sim_p6He_11635_310.root", "sim_p6He_11635_320.root", "sim_p6He_11635_330.root", "sim_p6He_11635_340.root", "sim_p6He_11635_350.root", "sim_p6He_11635_360.root", "sim_p6He_11635_370.root", "sim_p6He_11635_380.root", "sim_p6He_11635_390.root", "sim_p6He_11635_400.root", "sim_p6He_11635_410.root", "sim_p6He_11635_420.root", "sim_p6He_11635_430.root", "sim_p6He_11635_440.root", "sim_p6He_11635_450.root", "sim_p6He_11635_460.root", "sim_p6He_11635_470.root", "sim_p6He_11635_480.root", "sim_p6He_11635_490.root", "sim_p6He_11635_500.root",

//"sim_p6He_11645_220.root", "sim_p6He_11645_230.root", "sim_p6He_11645_240.root", "sim_p6He_11645_250.root", "sim_p6He_11645_260.root", "sim_p6He_11645_270.root", "sim_p6He_11645_280.root", "sim_p6He_11645_290.root", "sim_p6He_11645_300.root", "sim_p6He_11645_310.root", "sim_p6He_11645_320.root", "sim_p6He_11645_330.root", "sim_p6He_11645_340.root", "sim_p6He_11645_350.root", "sim_p6He_11645_360.root", "sim_p6He_11645_370.root", "sim_p6He_11645_380.root", "sim_p6He_11645_390.root", "sim_p6He_11645_400.root", "sim_p6He_11645_410.root", "sim_p6He_11645_420.root", "sim_p6He_11645_430.root", "sim_p6He_11645_440.root", "sim_p6He_11645_450.root", "sim_p6He_11645_460.root", "sim_p6He_11645_470.root", "sim_p6He_11645_480.root", "sim_p6He_11645_490.root", "sim_p6He_11645_500.root",



//"sim_p6He_11685_250_v1.root", "sim_p6He_11685_260_v1.root", "sim_p6He_11685_270_v1.root", "sim_p6He_11685_280_v1.root", "sim_p6He_11685_290_v1.root", "sim_p6He_11685_300_v1.root", "sim_p6He_11685_305_v1.root", "sim_p6He_11685_310_v1.root", "sim_p6He_11685_315_v1.root",  "sim_p6He_11685_320_v1.root", "sim_p6He_11685_330_v1.root", "sim_p6He_11685_340_v1.root",  "sim_p6He_11685_350_v1.root", "sim_p6He_11685_360_v1.root", "sim_p6He_11685_370_v1.root", "sim_p6He_11685_380_v1.root", "sim_p6He_11685_390_v1.root", "sim_p6He_11685_400_v1.root",

//"sim_p6He_11690_250_v1.root", "sim_p6He_11690_260_v1.root", "sim_p6He_11690_270_v1.root", "sim_p6He_11690_280_v1.root", "sim_p6He_11690_290_v1.root", "sim_p6He_11690_300_v1.root", "sim_p6He_11690_305_v1.root", "sim_p6He_11690_310_v1.root", "sim_p6He_11690_315_v1.root",  "sim_p6He_11690_320_v1.root", "sim_p6He_11690_330_v1.root", "sim_p6He_11690_340_v1.root",  "sim_p6He_11690_350_v1.root", "sim_p6He_11690_360_v1.root", "sim_p6He_11690_370_v1.root", "sim_p6He_11690_380_v1.root", "sim_p6He_11690_390_v1.root", "sim_p6He_11690_400_v1.root"






//"sim_p6He_11695_320_v1.root", "sim_p6He_11690_320_v1.root", "sim_p6He_11685_320_v1.root", "sim_p6He_11680_320_v1.root", "sim_p6He_11675_320_v1.root", "sim_p6He_11670_320_v1.root", "sim_p6He_11665_320_v1.root", "sim_p6He_11660_320_v1.root", "sim_p6He_11655_320_v1.root", "sim_p6He_11650_320_v1.root", "sim_p6He_11645_320_v1.root", "sim_p6He_11640_320_v1.root", "sim_p6He_11635_320_v1.root", "sim_p6He_11630_320_v1.root", "sim_p6He_11625_320_v1.root", "sim_p6He_11620_320_v1.root", "sim_p6He_11615_320_v1.root", "sim_p6He_11610_320_v1.root", "sim_p6He_11605_320_v1.root", "sim_p6He_11600_320_v1.root", "sim_p6He_11595_320_v1.root", "sim_p6He_11590_320_v1.root", "sim_p6He_11585_320_v1.root", "sim_p6He_11580_320_v1.root", "sim_p6He_11575_320_v1.root", "sim_p6He_11570_320_v1.root", "sim_p6He_11565_320_v1.root", "sim_p6He_11560_320_v1.root", "sim_p6He_11555_320_v1.root", "sim_p6He_11550_320_v1.root", "sim_p6He_11545_320_v1.root", "sim_p6He_11540_320_v1.root", "sim_p6He_11535_320_v1.root", "sim_p6He_11530_320_v1.root", "sim_p6He_11525_320_v1.root", "sim_p6He_11520_320_v1.root", "sim_p6He_11515_320_v1.root", "sim_p6He_11510_320_v1.root", "sim_p6He_11505_320_v1.root", "sim_p6He_11500_320_v1.root"};

//"sim_p6He_11660_300_v1.root", "sim_p6He_11660_305_v1.root", "sim_p6He_11660_310_v1.root", "sim_p6He_11660_315_v1.root", "sim_p6He_11660_320_v1.root", "sim_p6He_11660_325_v1.root", "sim_p6He_11660_330_v1.root", "sim_p6He_11655_300_v1.root", "sim_p6He_11655_305_v1.root", "sim_p6He_11655_310_v1.root", "sim_p6He_11655_315_v1.root", "sim_p6He_11655_320_v1.root", "sim_p6He_11655_325_v1.root", "sim_p6He_11655_330_v1.root", "sim_p6He_11650_300_v1.root", "sim_p6He_11650_305_v1.root", "sim_p6He_11650_310_v1.root", "sim_p6He_11650_315_v1.root", "sim_p6He_11650_320_v1.root", "sim_p6He_11650_325_v1.root", "sim_p6He_11650_330_v1.root", "sim_p6He_11640_300_v1.root", "sim_p6He_11640_305_v1.root", "sim_p6He_11640_310_v1.root", "sim_p6He_11640_315_v1.root", "sim_p6He_11640_320_v1.root", "sim_p6He_11640_325_v1.root", "sim_p6He_11640_330_v1.root"};



//"sim_p6He_11695_180_v3.root", "sim_p6He_11695_185_v3.root", "sim_p6He_11695_190_v3.root", "sim_p6He_11695_195_v3.root", "sim_p6He_11695_200_v3.root", "sim_p6He_11695_205_v3.root", "sim_p6He_11695_210_v3.root", "sim_p6He_11695_215_v3.root", "sim_p6He_11695_220_v3.root", "sim_p6He_11695_225_v3.root", "sim_p6He_11695_230_v3.root", "sim_p6He_11695_235_v3.root", "sim_p6He_11695_240_v3.root", "sim_p6He_11695_245_v3.root", "sim_p6He_11695_250_v3.root", "sim_p6He_11695_255_v3.root", "sim_p6He_11695_260_v3.root", "sim_p6He_11695_265_v3.root", "sim_p6He_11695_270_v3.root", "sim_p6He_11695_275_v3.root", "sim_p6He_11695_280_v3.root", "sim_p6He_11695_285_v3.root", "sim_p6He_11695_290_v3.root", "sim_p6He_11695_295_v3.root", "sim_p6He_11695_300_v3.root", "sim_p6He_11695_305_v3.root", "sim_p6He_11695_310_v3.root", "sim_p6He_11695_315_v3.root", "sim_p6He_11695_320_v3.root", "sim_p6He_11695_325_v3.root", "sim_p6He_11695_330_v3.root", "sim_p6He_11695_335_v3.root", "sim_p6He_11695_340_v3.root", "sim_p6He_11695_345_v3.root", "sim_p6He_11695_350_v3.root", "sim_p6He_11695_355_v3.root", "sim_p6He_11695_360_v3.root", "sim_p6He_11695_365_v3.root", "sim_p6He_11695_370_v3.root", "sim_p6He_11695_375_v3.root", "sim_p6He_11695_380_v3.root", "sim_p6He_11695_385_v3.root", "sim_p6He_11695_390_v3.root", "sim_p6He_11695_395_v3.root", "sim_p6He_11695_400_v3.root", };


  std::ofstream chilog;
  chilog.open("record_grid3.txt", std::ios_base::app); // append instead of overwrite
  

  fit Fit(N,NN,Npeaks,para2,useit);

  double chisq_min;
  bool isStop = false;
  for (int j=0; j<n; j++)
  {
    for (int i=0; i<m;i++)
    {
      string filename1 = "/home/Li7star/li7sim/rootout/" + files1[j];
      string filename2 = "/home/Li7star/li7sim/rootout/" + files2[i];

      Fit.read2new(filename1, filename2);
      cout << Fit.ND << endl;
      
      chisq_min = Fit.powell(para,xi,ftol);
      cout << "chisq_min = " << chisq_min << endl;

      chilog << chisq_min << "\t" << "f1:" << files1[j] << "    f2:" << files2[i] << endl;

      //if( i == 6 && j == 7){ isStop = true; break;}

    }
    if (isStop) break;
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
