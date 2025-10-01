//Macro which accepts SRIM output loss files and returns the correct format
//For use in the WashU loss functions.
//This code currently only rewrites loss files for elements with stable isotopes

//The output file name is (projectile Z)_(Target Z(s)).loss
//So for Be in H and C (regardless of stoichiometry), the output is Be_HC.loss
//This stoichiometry could be programmed in or renamed after conversion.
//For some stoichiometries, you might want a special name (ex. BC400)
//I will leave that up to the user.

//Run in a terminal with ROOT.
//Easiest to do <root -q 'SRIMtoLoss.C("<SRIM File.txt>")'>
//Currently pulls from the "SRIMFiles" directory, could be programmed differently.
//I think it's nice to separate the LossFile and SRIMFiles directories.

void SRIMtoLoss(string fname)
{

  //Want to build a library of stable nuclei, using the most naturally abundant isotope for E/A
  //Read in wallet card data
  string filename = "NNDC_WalletCards2023.txt";
  ifstream file(filename.c_str());
  if (file.is_open() != 1)
  {
    cout << " could not open file " << filename << endl;
    abort();
  }

  //Data format
  // Z(int)  A(int)  X(string)  HL(string)  HL(unit)  Abundance(%)  ME(float)  ME(unit)  

  //Sometimes the HL will be read in as a range (#%-#%).
  //Read HL in as a string and take lowest value in range and convert to float.

  //Also skip the 1 line header
  file.ignore(200,'\n');

  string ignore1; //HL unit

  vector<int> Z;
  vector<int> A;
  vector<string> name;
  vector<int> HL;
  vector<string> abundance;
  vector<float> abunf; //abundance as a float
  int Z_temp = 0;
  int A_temp = 0;
  string HL_temp;
  string name_temp;
  string abun_temp;

  //For saving stable nuclei and the highest elemental abundance
  //118 elements, but want to match Z to index
  int stabZ[119] = {0};
  int stabA[119] = {0};
  string stabname[119] = {"ignore"};
  float stababun[119] = {0};

  //Strings can be tab delimited and contain whitespace, use getLine
  string line;
  while (getline(file, line)) {

    istringstream ss(line);
    ss >> Z_temp;
    ss >> A_temp;
    ss >> name_temp;
    ss >> HL_temp;
    ss >> ignore1;

    int counter = 0;
    //Get string with whitespace
    while (getline(ss, abun_temp, '\t')) {
      //if (Z_temp == 1 && A_temp == 1) std::cout << "Extracted segment: [" << abun_temp << "]" << std::endl;
      counter++;
      if (counter == 2) break; //Only take the second
    }

    file.ignore(100,'\n');

    //Break at E limit or end of file
    if (file.eof() || file.bad()) break;

    //SMEx.push_back(SMEx_temp);
    Z.push_back(Z_temp);
    A.push_back(A_temp);
    name.push_back(name_temp);
    //cout << "HL temp " << HL_temp << endl;
    if (HL_temp == "STABLE") {
      HL.push_back(-1); //Assign -1 for stable nucleus;
    }
    else HL.push_back(0); //Assign everything else zero

    //Check if string is == -1
    if (abun_temp == "-1") {
      abunf.push_back(stof(abun_temp));
    }
    else {
      //Look at abundance and take lowest value if there is a range
      char targetchar = '-';
      size_t charpos = abun_temp.find(targetchar);

      //Check if it found the target character
      if (charpos != string::npos) {
        string beforechar = abun_temp.substr(0,charpos);
        abundance.push_back(beforechar);
        //cout << abun_temp << " " << beforechar << endl;
        abunf.push_back(stof(beforechar));
      }
      else {
        abunf.push_back(stof(abun_temp));
      }
    }
  }

  //Loop over vectors and find most abundant stable isotopes
  //Add to arrays
  for (int i=0;i<Z.size();i++) {
    if (HL[i] == -1) {
      if (abunf[i] > stababun[Z[i]]) {
        stabZ[Z[i]] = Z[i];
        stabA[Z[i]] = A[i];
        stabname[Z[i]] = name[i];
        stababun[Z[i]] = abunf[i];

        //cout << "here" << endl;
      }
    }
    else continue;
  }

  //for (int i=0;i<119;i++) cout << stabZ[i] << " " << stabA[i] << " " << stabname[i] << " " << stababun[i] << endl;

  //Read in SRIM file
  string sname = "SRIMFiles/" + fname;
  ifstream sfile(sname.c_str());
  if (sfile.is_open() != 1)
  {
    cout << " could not open file " << sname << endl;
    abort();
  }

  //Need to know nucleus type to divide by the right A.
  string nucleus;
  float nuclA;

  //SRIM files start with 24 lines of header info, skip the first 7
  for (int i=0;i<7;i++) sfile.ignore(100,'\n');

  string projline;
  string sprojZ;
  int projZ;
  bool foundpZ = false;
  //Read in projectile info and take Z. Look for [Z]
  for (int i=0;i<10;i++) {
    sfile >> projline;

    //Look for brackets
    size_t starbrak = projline.find('[');
    size_t endbrak = projline.find(']');
    if (starbrak == string::npos || endbrak == string::npos) continue;
    else foundpZ = true;

    if (foundpZ == true) {
      sprojZ = projline.substr(starbrak + 1, endbrak - starbrak - 1);
      projZ = stoi(sprojZ);
      sfile.ignore(100,'\n');
      break;
    }
  }

  if (foundpZ == false) {
    cout << "Could not find projectile Z, abort" << endl;
    abort();
  }
  else { 
    nuclA = stabA[projZ];
    cout << "Projectile Z: " << projZ << " A: " << nuclA << endl;
  }

  //Skip 6 lines before reading in target
  for (int i=0;i<6;i++) sfile.ignore(100,'\n');

  //Read in target. Can be any number of elements. End row search once you reach ======...
  string targline;  
  string stargZ;
  vector <int> targZ;
  vector <float> targAP;
  for (;;) {
  
  //Get line and break if you reach "====..."
  string line1;
  sfile >> line1;
  char item1 = line1[0];

  if (item1 == '=') {
    sfile.ignore(100,'\n');
    break;
  }   

    for (int i =1;i<3;i++) {
      //Easy search. Always goes "X" "Z" "%" etc
      sfile >> stargZ;
      if (i == 1) {
        targZ.push_back(stoi(stargZ));
      }
    
      if (i == 2) {
        targAP.push_back(stof(stargZ));
      }
    }

    sfile.ignore(100,'\n');

  }

  if (targZ.size() == 0) {
    cout << "No target" << endl;
    abort();
  }

  for (int i=0;i<targZ.size();i++) {
    cout << "Target Z: " << targZ[i] << " atomic %: " << targAP[i] << endl;
  }

  // Skip remaining header
  for (int i=0;i<7;i++) sfile.ignore(100,'\n');

  //Read in incident energy, energy unit, and Coulomb energy loss, ignore other stuff
  vector<float> energies;
  vector<float> EPA;
  vector<string> eunit;
  vector<float> Celoss;

  string tempe;
  string tempunit;
  string temploss;
  int slinecount = 0;

  //Read in values

  for (;;) {

    sfile >> tempe >> tempunit >> temploss;

    if (tempe[0] == '-') break; //after reaching bottom (----------...)

    eunit.push_back(tempunit);
    energies.push_back(stof(tempe));
    Celoss.push_back(stof(temploss));

    //Want units of MeV
    if (eunit[slinecount] == "eV") energies[slinecount] = energies[slinecount]/1000000.;
    else if (eunit[slinecount] == "keV") energies[slinecount] = energies[slinecount]/1000.;
    else if (eunit[slinecount] == "GeV") energies[slinecount] = energies[slinecount]*1000.;

    //Want in EPA
    EPA.push_back(energies[slinecount]/nuclA);

    //Write out to loss file
    sfile.ignore(100,'\n');
    slinecount++;

  }

  //Now write out loss file in the following format
  //Header with "energy loss of <Zproj> in <Ztarg> E/A MeV/mg/cm2 
  //Number of E and Eloss pairs
  //Rows as <E> <Eloss>
  
  //Filename should be <Zproj>_Z<targ>.loss
  string outname;
  string targtotal;
  for (int i=0;i<targZ.size();i++) {
    targtotal += stabname[targZ[i]];
  }
  outname = stabname[projZ] + "_" + targtotal + ".loss";  

  cout << "Writing to SRIMFiles/" << outname << endl;

  ofstream fileout;
  fileout.open("./SRIMFiles/" + outname);

  if (!fileout.is_open()) {
    cout << "Error: Could not open the file 'output.txt' for writing." << endl;
    abort(); // Indicate an error
  }

  fileout << "energy loss of " << stabname[projZ] << " in " << targtotal << " E/A Mev/mg/cm2" << endl;
  fileout << EPA.size() << endl;

  for (int i=0;i<EPA.size();i++) {
    fileout << EPA[i] << "\t" << Celoss[i] << endl;
  }

  //Close files
  file.close();
  sfile.close();
  fileout.close();
}
