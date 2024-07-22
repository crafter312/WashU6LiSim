void plot_dedx(){

 
  
  ifstream File("loss/H_in_CH2.loss");
  
  if (File.is_open() != 1){
    cout << " could not open loss file " << endl;
    return;
  }

  char line[80];
  File.getline(line,80);
  //print the first line of the file
  cout << line << endl;
  
  int N;
  File >> N;
  cout << N << endl;
  
  float *Ein;
  float *dedx;
  Ein = new float [N];
  dedx = new float [N];

  for (int i=0;i<N;i++){
    File >> Ein[i] >> dedx[i];
    cout << Ein[i] << " " << dedx[i] << endl;
  }
 
  File.close();
  
  TCanvas *c1 = new TCanvas();
  c1->SetLogy();
  
  
  TGraph *gr = new TGraph(N, Ein, dedx);

  
  gr->GetXaxis()->SetTitle("Ein");
  gr->GetYaxis()->SetTitle("dedx");

  gr->Draw();
  
  
  
  delete []Ein;
  delete []dedx;
  
}
  

  
