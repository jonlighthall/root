/*-------------------Emacs PostScript "pretty-print" page width (97 columns)-------------------*/
/* Program: generator.cc
 *       Created  by Jon Lighthall
 * Purpose: 
 *       A set of macros simulating data
 */

void testRandom(Int_t nrEvents=500000000)
{
  TRandom *r1=new TRandom();
  TRandom2 *r2=new TRandom2();
  TRandom3 *r3=new TRandom3();
  TH1D *h1=new TH1D("h1","TRandom",500,0,1);
  TH1D *h2=new TH1D("h2","TRandom2",500,0,1);
  TH1D *h3=new TH1D("h3","TRandom3",500,0,1);
  TStopwatch *st=new TStopwatch();
  st->Start();
  for (Int_t i=0; i<nrEvents; i++) { h1->Fill(r1->Uniform(0,1)); }
  st->Stop(); cout << "Random: " << st->CpuTime() << endl;
  st->Start();
  for (Int_t i=0; i<nrEvents; i++) { h2->Fill(r2->Uniform(0,1)); }
  st->Stop(); cout << "Random2: " << st->CpuTime() << endl;
  st->Start();
  for (Int_t i=0; i<nrEvents; i++) { h3->Fill(r3->Uniform(0,1)); }
  st->Stop(); cout << “Random3: “ << st->CpuTime() << endl;
}

void testRandom2(Int_t nrEvents=10e+08, double mean = 0, double sigma = 100)
{
  TRandom *r1=new TRandom();
  TRandom *r2=new TRandom();
  TRandom *r3=new TRandom();
  
  if ((TH2F *) gROOT->FindObject("h1")) {
    gROOT->FindObject("h1")->Delete();  
    printf("Histogram \"h1\" already exists. Deleting old histogram.\n");
  } 
  if ((TH2F *) gROOT->FindObject("h2")) {
    gROOT->FindObject("h2")->Delete();  
    printf("Histogram \"h2\" already exists. Deleting old histogram.\n");
  } 
  if ((TH2F *) gROOT->FindObject("h3")) {
    gROOT->FindObject("h3")->Delete();  
    printf("Histogram \"h3\" already exists. Deleting old histogram.\n");
  } 
  
  TH1D *h1=new TH1D("h1","TRandom, Uniform",500,mean-sigma*4,mean+sigma*4);
  TH1D *h2=new TH1D("h2","TRandom, Gaus",500,mean-sigma*4,mean+sigma*4);
  TH1D *h3=new TH1D("h3","TRandom, Rndm",5000,0,1);
  
  for (Int_t i=0; i<nrEvents; i++) {
    h1->Fill(r1->Uniform(mean-sigma*4,mean+sigma*4));
    h2->Fill(r2->Gaus(mean,sigma));
    h3->Fill(r3->Rndm());
    if(i%500000==0)
      printf("%5.1f\%: %d events generated\n",(double)i/nrEvents*100,i);
  }
  
  plotall("h");
  cFit->cd(2);
  gfit("h2");
}

void source(Int_t nevents=1000)
{
  TRandom *r1=new TRandom();//for angle 
  
  TH1D *htheta=new TH1D("htheta","TRandom, Uniform",500,mean-sigma*4,mean+sigma*4);
  
  for (Int_t i=0; i<nrEvents; i++) {
    h1->Fill(r1->Uniform(mean-sigma*4,mean+sigma*4));
    h2->Fill(r2->Gaus(mean,sigma));
    h3->Fill(r3->Rndm());
    if(i%500000==0)
      printf("%5.1f\%: %d events generated\n",(double)i/nrEvents*100,i);
  }
  
  plotall("h");
  cFit->cd(2);
  gfit("h2");
}
