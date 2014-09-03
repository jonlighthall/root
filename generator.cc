/*-------------------Emacs PostScript "pretty-print" page width (97 columns)-------------------*/
/* Program: generator.cc
 *       Created  by Jon Lighthall
 * Purpose: 
 *       A set of macros simulating data
 */

#include "TMath.h"

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


Double_t sigma_x=0.607956845;//for 90% in a 1mm radius
Double_t sigma_y=sigma_x;
Double_t offset_x=0, offset_y=0;
//TH2D *hbeam=0;
//TH2D *htheta=0;

void createbeam(Double_t set_offset_x=0, Double_t set_offset_y=0, Double_t set_sigma_x=0.607956845, Double_t set_sigma_y=0.607956845)
{
  sigma_x=set_sigma_x;
  sigma_y=set_sigma_y;
  offset_x=set_offset_x;
  offset_y=set_offset_y;

  if (gROOT->FindObject("hbeam"))
    gROOT->FindObject("hbeam")->Delete();   
  TH2D *hbeam=new TH2D("hbeam","Beam Spot",500,-10,10,500,-10,10);
  hbeam->SetXTitle("x-position (mm)");
  hbeam->SetYTitle("y-position (mm)");
  source(1e6);
  drpjxy("hbeam");
}

Float_t z_mask=637.0;
Float_t z_window=z_mask+3.25;
Float_t z_A2=z_window+6.35+12.7+23.6;
Float_t z_A1=z_A2+36.8;//anode separation
Float_t delta_z=0.125*25.4;//plane separation

Float_t z_Y2=z_A2-delta_z;
Float_t z_X2=z_A2+delta_z;
Float_t z_Y1=z_A1-delta_z;
Float_t z_X1=z_A1+delta_z;


Double_t x=0,y=0;//point on target
Double_t theta=0;//scattering angle
Double_t phi=0;//azimuthal angle
Double_t X=0, Y=0, Z=0;//position
Double_t r=0,rho=0;//radius
Bool_t doprint=kFALSE;


void trace_r(Double_t z_start,Double_t theta_start, Double_t phi_start)
{//given z, theta, phi; calcuate r
  r=z_start/(TMath::Cos(theta_start*TMath::DegToRad()));
  
  
  if(doprint) {
    printf(" Emmission angle (theta,phi)=(%f,%f)\n",theta_start,phi_start);
    printf("  Current position is Z=%f\n",Z);
    printf("  Radius is           r=%f\n",r);
  }
 }

void trace_x(Double_t x_start, Double_t theta_start, Double_t phi_start)
{
X=r*(TMath::Sin(theta_start*TMath::DegToRad()))*(TMath::Cos(phi_start*TMath::DegToRad()));
  
  if(doprint) {
printf("  X-position is %f (relative), with offset %f\n",X,x);
  }
  X+=x;  

if(doprint) {
printf("  X-position is %f (absolute)\n",X);
  }
}

void trace_y(Double_t y_start, Double_t theta_start, Double_t phi_start)
{
    Y=r*(TMath::Sin(theta_start*TMath::DegToRad()))*(TMath::Sin(phi_start*TMath::DegToRad()));
  
  if(doprint) {
printf("  Y-position is %f (relative), with offset %f\n",Y,y);
  }
  Y+=y;  

if(doprint) {
printf("  Y-position is %f (absolute)\n",Y);
 }
}




void source(Int_t nevents=1000)
{
  //TObjString histnames[20]={"hbeam","htheta"}
  //beam spot
  TRandom *rx=new TRandom();//for x-position of beam spot
  TRandom *ry=new TRandom();//for y-position of beam spot
 
 
  //scattering angle
  TRandom *rtheta=new TRandom();//for scattering angle 
  if (gROOT->FindObject("htheta"))
    gROOT->FindObject("htheta")->Delete();    
  TH1D *htheta=new TH1D("htheta","Emission Angle",500,0,180);
  
  Double_t theta_min=7;
  Double_t theta_max=-theta_min;
  Double_t theta_center=30;
  theta_min+=theta_center;
  theta_max+=theta_center;
 
  //azimuthal angle
  TRandom *rphi=new TRandom();
  
  Double_t phi_min=0;
  Double_t phi_max=360;
  if (gROOT->FindObject("hphi"))
    gROOT->FindObject("hphi")->Delete();    
  TH1D *hphi=new TH1D("hphi","Emission Angle",500,0,180);

  //Ray
  
  Bool_t hit=kFALSE;

  for (Int_t i=0; i<nevents; i++) {
    x=rx->Gaus(offset_x,sigma_x);
    y=rx->Gaus(offset_y,sigma_y);
    hbeam->Fill(x,y);
    theta=rtheta->Uniform(theta_min,theta_max);
    htheta->Fill(theta);
    phi=rphi->Uniform(phi_min,phi_max);
    hphi->Fill(phi);
    
    if(i%500000==0){
      printf("%5.1f\%: %d events generated\n",(double)i/nevents*100,i);
      doprint=kTRUE;
    }
    else
      doprint=kFALSE;
    
    Z=z_mask;
    trace_r(Z,theta,phi);
    trace_x(x,theta,phi);
    trace_y(y,theta,phi);
    rho=TMath::Sqrt((TMath::Power(X,2))+(TMath::Power(X,2)));
    if(doprint)
      printf("  rho=%f\n",rho);
    
   
}
}
