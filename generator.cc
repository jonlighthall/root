/*-------------------Emacs PostScript "pretty-print" page width (97 columns)-------------------*/
/* Program: generator.cc
 *       Created  by Jon Lighthall
 * Purpose: 
 *       A set of macros simulating data
 */

#include "TMath.h"
#include "TRandom3.h"


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
  //source(1e6);
  // drpjxy("hbeam");
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
{//given z, theta, phi; calculate r
r=z_start/(TMath::Cos(theta_start*TMath::DegToRad()));
    
  if(doprint) {
    printf(" Emmission angle (theta,phi)=(%f,%f)\n",theta_start,phi_start);
    printf("  Current position is Z=%7.2f\n",Z);
    printf("  Radius is           r=%7.2f\n",r);
  }
 }

void trace_x(Double_t x_start, Double_t theta_start, Double_t phi_start)
{
X=r*(TMath::Sin(theta_start*TMath::DegToRad()))*(TMath::Cos(phi_start*TMath::DegToRad()));
if(doprint) {
printf("  X-position is %7.2f (relative), with offset %7.2f\n",X,x_start);
  }
  X+=x_start;  

if(doprint) {
printf("  X-position is %7.2f (absolute)\n",X);
  }
}

void trace_y(Double_t y_start, Double_t theta_start, Double_t phi_start)
{
    Y=r*(TMath::Sin(theta_start*TMath::DegToRad()))*(TMath::Sin(phi_start*TMath::DegToRad()));
  
  if(doprint) {
printf("  Y-position is %7.2f (relative), with offset %7.2f\n",Y,y_start);
  }
  Y+=y_start;  

if(doprint) {
printf("  Y-position is %7.2f (absolute)\n",Y);
 }
}

void source(Int_t nevents=1000, Double_t theta_min=0, Double_t theta_max=180, Double_t phi_min=0, Double_t phi_max=360)
{
  //beam spot-------------------------------------
  TRandom3 *rx=new TRandom3();//for x-position of beam spot
  TRandom3 *ry=new TRandom3();//for y-position of beam spot
  rx->SetSeed(0);
  ry->SetSeed(0);
 
  //polar angle (scattering angle)----------------
  TRandom3 *rtheta=new TRandom3();//for scattering angle 
  rtheta->SetSeed(0);
  if (gROOT->FindObject("htheta"))
    gROOT->FindObject("htheta")->Delete();    
  TH1D *htheta=new TH1D("htheta","Polar Angle",500,0,180);
 
  Double_t theta_center=30;
  //theta_min+=theta_center;
  //theta_max+=theta_center;
 
  //azimuthal angle-------------------------------
  TRandom3 *rphi=new TRandom3();
  rphi->SetSeed(0);
 
  if (gROOT->FindObject("hphi"))
    gROOT->FindObject("hphi")->Delete();    
  TH1D *hphi=new TH1D("hphi","Azimuth Angle",500,0,360);

  if (gROOT->FindObject("hangles"))
    gROOT->FindObject("hangles")->Delete();
  TH2D *hangles=new TH2D("hangles","Maks Plane",500,0,180,500,0,360);
  hangles->SetXTitle("theta - polar angle (deg)");
  hangles->SetYTitle("phi - azimuth angle (deg)");

  //X, Y positions (ray-tracing)------------------
  Bool_t hit=kFALSE;
   if (gROOT->FindObject("hmask"))
    gROOT->FindObject("hmask")->Delete();
   Double_t x_min=z_mask*(TMath::Tan(theta_min*TMath::DegToRad()));
   Double_t x_max=z_mask*(TMath::Tan(theta_max*TMath::DegToRad()));
   //x_min=0-118-6;
   //x_max=154-118+6;
   printf("xmin=%f, xmax=%f\n",x_min,x_max);
   //   TH2D *hmask=new TH2D("hmask","Mask Plane",500,x_max,x_max,500,x_min,x_max);
   TH2D *hmask=new TH2D("hmask","Mask Plane",500,-80,-80,500,-80,80);
   if (gROOT->FindObject("hmaskg"))
     gROOT->FindObject("hmaskg")->Delete();
   TH2D *hmaskg=new TH2D("hmaskg","Mask Plane (gated)",500,-80,-80,500,-80,80);

   //Other correlation plots----------------------
   TH2D *hxtheta=new TH2D("hxtheta","Theta (polar) vs. X",100,-80,-80,100,0,180);
   hxtheta->SetYTitle("theta - polar angle (deg)");
   TH2D *hytheta=new TH2D("hytheta","Theta (polar) vs. Y",100,-80,-80,100,0,180);
   hytheta->SetYTitle("theta - polar angle (deg)");
   TH2D *hxphi=new TH2D("hxphi","Phi (azimuthal) vs. X",100,-80,-80,500,0,360);
   hxphi->SetYTitle("phi - azimuth angle (deg)");
   TH2D *hyphi=new TH2D("hyphi","Phi (azimuthal) vs. Y",100,-80,-80,500,0,360);
   hyphi->SetYTitle("phi - azimuth angle (deg)");

   for (Int_t i=0; i<nevents; i++){
    hit=kFALSE;
    //Origin position (beam spot)-----------------
    // x=rx->Gaus(offset_x,sigma_x);
    //y=rx->Gaus(offset_y,sigma_y);
    x=rx->Uniform(-sigma_x,sigma_x);
    y=ry->Uniform(-sigma_y,sigma_y);
    hbeam->Fill(x,y);
    
    //Emmission angle-----------------------------
    // Polar angle----------------------
    //theta=rtheta->Uniform(theta_min,theta_max);
    theta=(TMath::ACos(rtheta->Uniform(-1,1))*(TMath::RadToDeg()));
    htheta->Fill(theta);
    // Azimuthal angle------------------
    phi=rphi->Uniform(phi_min,phi_max);
    hphi->Fill(phi);
    hangles->Fill(theta,phi);
    
    if(i%5000==0) {
      printf("%5.1f\%: %d events generated\n",(double)i/nevents*100,i);
      doprint=kTRUE;
    }
    else
      doprint=kFALSE;
    
    //calculate positions at mask
    //theta-=theta_center;
    Z=z_mask;
    trace_r(Z,theta,phi);
    trace_x(x,theta,phi);
    trace_y(y,theta,phi);
    rho=TMath::Sqrt((TMath::Power(X,2))+(TMath::Power(Y,2)));
    if(doprint)
      printf("  rho=%f, z=%f\n",rho,TMath::Sqrt(r*r-rho*rho));
    if(rho>(99.5/2)) {
      hit=kTRUE;
      if(doprint)
	printf("  Hit!\n");
    }
    //if(((X>-31)&&(X<-26))||((X>19)&&(X<24)))
    //hit=kTRUE;
    // if(X>0)
    //hit=kTRUE;
    hxtheta->Fill(X,theta);
    hytheta->Fill(Y,theta);
    hxphi->Fill(X,phi);
    hyphi->Fill(Y,phi);



    hmask->Fill(X,Y);
    if(!hit)
      hmaskg->Fill(X,Y);
}
}
