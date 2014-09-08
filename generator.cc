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
Bool_t is_gaussian=kTRUE;
TH2D *hbeam;

void setbeam(Double_t set_offset_x=0, Double_t set_offset_y=0, Double_t set_sigma_x=0.607956845, Double_t set_sigma_y=0.607956845, Bool_t set_is_gaussian=kTRUE)
{
  sigma_x=set_sigma_x;
  sigma_y=set_sigma_y;
  offset_x=set_offset_x;
  offset_y=set_offset_y;
  is_gaussian=set_is_gaussian;
  TString distrib;
  if(is_gaussian)
    distrib="Gaussian";
  else
    distrib="Uniform";
  printf("Beam created at x=%5.2f, y=%5.2f\n with sigma_x=%5.2f, sigma_y=%5.2f;\n distribution is %s.\n",offset_x,offset_y,sigma_x,sigma_y,distrib.Data());

  clearhists();  
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
Bool_t doprint=1;//kFALSE;
Bool_t diag=kFALSE;

void trace_r(Double_t z_start,Double_t theta_start, Double_t phi_start)
{//given z, theta, phi; calculate r
  r=z_start/(TMath::Cos(theta_start*TMath::DegToRad()));
    
  if(iprint) {
    printf(" Emmission angle (theta,phi)=(%f,%f)\n",theta_start,phi_start);
    printf("  Current position is Z=%7.2f\n",Z);
    printf("  Radius is           r=%7.2f\n",r);
  }
}

void trace_x(Double_t x_start, Double_t theta_start, Double_t phi_start)
{
  X=r*(TMath::Sin(theta_start*TMath::DegToRad()))*(TMath::Cos(phi_start*TMath::DegToRad()));
  if(iprint) {
    printf("  X-position is %7.2f (relative), with offset %7.2f\n",X,x_start);
  }
  X+=x_start;  

  if(iprint) {
    printf("  X-position is %7.2f (absolute)\n",X);
  }
}

void trace_y(Double_t y_start, Double_t theta_start, Double_t phi_start)
{
  Y=r*(TMath::Sin(theta_start*TMath::DegToRad()))*(TMath::Sin(phi_start*TMath::DegToRad()));
  
  if(iprint) {
    printf("  Y-position is %7.2f (relative), with offset %7.2f\n",Y,y_start);
  }
  Y+=y_start;  

  if(iprint) {
    printf("  Y-position is %7.2f (absolute)\n",Y);
  }
}

Bool_t hit_mask(void)
{
  if(rho>(99.5/2)) 
    return kTRUE;//hits radius
  if((X>-31)&&(X<-26))
    return kTRUE;//hist left strip
  if((X>19)&&(X<24))
    return kTRUE;//hits right strip
  if((Y>-15)&&(Y<-10))
    return kTRUE;//hits bottom strip
  if((Y>10)&&(Y<15))
    return kTRUE;//hits top strip
  return kFALSE;
}

Bool_t hit_shield(void)
{
  if((X<-118)||(X>36))
    return kTRUE;//hits shield (x)
  if((Y<-27)||(Y>27))
    return kTRUE;//hist shield (y)
  return kFALSE;
}

Bool_t hit_window(void)
{
  if(X>34.0)
    return kTRUE;
  return kFALSE;
}

Double_t theta_min=0;
Double_t theta_max=180;
Double_t phi_min=0;
Double_t phi_max=360;

TH1D *htheta;
TH1D *hphi;
TH2D *hangles;

TH2D *hxtheta;
TH2D *hytheta;
TH2D *hxphi;
TH2D *hyphi;

TH2D *hmask;
TH2F *hmaskg;
TH1F *hhit;
TH2F *hwin;
TH2F *hwing;

TH1F *hx[4];
TH1F *hxg[4];
TH2F *hyx[2];
TH2F *hyxg[2];
TH2F *hyxgm[2];

void setangles(Double_t set_theta_min=0, Double_t set_theta_max=180, Double_t set_phi_min=0, Double_t set_phi_max=360)
{//define angles and build histograms
  theta_min=set_theta_min;
  theta_max=set_theta_max;
  phi_min=set_phi_min;
  phi_max=set_phi_max;
  printf("Emmission angles defined over:\n Theta %f to %f\n Cos(theta) %f to %f\n Phi %f to %f\n",theta_min,theta_max,TMath::Cos(theta_min*TMath::DegToRad()),TMath::Cos(theta_max*TMath::DegToRad()),phi_min,phi_max);
  clearhists();
}

Double_t x_min=0;
Double_t x_max=100;

void definehists()
{
  printf("Defining histograms...\n");
  hbeam=new TH2D("hbeam","Beam Spot",500,-10,10,500,-10,10);
  hbeam->SetXTitle("x-position (mm)");
  hbeam->SetYTitle("y-position (mm)"); 

  htheta=new TH1D("htheta","Theta - Polar Angle",500,0,180);
  hphi=new TH1D("hphi","Phi - Azimuthal Angle",500,0,360);
  hangles=new TH2D("hangles","Phi vs. Theta",500,0,180,500,0,360);
  hangles->SetXTitle("Theta - Polar Angle (deg)");
  hangles->SetYTitle("Phi - Azimuthal Angle (deg)");
  
  x_min=z_mask*(TMath::Tan(theta_min*TMath::DegToRad()));
  x_max=z_mask*(TMath::Tan(theta_max*TMath::DegToRad()));
  Double_t x_Max=124;
  if((x_max>x_Max)||(x_max<0))
    x_max=x_Max;
  x_min=-x_max;
  
  if(doprint)
    printf(" Mask histogram limits are xmin=%f, xmax=%f\n",x_min,x_max);
  hmask=new TH2D("hmask" ,"Mask Plane",500,x_min,x_max,500,x_min,x_max);
  hmaskg=new TH2F("hmaskg","Mask Plane (gated)",500,x_min,x_max,500,x_min,x_max);
  hhit=new TH1F("hhit","hhit",7,-1,8);
  hwin=new TH2F("hwin" ,"Window Plane",500,x_min,x_max,500,x_min,x_max);
  hwing=new TH2F("hwing","Window Plane (gated)",500,x_min,x_max,500,x_min,x_max);

  if(diag) {
    //Other correlation plots----------------------
    hxtheta=new TH2D("hxtheta","Theta (polar) vs. X",100,x_min,x_max,100,0,180);
    hxtheta->SetYTitle("theta - polar angle (deg)");
    hytheta=new TH2D("hytheta","Theta (polar) vs. Y",100,x_min,x_max,100,0,180);
    hytheta->SetYTitle("theta - polar angle (deg)");
    hxphi=new TH2D("hxphi","Phi (azimuthal) vs. X",100,x_min,x_max,500,0,360);
    hxphi->SetYTitle("phi - azimuth angle (deg)");
    hyphi=new TH2D("hyphi","Phi (azimuthal) vs. Y",100,x_min,x_max,500,0,360);
    hyphi->SetYTitle("phi - azimuth angle (deg)");
  }

  hx[0] = new TH1F("hx0","X1 Position",500,x_min,x_max);
  hx[1] = new TH1F("hx1","Y1 Position",500,x_min,x_max);
  hx[2] = new TH1F("hx2","X2 Position",500,x_min,x_max);
  hx[3] = new TH1F("hx3","Y2 Position",500,x_min,x_max);
  hxg[0] = new TH1F("hxg0","X1 Position (gated)",500,x_min,x_max);
  hxg[1] = new TH1F("hxg1","Y1 Position (gated)",500,x_min,x_max);
  hxg[2] = new TH1F("hxg2","X2 Position (gated)",500,x_min,x_max);
  hxg[3] = new TH1F("hxg3","Y2 Position (gated)",500,x_min,x_max);
  /* for(int i = 0; i < 4; i++) {//1D position plots
     hx[i]->SetXTitle("Relative Position");
     hx[i]->SetYTitle("Number of Entries");
     }*/
  hyx[0] = new TH2F("hyx0","Y vs X Positions",500,x_min,x_max,500,x_min,x_max);
  hyx[1] = new TH2F("hyx1","Y vs X Positions",500,x_min,x_max,500,x_min,x_max);
  hyxg[0] = new TH2F("hyxg0","Y vs X Positions (gated)",500,x_min,x_max,500,x_min,x_max);
  hyxg[1] = new TH2F("hyxg1","Y vs X Positions (gated)",500,x_min,x_max,500,x_min,x_max);
  hyxgm[0] = new TH2F("hyxgm0","Y vs X Positions (gated), measured",500,x_min,x_max,500,x_min,x_max);
  hyxgm[1] = new TH2F("hyxgm1","Y vs X Positions (gated), measured",500,x_min,x_max,500,x_min,x_max);
}

void clearhists()
{
  printf("Clearing histograms...\n");
  const int Nhists = 31;
  TString histnames[Nhists]={"hbeam","htheta","hphi","hangles","hmask","hmaskg","hxtheta","hytheta","hxphi","hyphi","hhit","hx0","hx1","hx2","hx3","hwin","hwing","hx0","hx1","hx2","hx3","hxg0","hxg1","hxg2","hxg3","hyx0","hyx1","hyxg0","hyxg1","hyxgm0","hyxgm1"};

  for(int i=0; i < Nhists; i++) {
    if (gROOT->FindObject(histnames[i])) {
      if(doprint)
	printf(" Histogram %10s already exits.  Clearing...\n",histnames[i].Data());
      gROOT->FindObject(histnames[i])->Clear();    
    }
  }
}

Bool_t iprint=kFALSE; //doprint;  
void source(Int_t nevents=1000, Bool_t set_doprint=kFALSE)
{
  if(!(gROOT->FindObject("hbeam")))
    definehists();
  printf("Generating %.0f events...\n",nevents);

  doprint=set_doprint;
  //beam spot-------------------------------------
  TRandom3 *rx=new TRandom3();//for x-position of beam spot
  TRandom3 *ry=new TRandom3();//for y-position of beam spot
  rx->SetSeed(0);
  ry->SetSeed(0);
 
  //polar angle-----------------------------------
  TRandom3 *rtheta=new TRandom3();//for scattering angle 
  rtheta->SetSeed(0);
  Double_t cos_theta_min=(TMath::Cos(theta_min*TMath::DegToRad()));
  Double_t cos_theta_max=(TMath::Cos(theta_max*TMath::DegToRad()));
  Double_t theta_center=30;
  //theta_min+=theta_center;
  //theta_max+=theta_center;
 
  //azimuthal angle-------------------------------
  TRandom3 *rphi=new TRandom3();
  rphi->SetSeed(0);
 
  //X, Y positions (ray-tracing)------------------
  Bool_t hit=kFALSE;
  Bool_t miss=kTRUE;
   
  int step=(int) (nevents/10);
  int step_max=5e4;
  if(step>step_max)
    step=step_max;  
    
  for (Int_t i=0; i<nevents; i++) {
    if(i%step==0) {
      printf("%5.1f%%: %d events generated\n",(double)i/nevents*100,i);
      iprint=doprint*kTRUE;
    }
    else
      iprint=kFALSE;  
    hit=kFALSE;//hit means an obstruction has been hit
    miss=kTRUE;//miss means the trajectory is unobstructed
    //Position of origin (beam spot)-----------------
    if(is_gaussian) {
      x=rx->Gaus(offset_x,sigma_x);
      y=rx->Gaus(offset_y,sigma_y);
    }
    else {
      x=rx->Uniform(-sigma_x,sigma_x);
      y=ry->Uniform(-sigma_y,sigma_y);
    }
    hbeam->Fill(x,y);
    
    //Emmission angle-----------------------------
    // Polar angle----------------------
    theta=(TMath::ACos(rtheta->Uniform(cos_theta_min,cos_theta_max))*(TMath::RadToDeg()));
    htheta->Fill(theta);
    // Azimuthal angle------------------
    phi=rphi->Uniform(phi_min,phi_max);
    hphi->Fill(phi);
    hangles->Fill(theta,phi);
     
    //calculate positions at mask-------
    Z=z_mask;
    trace_r(Z,theta,phi);
    trace_x(x,theta,phi);
    trace_y(y,theta,phi);
    rho=TMath::Sqrt((TMath::Power(X,2))+(TMath::Power(Y,2)));
    if(iprint)
      printf("  rho=%f, z=%f\n",rho,TMath::Sqrt(r*r-rho*rho));    

    //diagnostics-----------------------
    if(diag){
      hxtheta->Fill(X,theta);
      hytheta->Fill(Y,theta);
      hxphi->Fill(X,phi);
      hyphi->Fill(Y,phi);
    }

    //mask hit--------------------------
    hmask->Fill(X,Y);   
    miss*=!(hit_mask());
    
    if(miss) {
      hhit->Fill(0);
      hmaskg->Fill(X,Y);
    }
    else {
      if(iprint)
	printf("  Hit mask!\n");
      hhit->Fill(1);
    }
  
    //calculate positions at window-----
    Z=z_window;
    trace_r(Z,theta,phi);
    trace_x(x,theta,phi);
    trace_y(y,theta,phi);
 
    if(miss)
      hwin->Fill(X,Y);
    miss*=!(hit_window());
    if(miss){//passed through
      hwing->Fill(X,Y);   
    }
    else{//hit obstruction
      if(hit_window()) {
	if(iprint)
	  printf("  Hit window!\n");
	hhit->Fill(2);
      }
    }

    //calculate positions at Y2 shield
    Z=z_A2-delta_z/2;
    trace_r(Z,theta,phi);
    trace_x(x,theta,phi);
    trace_y(y,theta,phi);
 
    miss*=!(hit_shield());
    if(miss){
      hx[3]->Fill(Y);
    }
    else{//hit obstruction
      if(hit_shield()) {
	if(iprint)
	  printf("  Hit Y2 shield!\n");
	hhit->Fill(3);
      }
    }
   
    //calculate positions at Y2 (step back)
    Z=z_Y2; 
    trace_r(Z,theta,phi);
    trace_x(x,theta,phi);
    trace_y(y,theta,phi);
    if(miss)
      hxg[3]->Fill(Y);

    //calculate positions at X2 shield
    Z=z_A2+delta_z/2;
    trace_r(Z,theta,phi);
    trace_x(x,theta,phi);
    trace_y(y,theta,phi);
 
    miss*=!(hit_shield());
    if(miss){
      hx[2]->Fill(X);
    }
    else{//hit obstruction
      if(hit_shield()) {
	if(iprint)
	  printf("  Hit X2 shield!\n");
	hhit->Fill(4);
      }
    }
   
    //calculate positions at X2
    Z=z_X2; 
    trace_r(Z,theta,phi);
    trace_x(x,theta,phi);
    trace_y(y,theta,phi);
    if(miss)
      hxg[2]->Fill(X);
    
  }//end of generator loop
}

void masks()
{
  dr("hmask");
  odr("hmaskg");
  Double_t nhit=0;
  Double_t nnohit=0;
  nhit=hmaskg->GetSum();
  nnohit=hmask->GetSum();
  printf("Percentage of tragectories passing through mask = %.2f%%\n",nhit/nnohit);
  TEllipse *ellipse = new TEllipse(0,0,99.5/2.);
  ellipse->SetLineColor(2);
  ellipse->SetLineWidth(2);
  ellipse->SetLineStyle(4);
  ellipse->SetFillStyle(0);
  ellipse->Draw();
  Float_t range=0;
  range=hmask->GetXaxis()->GetXmax();
  setvlines(-range,range);
  plotvlines(-31,0,0,0,2);
  plotvlines(-26,0,0,0,2);
  plotvlines(19,0,0,0,2);
  plotvlines(24,0,0,0,2);
  plothlines(-15,0,0,-range,range,2);
  plothlines(-10,0,0,-range,range,2);
  plothlines(10,0,0,-range,range,2);
  plothlines(15,0,0,-range,range,2);
}

void windows()
{
  dr("hwin");
  odr("hwing");
  Double_t nhit=0;
  Double_t nnohit=0;
  nhit=hwing->GetSum();
  nnohit=hwin->GetSum();
  printf("Percentage of tragectories passing through window = %.2f%%\n",nhit/nnohit);
  TEllipse *ellipse = new TEllipse(0,0,99.5/2.);
  ellipse->SetLineColor(2);
  ellipse->SetLineWidth(2);
  ellipse->SetLineStyle(4);
  ellipse->SetFillStyle(0);
  ellipse->Draw();
  Float_t range=0;
  range=hmask->GetXaxis()->GetXmax();
  setvlines(-range,range);
  plotvlines(-31,0,0,0,2);
  plotvlines(-26,0,0,0,2);
  plotvlines(19,0,0,0,2);
  plotvlines(24,0,0,0,2);
  plotvlines(0,0,34,0,2);
  plothlines(-15,0,0,-range,range,2);
  plothlines(-10,0,0,-range,range,2);
  plothlines(10,0,0,-range,range,2);
  plothlines(15,0,0,-range,range,2);
}
