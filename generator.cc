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
  TH1F *h1=new TH1F("h1","TRandom",500,0,1);
  TH1F *h2=new TH1F("h2","TRandom2",500,0,1);
  TH1F *h3=new TH1F("h3","TRandom3",500,0,1);
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

void testRandom2(Int_t nrEvents=10e+08, Float_t mean = 0, Float_t sigma = 100)
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
  
  TH1F *h1=new TH1F("h1","TRandom, Uniform",500,mean-sigma*4,mean+sigma*4);
  TH1F *h2=new TH1F("h2","TRandom, Gaus",500,mean-sigma*4,mean+sigma*4);
  TH1F *h3=new TH1F("h3","TRandom, Rndm",5000,0,1);
  
  for (Int_t i=0; i<nrEvents; i++) {
    h1->Fill(r1->Uniform(mean-sigma*4,mean+sigma*4));
    h2->Fill(r2->Gaus(mean,sigma));
    h3->Fill(r3->Rndm());
    if(i%500000==0)
      printf("%5.1f\%: %d events generated\n",(Float_t)i/nrEvents*100,i);
  }
  
  plotall("h");
  cFit->cd(2);
  gfit("h2");
}

//Beam variables
Float_t sigma_x=0.607956845;//for 90% in a 1mm radius
Float_t sigma_y=sigma_x;
Float_t offset_x=0, offset_y=0;
Bool_t is_gaussian=kTRUE;
TH2F *hbeam;

void setbeam(Float_t set_offset_x=0, Float_t set_offset_y=0, Float_t set_sigma_x=0.607956845, Float_t set_sigma_y=0.607956845, Bool_t set_is_gaussian=kTRUE)
{
  sigma_x=set_sigma_x;
  sigma_y=set_sigma_y;
  offset_x=set_offset_x;
  offset_y=set_offset_y;
  is_gaussian=set_is_gaussian;
  TString beam_distrib;
  if(is_gaussian)
    beam_distrib="Gaussian";
  else
    beam_distrib="Uniform";
  printf("Beam created at x=%5.2f, y=%5.2f\n with sigma_x=%5.2f, sigma_y=%5.2f;\n distribution is %s.\n",offset_x,offset_y,sigma_x,sigma_y,beam_distrib.Data());
  clearhists();  
}

//Detector geometry
Float_t z_mask=637.0;//position of mask
Float_t z_window=z_mask+3.25+6.35;//position of window (back edge)
Float_t z_A2=z_window+12.7+23.6;//position of A2
Float_t z_A1=z_A2+36.8;//anode separation
Float_t delta_z=0.125*25.4;//plane separation

Float_t z_Y2=z_A2-delta_z;
Float_t z_X2=z_A2+delta_z;
Float_t z_Y1=z_A1-delta_z;
Float_t z_X1=z_A1+delta_z;
Float_t z_Y2s=z_A2-delta_z/2;
Float_t z_X2s=z_A2+delta_z/2;
Float_t z_Y1s=z_A1-delta_z/2;
Float_t z_X1s=z_A1+delta_z/2;

Float_t x=0,y=0;//point on target
Float_t theta=0;//scattering angle
Float_t phi=0;//azimuthal angle
Float_t X=0, Y=0, Z=0, Rho=0;//position
Float_t Xm=0, Ym=0, Xr=0, Yr=0;
Float_t r=0,rho=0;//radii
Bool_t doprint=1;//kFALSE;
Bool_t diag=kFALSE;

void rotate_axis(Float_t z_start,Float_t theta_start, Float_t phi_start)
{//rotate system about the y-axis by theta_center
  //first, calculate positions
  trace_r(z_start,theta_start,phi_start);
  trace_x(x,theta_start,phi_start);
  trace_y(y,theta_start,phi_start);
  //next, calculate the transformed coordinates
  if(iprint) 
    printf("  Rotating by %7.2f...\n",theta_center);
  Float_t XX=0,YY=0,ZZ=0;
  XX=X*TMath::Cos(theta_center*TMath::DegToRad())-z_start*TMath::Sin(theta_center*TMath::DegToRad());
  YY=Y;
  ZZ=X*TMath::Sin(theta_center*TMath::DegToRad())+z_start*TMath::Cos(theta_center*TMath::DegToRad());
  r=TMath::Sqrt(XX*XX+YY*YY+ZZ*ZZ);
  if(iprint) 
    printf("   New positions are X=%7.2f Y=%7.2f Z=%7.2f r=%7.2f\n",XX,YY,ZZ,r);  
  //export position to global variables
  X=XX;
  Y=YY;
  Z=ZZ; 
  //then, calculate the new trajectories
  Float_t old_theta=theta_start;
  Float_t old_phi=phi_start;
  theta=TMath::RadToDeg()*TMath::ACos(ZZ/r);
  phi=TMath::RadToDeg()*TMath::ATan(YY/XX);
  if(XX<0) phi+=180;
  if(YY<0&&XX>0) phi+=360;
  if(iprint) {
    printf("   New angles are (theta,phi)=(%7.2f,%7.2f)\n",theta,phi);  
    if(fabs(old_theta-theta)>1e-2)
      printf("---Theta mis-match! old_theta=%7.2f theta=%7.2f\n",old_theta,theta);
    if(fabs(old_phi-phi)>1e-2)
      printf("---Phi mis-match! old_phi=%7.2f phi=%7.2f\n",old_phi,phi);
  }
}

void trace_r(Float_t z_start,Float_t theta_start, Float_t phi_start)
{//given z, theta, phi; calculate r
  r=z_start/(TMath::Cos(theta_start*TMath::DegToRad()));
  rho=r*(TMath::Sin(theta_start*TMath::DegToRad()));
  if(iprint) {
    printf("  Emmission angle (theta,phi)=(%7.2f,%7.2f)\n",theta_start,phi_start);
    printf("   Current position is Z=%7.2f\n",z_start);
    printf("   Spherical radius is r=%7.2f\n",r);
    printf("   Cylindrical rad. is rho=%7.2f\n",rho);
  }
}

void trace_x(Float_t x_start, Float_t theta_start, Float_t phi_start)
{
  X=rho*(TMath::Cos(phi_start*TMath::DegToRad()));
  if(iprint) {
    printf("   X-position is %7.2f (relative), with offset %7.2f\n",X,x_start);
  }
  X+=x_start*TMath::Cos(theta_center*TMath::DegToRad());  
  if(iprint) {
    printf("   X-position is %7.2f (absolute)\n",X);
  }
}

void trace_y(Float_t y_start, Float_t theta_start, Float_t phi_start)
{
  Y=rho*(TMath::Sin(phi_start*TMath::DegToRad()));
  if(iprint) {
    printf("   Y-position is %7.2f (relative), with offset %7.2f\n",Y,y_start);
  }
  Y+=y_start;  
  if(iprint) {
    printf("   Y-position is %7.2f (absolute)\n",Y);
  }
}

Float_t slit_width=0.5;
Float_t slits[5]={-20,-10,0,5,10};
//the position of the slits are not needed outside fo the function, but delcaring the variables each run is very slow
Bool_t hit_slits(void)
{
  for(int i=0; i<5; i++) {
    if((X>(slits[i]-slit_width/2))&&(X<(slits[i]+slit_width/2)))
      return kFALSE;
  }
  return kTRUE;
}

Bool_t hit_mask(void)
{
  Rho=TMath::Sqrt((TMath::Power(X,2))+(TMath::Power(Y,2)));
  if(Rho>(99.5/2)) 
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

Float_t y_shield=28;
Bool_t hit_shield(void)
{
  if((X<-118)||(X>36))
    return kTRUE;//hits shield (x)
  if((Y<-y_shield)||(Y>y_shield))
    return kTRUE;//hist shield (y)
  return kFALSE;
}

Float_t x_win=32;//34.0073412589711;//nominally 34.01
Float_t y_win=25;
void setwin(Float_t set_x_win=34.0073412589711,Float_t set_y_win=100)
{
  x_win=set_x_win;
  y_win=set_y_win;
  x_feature[6]=x_win;
  clearhists();  
}

Bool_t hit_window(void)
{
  if(X>x_win||Y>y_win||Y<-y_win)
    return kTRUE;
  return kFALSE;
}

Float_t theta_min=0;
Float_t theta_max=180;
Float_t theta_center=30;
Float_t phi_min=0;
Float_t phi_max=360;
Bool_t is_rutherford=kFALSE;

void setangles(Float_t set_theta_min=0, Float_t set_theta_max=180, Float_t set_phi_min=0, Float_t set_phi_max=360, Float_t set_theta_center =0, Bool_t set_is_rutherford=kFALSE)
{//define angles and build histograms
  theta_min=set_theta_min;
  theta_max=set_theta_max;
  phi_min=set_phi_min;
  phi_max=set_phi_max;
  theta_center=set_theta_center;
  printf("Emmission angles defined over:\n Theta %f to %f\n Cos(theta) %f to %f\n Phi %f to %f\n",theta_min,theta_max,TMath::Cos(theta_min*TMath::DegToRad()),TMath::Cos(theta_max*TMath::DegToRad()),phi_min,phi_max);
  is_rutherford=set_is_rutherford;
  TString ang_distrib;
  if(is_rutherford) {
    ang_distrib="Rutherford scattering cross section";
    rutherdef(theta_min,theta_max);
  }
  else
    ang_distrib="uniform (over unit sphere)";
  printf(" Angular distribution is %s.\n",ang_distrib.Data());
  printf(" Central angle is %f.\n",theta_center);
  clearhists();
}

TH1F *htheta;
TH1F *hphi;
TH2F *hphitheta;
TH1F *hcostheta;
TH2F *hphicostheta;

TH2F *hxtheta;
TH2F *hytheta;
TH2F *hxphi;
TH2F *hyphi;

TH2F *hmask;
TH2F *hmaskg;
TH1F *hhit;
TH1F *hmiss;
TH1F *hnewhit;

TH2F *hwin;
TH2F *hwing;

TH1F *hx[4];
TH1F *hxg[4];
TH1F *hxgm[4];
TH1F *hxgmr[4];
TH2F *hyx[4];
TH2F *hyxg[4];
TH2F *hyxgm[2];
TH2F *hyxgmr[2];

Float_t x_min=0;
Float_t x_max=100;
Float_t x_cal=118;//location of central beam axis, relative to edge of shield; not to be confused with offset_x
Float_t y_cal=54/2;// not to be confused with offset_y

void definehists()
{
  printf("Defining histograms...\n");
  Float_t b_max=0;
  Float_t b_width=7;
  Float_t b_bin=500;//beam spot bin
  Float_t x_bin=b_bin;//position bin
  Float_t a_bin=100;//angle bin

  if((sigma_x==0)&&(sigma_y==0)){
    b_max=0;
    b_bin=3;
  }
  else
    if(sigma_x>sigma_y)
      b_max=sigma_x;
    else
      b_max=sigma_y;
  b_max*=b_width;
  if((fabs(offset_x))>(fabs(offset_y)))
    b_max+=fabs(offset_x);
  else
    b_max+=fabs(offset_y);
  hbeam=new TH2F("hbeam","Beam Spot",b_bin,-b_max,b_max,b_bin,-b_max,b_max);
  hbeam->SetXTitle("x-position (mm)");
  hbeam->SetYTitle("y-position (mm)"); 
  printf(" Beam histogram limits are bmin=%7.2f, bmax=%7.2f\n",-b_max,b_max);
  
  htheta=new TH1F("htheta","#theta (Polar Angle)",a_bin,0,180);
  hphi=new TH1F("hphi","#phi (Azimuthal Angle)",a_bin,0,360);
  hphitheta=new TH2F("hphitheta","#phi vs. #theta",a_bin,0,180,a_bin,0,360);
  hphitheta->SetXTitle("#theta - Polar Angle (deg)");
  hphitheta->SetYTitle("#phi - Azimuthal Angle (deg)");
  hcostheta=new TH1F("hcostheta","Cos(#theta) - Cos Polar Angle",a_bin,-1,1);
  hphicostheta=new TH2F("hphicostheta","#phi vs. cos(#theta)",a_bin,-1,1,a_bin,0,360);

  x_max=z_X1*(TMath::Tan(theta_max*TMath::DegToRad()));
  Float_t offset_max=0;
  if((fabs(offset_x))>(fabs(offset_y)))
    offset_max+=fabs(offset_x);
  else
    offset_max+=fabs(offset_y);
  x_max+=offset_max;
  x_max*=1.05;
  x_max=ceil(x_max);
  Float_t x_Max=124;
  if((x_max>x_Max)||(x_max<0))
    x_max=x_Max;
  x_min=-x_max;
  
  printf(" Position histogram limits are xmin=%7.2f, xmax=%7.2f\n",x_min,x_max);
  printf(" to cover +/- %5.1f deg + %f mm (x1.05)\n",theta_max,offset_max);

  hmask=new TH2F("hmask" ,"Mask Plane",x_bin,x_min,x_max,x_bin,x_min,x_max);
  hmaskg=new TH2F("hmaskg","Mask Plane (gated)",x_bin,x_min,x_max,x_bin,x_min,x_max);
  hhit=new TH1F("hhit","hhit",7,-0.5,6.5);
  hmiss=new TH1F("hmiss","hmiss",7,-0.5,6.5);
  hnewhit=new TH1F("hnewhit","hnewhit",7,-0.5,6.5);
  hwin=new TH2F("hwin" ,"Window Plane",x_bin,x_min,x_max,x_bin,x_min,x_max);
  hwing=new TH2F("hwing","Window Plane (gated)",x_bin,x_min,x_max,x_bin,x_min,x_max);

  if(diag) {
    //Other correlation plots----------------------
    hxtheta=new TH2F("hxtheta","#theta (polar) vs. X",x_bin,x_min,x_max,a_bin,0,180);
    hxtheta->SetYTitle("#theta - polar angle (deg)");
    hytheta=new TH2F("hytheta","#theta (polar) vs. Y",x_bin,x_min,x_max,a_bin,0,180);
    hytheta->SetYTitle("#theta - polar angle (deg)");
    hxphi=new TH2F("hxphi","#phi (azimuthal) vs. X",x_bin,x_min,x_max,a_bin,0,360);
    hxphi->SetYTitle("#phi - azimuth angle (deg)");
    hyphi=new TH2F("hyphi","#phi (azimuthal) vs. Y",x_bin,x_min,x_max,a_bin,0,360);
    hyphi->SetYTitle("#phi - azimuth angle (deg)");
  }

  hx[0] = new TH1F("hx0","X1 Position",x_bin,x_min,x_max);
  hx[1] = new TH1F("hx1","Y1 Position",x_bin,x_min,x_max);
  hx[2] = new TH1F("hx2","X2 Position",x_bin,x_min,x_max);
  hx[3] = new TH1F("hx3","Y2 Position",x_bin,x_min,x_max);
  hxg[0] = new TH1F("hxg0","X1 Position (gated)",x_bin,x_min,x_max);
  hxg[1] = new TH1F("hxg1","Y1 Position (gated)",x_bin,x_min,x_max);
  hxg[2] = new TH1F("hxg2","X2 Position (gated)",x_bin,x_min,x_max);
  hxg[3] = new TH1F("hxg3","Y2 Position (gated)",x_bin,x_min,x_max);
  hxgm[0] = new TH1F("hxgm0","X1 Position (gated)",x_bin,x_min+x_cal,x_max+x_cal);
  hxgm[1] = new TH1F("hxgm1","Y1 Position (gated)",x_bin,x_min+y_cal,x_max+y_cal);
  hxgm[2] = new TH1F("hxgm2","X2 Position (gated)",x_bin,x_min+x_cal,x_max+x_cal);
  hxgm[3] = new TH1F("hxgm3","Y2 Position (gated)",x_bin,x_min+y_cal,x_max+y_cal);
  hxgmr[0] = new TH1F("hxgmr0","X1 Position (gated)",x_bin,x_min+x_cal,x_max+x_cal);
  hxgmr[1] = new TH1F("hxgmr1","Y1 Position (gated)",x_bin,x_min+y_cal,x_max+y_cal);
  hxgmr[2] = new TH1F("hxgmr2","X2 Position (gated)",x_bin,x_min+x_cal,x_max+x_cal);
  hxgmr[3] = new TH1F("hxgmr3","Y2 Position (gated)",x_bin,x_min+y_cal,x_max+y_cal);
  /* for(int i = 0; i < 4; i++) {//1D position plots
     hx[i]->SetXTitle("Relative Position");
     hx[i]->SetYTitle("Number of Entries");
     }*/
  hyx[0] = new TH2F("hyx0","Y vs X at X1",x_bin,x_min,x_max,x_bin,x_min,x_max);
  hyx[1] = new TH2F("hyx1","Y vs X at Y1",x_bin,x_min,x_max,x_bin,x_min,x_max);
  hyx[2] = new TH2F("hyx2","Y vs X at X2",x_bin,x_min,x_max,x_bin,x_min,x_max);
  hyx[3] = new TH2F("hyx3","Y vs X at Y2",x_bin,x_min,x_max,x_bin,x_min,x_max);
  hyxg[0] = new TH2F("hyxg0","Y vs X at X1 (gated on X1 Shield)",x_bin,x_min,x_max,x_bin,x_min,x_max);
  hyxg[1] = new TH2F("hyxg1","Y vs X at Y1 (gated on Y1 Shield)",x_bin,x_min,x_max,x_bin,x_min,x_max);
  hyxg[2] = new TH2F("hyxg2","Y vs X at X2 (gated on X2 Shield)",x_bin,x_min,x_max,x_bin,x_min,x_max);
  hyxg[3] = new TH2F("hyxg3","Y vs X at Y2 (gated on Y2 Shield)",x_bin,x_min,x_max,x_bin,x_min,x_max);
  hyxgm[0] = new TH2F("hyxgm0","Y1 vs X1 Positions (gated), measured",x_bin,x_min+x_cal,x_max+x_cal,x_bin,x_min+y_cal,x_max+y_cal);
  hyxgm[1] = new TH2F("hyxgm1","Y2 vs X2 Positions (gated), measured",x_bin,x_min+x_cal,x_max+x_cal,x_bin,x_min+y_cal,x_max+y_cal);
  hyxgmr[0] = new TH2F("hyxgmr0","Y1 vs X1 Positions (gated), measured, blurred",x_bin,x_min+x_cal,x_max+x_cal,x_bin,x_min+y_cal,x_max+y_cal);
  hyxgmr[1] = new TH2F("hyxgmr1","Y2 vs X2 Positions (gated), measured, blurred",x_bin,x_min+x_cal,x_max+x_cal,x_bin,x_min+y_cal,x_max+y_cal);
}

void clearhists()
{
  printf(" Clearing histograms...\n");
  const int Nhists = 45;
  TString histnames[Nhists]={"hbeam","htheta","hphi","hphitheta","hmask","hmaskg","hxtheta","hytheta","hxphi","hyphi","hhit","hx0","hx1","hx2","hx3","hwin","hwing","hxg0","hxg1","hxg2","hxg3","hyx0","hyx1","hyxg0","hyxg1","hyxgm0","hyxgm1","hmiss","hnewhit","hyx2","hyx3","hyxg2","hyxg3","hcostheta","hphicostheta","hyxgmr0","hyxgmr1","hxgm0","hxgm1","hxgm2","hxgm3","hxgmr0","hxgmr1","hxgmr2","hxgmr3"};

  for(int i=0; i < Nhists; i++) {
    if (gROOT->FindObject(histnames[i])) {
      if(doprint)
	printf(" Histogram %10s already exits.  Clearing...\n",histnames[i].Data());
      gROOT->FindObject(histnames[i])->Clear();    
    }
  }
}

Float_t xres=0, yres=0;
void setres(Float_t set_xres=0, Float_t set_yres=0)
{
  xres=set_xres;
  yres=set_yres;
  printf("Dectector resolution set to %f (x) and %f (y)\n",xres,yres);
  printf("                       %f FWHM (x) and %f FWHM (y)\n",xres*2.35482,yres*2.35482);
  clearhists();
}

Bool_t iprint=kFALSE; //doprint;  
Bool_t donewhit=kFALSE;
Bool_t doslits=kFALSE;
void source(Int_t nevents=5e4, Float_t weight=1, Bool_t set_doprint=kFALSE)
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
  Float_t cos_theta_min=(TMath::Cos(theta_min*TMath::DegToRad()));
  Float_t cos_theta_max=(TMath::Cos(theta_max*TMath::DegToRad()));
   
  //azimuthal angle-------------------------------
  TRandom3 *rphi=new TRandom3();
  rphi->SetSeed(0);
 
  //measurement-----------------------------------
  TRandom3 *rxres=new TRandom3();//x resolution
  TRandom3 *ryres=new TRandom3();//y resoultion
  rxres->SetSeed(0);
  ryres->SetSeed(0);
  //jucntion peaks------------
  Float_t yjunk1l=22-6-y_cal-.5;
  Float_t yjunk1r=44-6-y_cal-.3;
  Float_t yjunk2l=22-6-y_cal;
  Float_t yjunk2r=44-6-y_cal-.55;
  Float_t junk_width1=1.7;
  Float_t junk_width2=0.41;
  Float_t junk_res=.4;

  //X, Y positions (ray-tracing)------------------
  Bool_t hit=kFALSE;
  Bool_t miss=kTRUE;
   
  //step size to print updates 
  int step=(int) (nevents/10);
  int step_max=5e3;
  if(step>step_max)
    step=step_max;  
    
  //create events!
  for (Int_t i=0; i<nevents; i++) {
    if(i%step==0) {
      printf("%5.1f%%: %d events generated\n",(Float_t)i/nevents*100,i);
      iprint=doprint*kTRUE;
    }
    else
      iprint=kFALSE;  
    hit=kFALSE;//hit means an obstruction has been hit
    miss=kTRUE;//miss means the trajectory is unobstructed
    //Position of origin (beam spot)-----------------
    if(is_gaussian) {
      x=rx->Gaus(offset_x,sigma_x);
      y=ry->Gaus(offset_y,sigma_y);
    }
    else {
      x=rx->Uniform(-sigma_x,sigma_x);
      y=ry->Uniform(-sigma_y,sigma_y);
    }
    hbeam->Fill(x,y);
    
    //Emmission angle-----------------------------
    // Polar angle----------------------
    if(!is_rutherford)
      theta=(TMath::ACos(rtheta->Uniform(cos_theta_min,cos_theta_max))*(TMath::RadToDeg()));
    else {
         theta=ruther->GetRandom();
       }
    htheta->Fill(theta);
    hcostheta->Fill(TMath::Cos(theta*TMath::DegToRad()));
    // Azimuthal angle------------------
    phi=rphi->Uniform(phi_min,phi_max);
    hphi->Fill(phi);
    hphitheta->Fill(theta,phi);
    hphicostheta->Fill(TMath::Cos(theta*TMath::DegToRad()),phi);
    // Diagnostics----------------------
    if(diag){
      hxtheta->Fill(X,theta);
      hytheta->Fill(Y,theta);
      hxphi->Fill(X,phi);
      hyphi->Fill(Y,phi);
    } 
    hmiss->Fill(0);//initial number of particles
       
    //calculate positions at mask-------
    Z=z_mask;
    rotate_axis(Z,theta,phi);
    trace_r(Z,theta,phi);
    trace_x(x,theta,phi);
    trace_y(y,theta,phi);
   
    //mask hit--------------------------
    hmask->Fill(X,Y);//position at mask before hit
    miss*=!(hit_mask());
    
    if(miss) {//passes through mask
      hmiss->Fill(1);
      hmaskg->Fill(X,Y);
     
      //calculate positions at window-----
      Z=z_window;
      trace_r(Z,theta,phi);
      trace_x(x,theta,phi);
      trace_y(y,theta,phi);
      hwin->Fill(X,Y);
    }
    miss*=!(hit_window());
    if(miss) {//passes through mask, window
      hmiss->Fill(2);
      hwing->Fill(X,Y);   
      //------------------------------------------------------
      //Detector 2----------------------------------
      //calculate positions at Y2 shield before calculating hits
      Z=z_Y2s;
      trace_r(Z,theta,phi);
      trace_x(x,theta,phi);
      trace_y(y,theta,phi);
      
      hx[3]->Fill(Y);
      hyx[3]->Fill(X,Y);
    }
    miss*=!(hit_shield());
    if(miss){//passes through mask, window, Y2 shield
      hmiss->Fill(3); 
    
      //calculate positions at Y2 (step back)
      Z=z_Y2; 
      trace_r(Z,theta,phi);
      trace_x(x,theta,phi);
      trace_y(y,theta,phi);    
      hxg[3]->Fill(Y);
      hyxg[3]->Fill(X,Y);
    
      //calculate positions at X2 shield before calculating hits
      Z=z_X2s;
      trace_r(Z,theta,phi);
      trace_x(x,theta,phi);
      trace_y(y,theta,phi);   
      
      hx[2]->Fill(X);
      hyx[2]->Fill(X,Y);
    }
    miss*=!(hit_shield());
    if(miss){//passes through mask, window, Y2 shield, X2 shield (detected at A2)
      hmiss->Fill(4);
   
      //calculate positions at X2 (step forward)
      Z=z_X2; 
      trace_r(Z,theta,phi);
      trace_x(x,theta,phi);
      trace_y(y,theta,phi);     
      hxg[2]->Fill(X);
      hyxg[2]->Fill(X,Y);
   
      //calculate positions *measured* at anode (step back)
      Z=z_Y2s;
      trace_r(Z,theta,phi);
      trace_y(y,theta,phi);
      Z=z_X2s;
      trace_r(Z,theta,phi);
      trace_x(x,theta,phi);
      
    }
    if(doslits)
      miss*=!(hit_slits());
    if(miss) {//passes through slits (if pressent)
      hxgm[3]->Fill(Y+y_cal);
      hxgm[2]->Fill(X+x_cal);
      hyxgm[1]->Fill(X+x_cal,Y+y_cal);
      //calculate measured positions with given detector resolution
      Xr=rxres->Gaus(X,xres);
      hxgmr[2]->Fill(Xr+x_cal,weight);      
      if(((Y<yjunk2l+junk_width2/2)&&(Y>yjunk2l-junk_width2/2))||((Y<yjunk2r+junk_width2/2)&&(Y>yjunk2r-junk_width2/2))) {
       	for (Int_t j=0; j<2; j++) {
	  Yr=ryres->Gaus(Y,junk_res);
	  hxgmr[3]->Fill(Yr+y_cal,weight);
	  hyxgmr[1]->Fill(Xr+x_cal,Yr+y_cal);
	}
      }
      // else{
      Yr=ryres->Gaus(Y,yres);
      hxgmr[3]->Fill(Yr+y_cal,weight);
      hyxgmr[1]->Fill(Xr+x_cal,Yr+y_cal);
      //  }
      
      //------------------------------------------------------
      //Detector 1----------------------------------
      //calculate positions at Y1 shield before calculating hits
      Z=z_Y1s;
      trace_r(Z,theta,phi);
      trace_x(x,theta,phi);
      trace_y(y,theta,phi);
    
      hx[1]->Fill(Y);
      hyx[1]->Fill(X,Y);
    }
    miss*=!(hit_shield());
    if(miss){//passes through Y1 shield
      hmiss->Fill(5); 
   
      //calculate positions at Y1 (step back)
      Z=z_Y1; 
      trace_r(Z,theta,phi);
      trace_x(x,theta,phi);
      trace_y(y,theta,phi);   
      hxg[1]->Fill(Y);
      hyxg[1]->Fill(X,Y);
  
      //calculate positions at X1 shield before calculating hits
      Z=z_X1s;
      trace_r(Z,theta,phi);
      trace_x(x,theta,phi);
      trace_y(y,theta,phi);    
      hx[0]->Fill(X);
      hyx[0]->Fill(X,Y);
    }
    miss*=!(hit_shield());
    if(miss) {//passes through Y1 shield, X1 shield (detected at A1)
      hmiss->Fill(6);
    
      //calculate positions at X1 (step forward)
      Z=z_X1; 
      trace_r(Z,theta,phi);
      trace_x(x,theta,phi);
      trace_y(y,theta,phi);    
      hxg[0]->Fill(X);
      hyxg[0]->Fill(X,Y);

      //calculate positions *measured* at anode (step back)
      Z=z_Y1s;
      trace_r(Z,theta,phi);
      trace_y(y,theta,phi);
      Z=z_X1s;
      trace_r(Z,theta,phi);
      trace_x(x,theta,phi);

      hxgm[1]->Fill(Y+y_cal);
      hxgm[0]->Fill(X+x_cal);
      hyxgm[0]->Fill(X+x_cal,Y+y_cal);
      //calculate measured positions with given detector resolution
      Xr=rxres->Gaus(X,xres);
      hxgmr[0]->Fill(Xr+x_cal,weight);     
      if(((Y<yjunk1l+junk_width1/2)&&(Y>yjunk1l-junk_width1/2))||((Y<yjunk1r+junk_width1/2)&&(Y>yjunk1r-junk_width1/2))) {
       	for (Int_t j=0; j<1; j++) {
	  Yr=ryres->Gaus(Y,junk_res);
	  hxgmr[1]->Fill(Yr+y_cal,weight);
	  hyxgmr[0]->Fill(Xr+x_cal,Yr+y_cal);
	}
      }
      //   else{
      Yr=ryres->Gaus(Y,yres);
      hxgmr[1]->Fill(Yr+y_cal,weight);
      hyxgmr[0]->Fill(Xr+x_cal,Yr+y_cal);
      //   }
    }
  }//end of generator loop
}
Float_t range=0;
void mask(Char_t *histin1, Char_t *histin2)
{
  gate(histin1,histin2)
    TEllipse *ellipse = new TEllipse(0,0,99.5/2.);
  ellipse->SetLineColor(2);
  ellipse->SetLineWidth(2);
  ellipse->SetLineStyle(4);
  ellipse->SetFillStyle(0);
  ellipse->Draw();
  plotvlines(-31,0,0,0,2);
  plotvlines(-26,0,0,0,2);
  plotvlines(19,0,0,0,2);
  plotvlines(24,0,0,0,2);
  plothlines(-15,0,0,-range,range,2);
  plothlines(-10,0,0,-range,range,2);
  plothlines(10,0,0,-range,range,2);
  plothlines(15,0,0,-range,range,2);
}

void masks()
{
  gate("hmask","hmaskg");
  maskz(z_mask);
}

void windows()
{
  gate("hwin","hwing");
  maskz(z_window);
  windowz(z_window);
}

void battleship()
{
  mkCanvas2("cHits","cHits");
  cHits->cd();
  cHits->Divide(1,2);
  cHits->cd(1);
  hmiss->Draw();
  odr("hhit");
  cHits->cd(2);
  hnewhit->Draw();
}

void y2s(Float_t z_plane=z_Y2s)
{
  gate("hyx3","hyxg3");  
  shadowz(z_plane);
}

void x2s(Float_t z_plane=z_X2s)
{
  gate("hyx2","hyxg2");
  shadowz(z_plane);
}

void y1s(Float_t z_plane=z_Y1s)
{
  gate("hyx1","hyxg1"); 
  shadowz(z_plane);
}

void x1s(Float_t z_plane=z_X1s)
{
  gate("hyx0","hyxg0");
  shadowz(z_plane);
}

void shadowz(Float_t z_plane=0, Bool_t docal=kFALSE)
{
  printf("Calculating positions at z=%7.2f\n with calibration offsets of x=%f, y=%f\n",z_plane,x_cal,y_cal);
  maskz(z_plane,docal);
  //if(!docal)
    windowz(z_plane,docal);
  shieldz(z_plane,docal);
  printf("X-gaps are centered at %6.3f (%6.3f wide) and %6.3f (%6.3f wide)\n",(x_shadow[1]+x_shadow[2])/2,(x_shadow[2]-x_shadow[1]),(x_shadow[3]+x_shadow[4])/2,(x_shadow[4]-x_shadow[3]));
  printf("X-spans are %6.3f, %6.3f, and ",(x_shadow[1]-x_shadow[0]),(x_shadow[3]-x_shadow[2]));
  Float_t right_wide=0;
  right_wide=(x_shadow[5]-x_shadow[4]);
  if(right_wide>(x_shadow[6]-x_shadow[4]))
    right_wide=(x_shadow[6]-x_shadow[4]);
  printf("%6.3f wide (%6.3f total)\n",right_wide,(x_shadow[5]-x_shadow[0]));

  printf("Y-gaps are centered at %6.3f (%6.3f wide) and %6.3f (%6.3f wide)\n",(y_shadow[1]+y_shadow[2])/2,(y_shadow[2]-y_shadow[1]),(y_shadow[3]+y_shadow[4])/2,(y_shadow[4]-y_shadow[3]));
  printf("Y-spans are %6.3f, %6.3f, and %6.3f wide (%6.3f total)\n",(y_shadow[1]-y_shadow[0]),(y_shadow[3]-y_shadow[2]),(y_shadow[5]-y_shadow[4]),(y_shadow[5]-y_shadow[0]));
}

void shield(Char_t *histin1, Char_t *histin2)
{
  gate(histin1,histin2);
  plotvlines(-118,0,0,0,2);
  plotvlines(36,0,0,0,2);
  plothlines(-27,0,0,-range,range,2);
  plothlines(27,0,0,-range,range,2);
}

Float_t x_feature[7]={-49.75, -31, -26, 19, 24, 36, x_win};
Float_t y_feature[6]={-27, -15, -10, 10, 15, 27};
Float_t theta_proj[7]=0;
Float_t x_shadow[7]=0;
Float_t y_shadow[6]=0;

void maskz(Float_t z_plane, Bool_t docal=kFALSE)
{//calculate position of mask features at position z
  Float_t z_feature=z_mask;
  FILE * outfile_x, * outfile_y;
  outfile_y=fopen("temp_Y.lst","w");
  outfile_x=fopen("temp_X.lst","w");
 
  printf("Positions in y-direction of mask\n");
  for(int i=1; i<5; i++){
    printf(" %d z_feature = %5.1f ",i+1, z_feature);
    theta_proj[i]=TMath::ATan((y_feature[i]-offset_y)/z_feature);
    printf("theta_proj = %6.3f ",(TMath::RadToDeg()*theta_proj[i]));
    y_shadow[i]=z_plane*(TMath::Tan(theta_proj[i]))+offset_y;
    printf("y_shadow = %7.2f (%7.2f)\n", y_shadow[i],y_shadow[i]+y_cal);
    fprintf(outfile_y,"%f\n",y_shadow[i]+y_cal);
    if(docal)
      plothlines(y_shadow[i]+y_cal,0,0,range-span,range,2);   
    else
      plothlines(y_shadow[i],0,0,-range,range,2);   
  }
  fclose(outfile_y);
  printf("Positions in x-direction of mask\n");
  for(int i=0; i<5; i++) {
    printf(" %d z_feature = %5.1f ",i+1, z_feature);
    theta_proj[i]=TMath::ATan((x_feature[i]-offset_x)/z_feature);
    printf("theta_proj = %6.3f ",(TMath::RadToDeg()*theta_proj[i]));
    x_shadow[i]=z_plane*(TMath::Tan(theta_proj[i]))+offset_x;
    printf("x_shadow = %7.2f (%7.2f)\n", x_shadow[i],x_shadow[i]+x_cal);
    fprintf(outfile_x,"%f\n",x_shadow[i]+x_cal);    
    if(docal){
      if(i>0)
	plotvlines(x_shadow[i]+x_cal,0,0,0,2);
    }
    else
      plotvlines(x_shadow[i],0,0,0,2);
  }
  fclose(outfile_x);

  if(docal)
    TEllipse *emask = new TEllipse(0+x_cal,0+y_cal,-x_shadow[0]);
  else {
    TEllipse *emask = new TEllipse(0,0,-x_shadow[0]);
    //beam envelope
    Float_t temp_x=0;
    rotate_axis(z_mask,0,0);
    temp_x=X;
    rotate_axis(z_mask,theta_max,0);
    TEllipse *econe = new TEllipse(temp_x+offset_x,Y+offset_y,z_plane*TMath::Tan(theta*TMath::DegToRad()));
    econe->SetLineColor(3);
    econe->SetLineWidth(2);
    econe->SetLineStyle(4);
    econe->SetFillStyle(0);
    econe->Draw();
    rotate_axis(z_mask,theta_min,0);
    TEllipse *econe2 = new TEllipse(temp_x+offset_x,Y+offset_y,z_plane*TMath::Tan(theta*TMath::DegToRad()));
    econe2->SetLineColor(4);
    econe2->SetLineWidth(2);
    econe2->SetLineStyle(4);
    econe2->SetFillStyle(0);
    econe2->Draw();
  }
  emask->SetLineColor(2);
  emask->SetLineWidth(2);
  emask->SetLineStyle(4);
  emask->SetFillStyle(0);
  emask->Draw();
}

void windowz(Float_t z_plane, Bool_t docal=kFALSE)
{//plot the position of the window edge for a given z
  Float_t z_feature=z_window;
  Int_t i=6;
  printf(" %d z_feature = %5.1f ",i+1, z_feature);
  theta_proj[i]=TMath::ATan((x_feature[i]-offset_x)/z_feature);
  printf("theta_proj = %6.3f ",(TMath::RadToDeg()*theta_proj[i]));
  x_shadow[i]=z_plane*(TMath::Tan(theta_proj[i]))+offset_x;
  printf("x_shadow = %7.2f (%7.2f)\n", x_shadow[i],x_shadow[i]+x_cal);
  if(docal)
    plotvlines(0,x_shadow[i]+x_cal,0,0,4);
  else
    plotvlines(0,x_shadow[i],0,0,4);
}

void shieldz(Float_t z_plane, Bool_t docal=kFALSE, Int_t setlinecolor=3)
{
  Float_t z_feature=0;
  if(z_plane>z_A1)
    z_feature=z_A1+delta_z/2;
  else
    if(z_plane>z_Y1)
      z_feature=z_A1-delta_z/2;
    else
      if(z_plane>z_A2)
	z_feature=z_A2+delta_z/2;
      else
	z_feature=z_A2-delta_z/2;
  printf("Positions in y-direction of shields\n");  
  for(int i=0; i<6; i+=5) {
    printf(" %d z_feature = %5.1f ",i+1, z_feature);
    theta_proj[i]=TMath::ATan((y_feature[i]-offset_y)/z_feature);
    printf("theta_proj = %6.3f ",(TMath::RadToDeg()*theta_proj[i]));
    y_shadow[i]=z_plane*(TMath::Tan(theta_proj[i]))+offset_y;
    printf("y_shadow = %7.2f (%7.2f)\n", y_shadow[i],y_shadow[i]+y_cal);
    if(docal)   
      plothlines(0,0,y_shadow[i]+y_cal,range-span,range,setlinecolor);   
    else
      plothlines(y_shadow[i],0,0,-range,range,6);   
  }
  printf("Positions in x-direction of shields\n");  
  for(int i=5; i<6; i++) {
    printf(" %d z_feature = %5.1f ",i+1, z_feature);
    theta_proj[i]=TMath::ATan((x_feature[i]-offset_x)/z_feature);
    printf("theta_proj = %6.3f ",(TMath::RadToDeg()*theta_proj[i]));
    x_shadow[i]=z_plane*(TMath::Tan(theta_proj[i]))+offset_x;
    printf("x_shadow = %7.2f (%7.2f)\n", x_shadow[i], x_shadow[i]+x_cal);
    if(docal)
      plotvlines(0,0,x_shadow[i]+x_cal,0,setlinecolor);
    else
      plotvlines(x_shadow[i],0,0,0,6);
  }
}

void gate(Char_t *histin1, Char_t *histin2)
{
  TH2F *hist1=(TH2F *) gROOT->FindObject(histin1);
  TH2F *hist2=(TH2F *) gROOT->FindObject(histin2);
  dr(histin1);
  odr(histin2);
  Float_t nstart=0;
  Float_t npass=0;
  nstart=hist1->Integral();
  npass=hist2->Integral();
  printf("%9.2f entries in ungated spectrum\n",nstart);
  printf("%9.2f entries in gated spectrum\n",npass);
  printf("Percentage of tragectories passing through window = %.2f%%\n",100*npass/nstart);
  range=hist1->GetXaxis()->GetXmax();
  setvlines(-range,range);
}

Float_t span=0;
void shadowzc(Char_t *histin1,Float_t z_plane=0)
{
  TH2F *hist1=(TH2F *) gROOT->FindObject(histin1);
  dr(histin1);
  range=hist1->GetXaxis()->GetXmax();
  span=range-(hist1->GetXaxis()->GetXmin());
  setvlines(hist1->GetYaxis()->GetXmin(),hist1->GetYaxis()->GetXmax());
  shadowz(z_plane, kTRUE);
}

void setsim(Float_t set_res=0, Float_t set_beam=0.607956845, Float_t set_events=1e5, Float_t set_weight=1)
{
  setres(set_res,set_res);
  setbeam(0,0,set_beam,set_beam);
  source(set_events,set_weight);
}

void compY2(Float_t set_res=0, Float_t set_beam=0.607956845, Float_t set_events=1e5)
{
  if((xres==set_res)&&(yres==set_res))
    {
      printf("Resolutions match\n");
      if((sigma_x==set_beam)&&(sigma_y==set_beam))
	printf("...and beam spot sizes match.\n");
      else
	setsim(set_res,set_beam,set_events);  
    }
  else
    setsim(set_res,set_beam,set_events);
  dr("hxc3");
  odr("hxgmr3");
  leg = new TLegend(.78,.71,.86,.81);
  leg->AddEntry("hxc3","data");
  leg->AddEntry("hxgmr3","sim");
  leg->Draw();
}

void compX2(Float_t set_res=0, Float_t set_beam=0.607956845, Float_t set_events=1e5)
{
  if((xres==set_res)&&(yres==set_res))
    {
      printf("Resolutions match\n");
      if((sigma_x==set_beam)&&(sigma_y==set_beam))
	printf("...and beam spot sizes match.\n");
      else
	setsim(set_res,set_beam,set_events);  
    }
  else
    setsim(set_res,set_beam,set_events); 
  dr("hxc2");
  odr("hxgmr2");
  leg = new TLegend(.78,.71,.86,.81);
  leg->AddEntry("hxc2","data");
  leg->AddEntry("hxgmr2","sim");
  leg->Draw();
}
void compY1(Float_t set_res=0, Float_t set_beam=0.607956845, Float_t set_events=1e5)
{
  if((xres==set_res)&&(yres==set_res))
    {
      printf("Resolutions match\n");
      if((sigma_x==set_beam)&&(sigma_y==set_beam))
	printf("...and beam spot sizes match.\n");
      else
	setsim(set_res,set_beam,set_events);  
    }
  else
    setsim(set_res,set_beam,set_events);
  dr("hxc1");
  odr("hxgmr1");
  leg = new TLegend(.78,.71,.86,.81);
  leg->AddEntry("hxc1","data");
  leg->AddEntry("hxgmr1","sim");
  leg->Draw();
}
void compX1(Float_t set_res=0, Float_t set_beam=0.607956845, Float_t set_events=1e5)
{
  if((xres==set_res)&&(yres==set_res))
    {
      printf("Resolutions match\n");
      if((sigma_x==set_beam)&&(sigma_y==set_beam))
	printf("...and beam spot sizes match.\n");
      else
	setsim(set_res,set_beam,set_events);  
    }
  else
    setsim(set_res,set_beam,set_events);
  dr("hxc0");
  odr("hxgmr0");
  leg = new TLegend(.78,.71,.86,.81);
  leg->AddEntry("hxc0","data");
  leg->AddEntry("hxgmr0","sim");
  leg->Draw();
}
void compX(Float_t set_res=0, Float_t set_beam=0.607956845, Float_t set_events=1e5, Float_t set_weight=1)
{
  if((xres==set_res)&&(yres==set_res))
    {
      printf("Resolutions match\n");
      if((sigma_x==set_beam)&&(sigma_y==set_beam))
	printf("...and beam spot sizes match.\n");
      else
	setsim(set_res,set_beam,set_events,set_weight);  
    }
  else
    setsim(set_res,set_beam,set_events,set_weight);  
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();      
  cFit->Clear();
  cFit->Divide(1,2);
  cFit->cd(1);
  hxc0->Draw();
  odr("hxgmr0");
  leg = new TLegend(.78,.71,.86,.81);
  leg->AddEntry("hxc0","data");
  leg->AddEntry("hxgmr0","sim");
  leg->Draw();
  cFit->cd(2);
  hxc2->Draw();
  odr("hxgmr2");
  leg = new TLegend(.78,.71,.86,.81);
  leg->AddEntry("hxc2","data");
  leg->AddEntry("hxgmr2","sim");
  leg->Draw();
}
void compY(Float_t set_res=0, Float_t set_beam=0.607956845, Float_t set_events=1e5, Float_t set_weight=1)
{
  if((xres==set_res)&&(yres==set_res))
    {
      printf("Resolutions match\n");
      if((sigma_x==set_beam)&&(sigma_y==set_beam))
	printf(" ...and beam spot sizes match.\n");
      else
	setsim(set_res,set_beam,set_events,set_weight);  
    }
  else
    setsim(set_res,set_beam,set_events,set_weight);  
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();      
  cFit->Clear();
  cFit->Divide(1,2);
  cFit->cd(1);
  hxc1->Draw();
  odr("hxgmr1");
  leg = new TLegend(.78,.71,.86,.81);
  leg->AddEntry("hxc1","data");
  leg->AddEntry("hxgmr1","sim");
  leg->Draw();
  cFit->cd(2);
  hxc3->Draw();
  odr("hxgmr3");
  leg = new TLegend(.78,.71,.86,.81);
  leg->AddEntry("hxc3","data");
  leg->AddEntry("hxgmr3","sim");
  leg->Draw();
}

void show_angles(Float_t set_theta_center=30,Bool_t set_rutherford=kFALSE,
		 Bool_t doset=kFALSE)
{
  theta_center=set_theta_center;
  iprint=kFALSE;
  Float_t set_theta_min=0,set_theta_max=0,set_phi_max=0;
  printf("Calculating angle ranges...\n");
  rotate_axis(z_mask,TMath::ATan((99.5/2)/z_mask)*TMath::RadToDeg(),0);
  printf(" The minimum polar angle is %7.2f\n",theta);  
  set_theta_min=theta;
  rotate_axis(z_mask,TMath::ATan((99.5/2)/z_mask)*TMath::RadToDeg(),180);
  printf(" The maximum polar angle is %7.2f\n",theta);
  set_theta_max=theta;
  // rotate_axis(z_mask,TMath::ATan((99.5/2)/z_mask)*TMath::RadToDeg(),90);
  printf(" (phi = %7.2f",phi);
  if(phi>140)
    set_phi_max=180-phi;
  else
    set_phi_max=phi/2;
  set_phi_max=TMath::RadToDeg()*TMath::ASin((99.5/2)/(X+(99.5/2)));
  printf(" set_phi_max=%7.2f)\n",set_phi_max);
  printf(" The aximuthal angle range is %7.2f (%7.2f to %7.2f)\n",set_phi_max*2,-set_phi_max,set_phi_max);

 if(doset)
   setangles(set_theta_min,set_theta_max,-set_phi_max,set_phi_max,set_theta_center,set_rutherford);
}

