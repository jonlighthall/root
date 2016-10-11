/*-------------------Emacs PostScript "pretty-print" page width (97 columns)-------------------*/
/* Program: load_and_plot.cc (formerly load_psd.cc)
 *       Created  by Jon Lighthall
 * Purpose: 
 *       A set of macros for loading and analyzing data from the:
 *          + 2006-2007 HELIOS PSD characterization tests at  WMU
 *          + August 2008 29Si(d,p) HELIOS commissioning experiment (including simulations)
 *          + 12B(d,p) experiment
 *       ...and plots figures for the HELIOS NIM paper and the Lighthall doctoral thesis.
 * Requires:
 *       fit.cc (most)
 *       psd_analyze.cc (some)
 */
Bool_t bOrig=0;//Turn on/off settings for 1.91T
Bool_t bw=0;//Turn on/off black and white figures
Bool_t bThesis=1;//Turn on/off changes between thesis and NIM paper
Bool_t roman=1;
Bool_t bModern=1;
Int_t pal=11;
Float_t B0=1.91585;   
Bool_t bslope=1;
Bool_t bsimQ=kTRUE;
TFile *_filename0=0;
TFile *_filename1=0;
TFile *_filename2=0;
TFile *_filename3=0;
TFile *_filename4=0;
TFile *_filename5=0;
TFile *_filename6=0;
TFile *_filename7=0;
TFile *_filename8=0;
TFile *_filename9=0;
TFile *_filename10=0;
TFile *_filename11=0;
Int_t reaction=0;

TF1 *f2,*f3;//for plotlabangle()
TF1 *f9,*f10;//for func()
TF1 *f1,*f4,*f5,*f6;//for plotbore()
TF1 *fEDiff1,*fEDiff2,*fXFXN;

Float_t pi=4.0*atan(1.0); 
//Float_t e=1.602176487E-19;//NIST value
//Float_t MeV=e*1E6; //J/MeV

UInt_t width=600;
Float_t ratio1=3./4.;
Float_t ratio2=9./16.;
Float_t ratio3=11./8.5;

Float_t minE=0;
Float_t maxE=10;
Float_t plot_minZ=-880;
Float_t plot_maxZ=-80;
TString ytitle="E (MeV)";
TString xtitle="z (mm)";
Float_t ang=18.5;
Float_t ang_fan=66.5;
Float_t zmax=0;
Float_t Or=0;//origin (need to be declared globally)
Float_t aspan=342.4;
Float_t vlines[3]={-129.25,-367.25,-513.6};//"in sort" position of detectors
Float_t bshift[3]={35,25,22};//corrections for B-field
Float_t edges[3][6]={{-436.67,-377.86,-319.1,-260.47,-202.11,-144.75},
		     {-684.67,-625.86,-567.1,-508.47,-450.11,-392.75},
		     {-834.02,-775.21,-716.45,-657.82,-599.46,-542.1}};
TString fname;

/* Load data-----------------------
 * 
 */

//-----------------------------------------------------------------------------------------------
//Macros specific to the 28Si(d,p) experiment:---------------------------------------------------
void loadfiles(Char_t *histin="_cut")
{//the scope of these files is limited to this script, so it isn't useful.
  hname="100";
  hname+=histin; 
  hname+=".root"; 
  printf("Input File is %s\n",hname.Data());
  _filename0 = new TFile(hname.Data());
  hname="350";
  hname+=histin; 
  hname+=".root"; 
  printf("Input File is %s\n",hname.Data());
  _filename1 = new TFile(hname.Data());
  hname="500";
  hname+=histin; 
  hname+=".root"; 
  printf("Input File is %s\n",hname.Data());
  _filename2 = new TFile(hname.Data());
}

void merge(Char_t *histin="hEZg",Bool_t DoRebin=1)
{//merges histograms from three files.  Adapted from addpositions.C
 // if(!(gROOT->FindObject("home")))
 //    gROOT->ProcessLine("TDirectory *home=gDirectory"); 

  loadfiles();
  _filename0->cd();
  hname=histin;
  TH2F * hInput=(TH2F *) gROOT->FindObject(hname.Data());
  hname+="_100"; 
  if(DoRebin)
    hInput->Rebin2D(hInput->GetNbinsX()/512,hInput->GetNbinsY()/512);
  hInput->Clone(hname.Data());
  hOutput=(TH2F *) gROOT->FindObject(hname.Data());
  hOutput->SetDirectory(home);
  
  _filename1->cd();
  hname=histin;
  TH2F *   hInput=(TH2F *) gROOT->FindObject(hname.Data());
  hname+="_350"; 
  if(DoRebin)
    hInput->Rebin2D(hInput->GetNbinsX()/512,hInput->GetNbinsY()/512);
  hInput->Clone(hname.Data());
  hOutput=(TH2F *) gROOT->FindObject(hname.Data());
  hOutput->SetDirectory(home);
  
  _filename2->cd();
  hname=histin;
  TH2F * hInput=(TH2F *) gROOT->FindObject(hname.Data());
  hname+="_500"; 
  if(DoRebin)
    hInput->Rebin2D(hInput->GetNbinsX()/512,hInput->GetNbinsY()/512);
  hInput->Clone(hname.Data());
  hOutput=(TH2F *) gROOT->FindObject(hname.Data());
  hOutput->SetDirectory(home);
  
  home->cd();
  hname=histin;
  hname+="_100"; 
  TH2F * hInput=(TH2F *) gROOT->FindObject(hname.Data());
  hname=histin;
  hname+="_all"; 
  hInput->Clone(hname.Data());
  hOutput=(TH2F *) gROOT->FindObject(hname.Data());//added
  
  hname=histin;//omitted
  add2(hname+"_100",hname+"_350",hname+"_all");
  add2(hname+"_500",hname+"_all",hname+"_all");
}

void merge2(Char_t *histin="hEZg")
{//same as merge(), but only combines two files
  _filename0->cd();
  hname=histin;
  TH2F * hInput=(TH2F *) gROOT->FindObject(hname.Data());
  hname+="_100";
  hInput->Rebin2D(hInput->GetNbinsX()/512,hInput->GetNbinsY()/512);
  hInput->Clone(hname.Data());
  hOutput=(TH2F *) gROOT->FindObject(hname.Data());
  hOutput->SetDirectory(home);

  _filename2->cd();
  hname=histin;
  TH2F * hInput=(TH2F *) gROOT->FindObject(hname.Data());
  hname+="_500"; 
  hInput->Rebin2D(hInput->GetNbinsX()/512,hInput->GetNbinsY()/512);
  hInput->Clone(hname.Data());
  hOutput=(TH2F *) gROOT->FindObject(hname.Data());
  hOutput->SetDirectory(home);
  
  home->cd();
  hname=histin;
  hname+="_100"; 
  TH2F * hInput=(TH2F *) gROOT->FindObject(hname.Data());
  hname=histin;
  hname+="_all"; 
  hInput->Clone(hname.Data());
  hname=histin;
  add2(hname+"_100",hname+"_500",hname+"_all");
}

void test(Int_t dz=12, Int_t pos=2,Char_t *names="_lineB.txt")
{//tests offset parameters (dz) against calculated E vs. z 
  loadfiles();
  
  switch(pos){
  case 0:
    _filename0->cd();
    TH2F * hInput=(TH2F *) gROOT->FindObject("hEZg");
    if((hInput->GetNbinsY())>512)
      hInput->Rebin2D(hInput->GetNbinsX()/512,hInput->GetNbinsY()/512);
    shiftx2("hEZg",dz);
    dr("hEZg_shift",vlines[0]-352.4+dz,vlines[0]+10+dz);
    cFit->SetLogz();  
    plotvlines(vlines[0]+dz);
    break;
  case 1:
    _filename1->cd();
    TH2F * hInput=(TH2F *) gROOT->FindObject("hEZg");
    if((hInput->GetNbinsY())>512)
      hInput->Rebin2D(hInput->GetNbinsX()/512,hInput->GetNbinsY()/512);
    shiftx2("hEZg",dz);
    dr("hEZg_shift",vlines[1]-352.4+dz,vlines[1]+10+dz);
    cFit->SetLogz();
    plotvlines(0,vlines[1]+dz);
    plotlabangle(ang_fan);
    plotlabangle(ang);
    break;
  case 2:
    _filename2->cd();
    TH2F * hInput=(TH2F *) gROOT->FindObject("hEZg");
    if((hInput->GetNbinsY())>512)
      hInput->Rebin2D(hInput->GetNbinsX()/512,hInput->GetNbinsY()/512);
    shiftx2("hEZg",dz);
    dr("hEZg_shift",vlines[2]-352.4+dz,vlines[2]+10+dz);
    cFit->SetLogz();
    plotvlines(0,0,vlines[2]+dz);
    break;
  default:
    printf("Position %d not recognized!\n",pos);
    break; 
  }
  plotlines(names,7);
}

//HELIOS-specific plotting utilities

void plotlabangle(Float_t angle=20, Float_t origin=0, Int_t reset=0, Float_t amu=1.007, Int_t q=1, Float_t B=1.91585, Int_t orbits=1)
{
  Float_t m=amu*1.661E-27; //Ejectile mass in kg
  //  printf("Mass is %g kg.\n",m);
  Float_t Tcyc=2*pi*m/(q*e*B)*orbits; //Cyclotron period in seconds
  // printf("2 pi M is %g.\n",2*pi*m);
  // printf("q*B is %g, e = %g, qeB=%g.\n",q*B,e,q*e*B);
  printf("Cyclotron period is %.2f ns.\n",Tcyc*pow(10,9));
  printf("Angle is %2.1f deg.  Origin is %2.2f mm.\n",angle,origin);
  if(reset)
    if ((TF1 *) gROOT->FindObject("f3")) {
      //cFit->cd();
      f3->Delete();
      printf("Function \"f3\" already exists.  Deleting...\n");
    }
        
  Or=origin;
  Float_t Zoff(Int_t x)
  {
    return (x-Or);
  }

  // f3 =new TF1("f3","([0]/2.)*pow((Zoff(x)/1000.)*(1./[1])*(1./cos([2])),2.)*(1/[3])",-1200,1200);
  f3 =new TF1("f3","([0]/2.)*pow(((x-[4])/1000.)*(1./[1])*(1./cos([2])),2.)*(1/[3])",-1200,1200);
  f3->SetParameter(4,origin);
  f3->SetParameter(0,m);
  f3->SetParameter(1,Tcyc);//Here "Tcyc" is actually the flight time
  f3->SetParameter(2,(angle/180.0)*pi);
  f3->SetParameter(3,MeV);
  
  f3->SetLineColor(1);
  f3->SetLineStyle(2);
  f3->Draw("same");
}

void func(Float_t r0=0.895*2.54/100./2.)
{//prototype function for plotting kinematic curves
  Float_t m=1.673E-27; //Ejectile mass in kg
  Float_t Tcyc=34.238E-9; //Cyclotron period in seconds
  Float_t Vcm=3.176E7; //Center of mass velocity in m/s
  Float_t V0=5.697E7; //Ejectile velocity in center of mass in m/s
  Int_t min=0,max=0,min2=0,max2=0;
  Float_t slopeEcm=m*Vcm/Tcyc/MeV/1000;

  f4 =new TF1("f4","acos(((x/(1000*[1]))-[2])/[3])",-1200,1200);//CoM angle over full range
  f4->SetParameter(1,Tcyc);
  f4->SetParameter(2,Vcm);
  f4->SetParameter(3,V0);
  for(int i=-1200;i<1200;i++){
    if(f4(i)==f4(i)){//tests if function is valid
      if(i<min)min=i;
      if(i>max)max=i;
    }
  }
  printf("minimum is %d, maximum is %d\n",min,max);
  if ((TH2F *) gROOT->FindObject("h2")) {
    gROOT->FindObject("h2")->Delete();  
    printf("Histogram \"h2\" already exists. Deleting old histogram.\n");
  }

  h2 = new TH2F("h2","Histogram",1000,min,max,100,0.,12);//sets valid range
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
  h2->Draw();
  //should just re-define f4 with new range!
  f5 = new TF1("f5","f4",min,max);//output is center-of-mass angle in radians; works!
  f5->SetLineColor(3);
  f5->SetLineStyle(2);
  f5->Draw("same");
  f6 = new TF1("f6","pi/2+atan(((x/1000)/[1])/([3]*sin(f5)))",min,max);//output is lab angle in radians
  f6->SetParameter(1,Tcyc);
  f6->SetParameter(3,V0);
  f6->SetLineColor(2);
  f6->SetLineStyle(2);
  f6->Draw("same");
  // TF1 *f7 = new TF1("f7","([3]*cos(f5)+[2])*([1]-[0]/([3]*sin(f5)))",min,max);
  TF1 *f7 = new TF1("f7","[3]*cos(f5)+[2]",min,max);//works! correctly produces v_para
  // TF1 *f7 = new TF1("f7","([1]-[0]/([3]*sin(f5)))",min,max);//works! correctly produces time
  printf("ro is %6.4f m\n",r0);
 
  f7->SetParameter(0,r0);
  f7->SetParameter(1,Tcyc);
  f7->SetParameter(2,Vcm);
  f7->SetParameter(3,V0);
  
  TF1 *f8 = new TF1("f8","[3]*sin(f5)",min,max);//calculates V_perp
  f8->SetParameter(3,V0);
  f9 = new TF1("f9","(1/2)*[0]*(f7*f7+f8*f8)/[4]",min,max);//calculates E, plots z for r0=0
  f9->SetParameter(0,m);
  f9->SetParameter(4,MeV);
  f9->SetLineStyle(2);
  f9->Draw("same");
  TF1 *f13 = new TF1("f13","[0]*x+[1]",-1200,1200);//defines f9() over larger range
  f13->SetParameter(0,slopeEcm);
  f13->SetParameter(1,11.676);
  f13->SetLineWidth(2);
  f13->SetLineStyle(9);
  //  TF1 *f10 = new TF1("f10","([0]/tan(f6))*1000",min,max);//offset from z
  f10 = new TF1("f10","([0]/tan(f6))*1000",min,max);//offset from z. works!

  f10->SetParameter(0,r0);
  f10->SetLineColor(3);
  f10->Draw("same");

  Float_t Zoff(Int_t x)
  {
    return x-f10(x);
    //return x-10;
  }

  Float_t Zpara(Int_t x)
  {
    return f9(Zoff(x));
  }

  TF1 *f11 = new TF1("f11","Zpara(x)",min,max);//-834
  
  for(int i=min;i<max;i++){
    if(f11(i)==f11(i)){//tests if function is valid
      if(i<min2)min2=i;
      if(i>max2)max2=i;
    }
  }
  printf("new min is %d, new max is %d\n",min2,max2);
  f11 = new TF1("f11","Zpara(x)",min2,max2);
  f11->SetLineColor(9);
  f11->SetLineStyle(2);
  f11->SetLineWidth(1);
  f11->Draw("same");

  TF1 *f12 = new TF1("f12","Zpara(x)",min2,0);
  f12->SetLineColor(2);
  f12->SetLineStyle(1);
  f12->Draw("same");
}

Float_t func4(Float_t x=0)
{//required for calling f4() at specific X. - moved to outside plotbore()
  return(f4(x));
} 

void plotbore(Float_t radius=36.4*2.54/100/2,Int_t amu=1.007, Int_t q=1, Float_t B=1.91585)
{//plots energy cut-off for a given bore radius
  Float_t m=amu*1.661E-27; //Ejectile mass in kg
  Float_t Tcyc=2*pi*m/(q*e*B); //Cyclotron period in seconds
  Float_t Vcm=3.180E7; //Center of mass velocity in m/s
  Float_t V0=5.701E7; //Ejectile velocity in center of mass in m/s
  switch (reaction){
  case 1:
    //values for d(B11,p) at 81.0 MeV
    Vcm=3.185E7; //Center of mass velocity in m/s
    V0=4.915E7; //Ejectile velocity in center of mass in m/s
    break;
  case 2:
    //values for d(B12,p) at 75.0 MeV
    Vcm=2.972E7; //Center of mass velocity in m/s
    V0=4.884E7; //Ejectile velocity in center of mass in m/s
    break;
  default:
    break;
  }
  Int_t min=0,max=0,min2=0,max2=0;
  printf("WARNING: This fuction requires the center-of-mass velocity Vcm and the ejectile velocity v0 to be input mannually in the code!  Default calculation is for 28Si(d,p)!\n");

  f4 =new TF1("f4","acos(((x/(1000*[1]))-[2])/[3])",-1200,1200);//CoM angle over full range
  f4->SetParameter(1,Tcyc);
  f4->SetParameter(2,Vcm);
  f4->SetParameter(3,V0);

  for(int i=-1200;i<1200;i++){
    if(gROOT->GetVersionInt()<52400){
      if(i==0)printf("Loop 1 - Deprecated ROOT Version.\n");
      if(f4==f4){if(i<min)min=i;if(i>max)max=i;}
    }
    else{
      if(i==0)printf("Loop 1 - Modern ROOT Version.\n");
      if(func4(i)==func4(i)){if(i<min)min=i;if(i>max)max=i;}
    }
  }
  printf("Center-of-mass angle is defined over %d,%d\n",min,max);

  f5 = new TF1("f5","f4",min,max);//output is center-of-mass angle in radians; works!
  f5->SetLineColor(3);
  f5->SetLineStyle(2);

  f6 = new TF1("f6","pi/2+atan(((x/1000)/[1])/([3]*sin(f5)))",min,max);//output is lab angle in radians
  f6->SetParameter(1,Tcyc);
  f6->SetParameter(3,V0);
  f6->SetLineColor(2);
  f6->SetLineStyle(2);

  printf("Limiting radius is %3.2f m\n",radius);
  //  Float_t q=1;
  //Float_t B=1.91585;
  Float_t wc=q*e*B/m;//Cyclotron frequency (radians/sec)

  TF1 *f1=new TF1("f1","([0]/2)*pow(([1]/2)*[2]*(1/sin(f6)),2)/[3]",min,max);
  //f1=f1(f6(f5(f4))), so range should be the same as f4.
  f1->FixParameter(0,m);
  f1->FixParameter(1,radius);
  f1->FixParameter(2,wc);
  f1->FixParameter(3,MeV);
  /*
    for(int i=-1000;i<1000;i++){
    if(gROOT->GetVersionInt()<52400){
    if(i==0)printf("Loop 2 - Deprecated ROOT Version.\n");
    if(f1==f1){if(i<min2)min2=i;if(i>max2)max2=i;}
    }  else{
    if(i==0)printf("Loop 2 - Modern ROOT Version.\n");
    if(f1(i)==f1(i)){if(i<min2)min2=i;if(i>max2)max2=i;}
    }
    }
    printf("For bore function, min is %d, max is %d\n",min2,max2);  
    f1=new TF1("f1","([0]/2)*pow(([1]/2)*[2]*(1/sin(f6)),2)/[3]",min2,max2);
    
    f1->FixParameter(0,m);
    f1->FixParameter(1,radius);
    f1->FixParameter(2,wc);
    f1->FixParameter(3,MeV);
  */ 
  f1->SetLineStyle(9);
  f1->Draw("same");
}

//General plotting utilities
void plotlines(Char_t *names="_line.txt", Int_t lines=6, Int_t linestyle=1)
{//Plots a set of lines from text files.  File names sequencally numbererd, starting with 0.
  TString fname;
  TString lname;
  for(Int_t i=0;i<lines;i++){
    fname="";
    fname+=i;
    fname+=names;
    lname="g";
    lname+=i;
    printf("File is %s.  Line name is %s: ",fname.Data(),lname.Data());
    drawline(fname.Data(),lname.Data(),0,linestyle);
  }
}

void tilt(Char_t * histin, Float_t tilt=0.000413, Int_t plot=1)
{//Tilts a 2D histogram by a given angle.
  zmax=0;
  if(plot>0) if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();    
  Float_t m=1*1.673E-27; //Ejectile mass in kg
  Float_t Vcm=3.174E7; //Center of mass velocity in m/s
  Float_t Tcyc=34.246E-9; //Cyclotron period in seconds
  Float_t slopeEcm=m*Vcm/Tcyc/MeV/1000;
  
  Float_t slope=slopeEcm+tilt;
  Float_t xmin,xmax; 
  Float_t ymin,ymax; 
  
  Int_t xbin;
  Int_t ybin,entry=0;
  Float_t x,y,z;
  if(gROOT->FindObject(histin))  
    TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  else
    printf("Histogram \"%s\" not recognized\n",histin);

  xbin=hInput->GetXaxis()->GetNbins();
  ybin=hInput->GetYaxis()->GetNbins();
  
  xmax=hInput->GetXaxis()->GetXmax();
  ymax=hInput->GetYaxis()->GetXmax();
  xmin=hInput->GetXaxis()->GetXmin();
  ymin=hInput->GetYaxis()->GetXmin();
  
  ymin=floor(ymin-xmax*slope);
  ymax=ceil(ymax-xmin*slope);

  hname=histin;
  hname+="_tilt"; 
  Char_t *htitle = hInput->GetTitle();
 
  if ((TH2F *) gROOT->FindObject(hname)) {
    gROOT->FindObject(hname)->Delete();  
    //printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
  }
  else
    printf("Output histogram is \"%s\"\n",hname.Data());
  TH2F * hOutput=new  TH2F(hname,htitle,xbin,xmin,xmax,ybin,ymin,ymax);
  //  printf("Output histogram is constructed as:\n TH2F(\"%s\",\"%s\",%d,%1.0f,%1.0f,%d,%1.0f,%1.0f)\n",hname.Data(),htitle,xbin,xmin,xmax,ybin,ymin,ymax);
  hname+="_py"; 
  if ((TH1F *) gROOT->FindObject(hname)) {
    gROOT->FindObject(hname)->Delete();  
    //printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
  }
  else
    printf("Output histogram is \"%s\"\n",hname.Data());
  TH1F * hProj=new  TH1F(hname,htitle,ybin,ymin,ymax);
  
  for(int i=1;i<(xbin+1);i++){
    for(int j=1;j<(ybin+1);j++){
      //Note: The 0 bin contains the underflow, so the loop starts at 1;
      //      and the max_bin+1 contains overflow, so the loop terminates at max_bin+1
      x=hInput->GetXaxis()->GetBinCenter(i);
      y=hInput->GetYaxis()->GetBinCenter(j);
      z=hInput->GetBinContent(i,j);
      //printf("i=%2d, j=%2d, z=%2.0f \n",i,j,z);
      hProj->Fill(y-slope*x,z);
      hOutput->Fill(x,y-slope*x,z);
    }
  }
  for(int j=1;j<(ybin+1);j++){
    z=hProj->GetBinContent(j);
    if(z>zmax)
      zmax=z;//locates largest bin
  }
  
  printf("For slope = %.6f MeV/mm, Maximum peak height = %.0f counts.\n",slope,zmax);
  //hProj->SetAxisRange(0,40000,"Y");  

  if(plot==2){
    cFit->Clear();
    cFit->Divide(1,2);
    cFit->cd(1);
    hOutput->Draw("col");
    cFit->SetLogz();   
    cFit->cd(2);
  }

  if(plot==3){
    cFit->Clear();
    cFit->Divide(1,3);
    cFit->cd(1);
    hInput->Draw("col");
    cFit->cd(2);
    hOutput->Draw("col");
    cFit->SetLogz();   
    cFit->cd(3);
  }

  if(plot>0) hProj->Draw();
  /*
    Float_t mu=12;
    Float_t range=.700; 
    hProj->Fit("gaus","Q","",mu-range/2,mu+range/2);
    Float_t sig=hProj->GetFunction("gaus")->GetParameter(2);
    Float_t fw=2*sqrt(2*log(2))*sig;
    printf("Sigma = %.3f keV, FWHM = %.3f keV\n",sig/1000,fw/1000);
  */
}

Float_t gettilt(Char_t * histin, Float_t tilt)
{
  tilt(histin,tilt);
  return zmax;
}

void findtilt(Char_t * histin="hEZg", Int_t steps=1)
{
  for(int i=(-1*steps);i<steps;steps++;){
    printf("%3d: ",i);
    //    gettilt(histin,i/10000.);
  }
}

//Macros for generating figures the HELIOS NIM paper and Lighthall thesis.--------------

void paperplots()
{
  gStyle->SetOptDate(0);
  if(bOrig){
    B0=1.91585;   
    printf("Detector positions are: %.2f, %.2f, %.2f\n",vlines[0],vlines[1],vlines[2]);
  }
  else{
    B0=2.0020;   
    printf("Detector positions are: %.2f, %.2f, %.2f\n",vlines[0]+bshift[0],vlines[1]+bshift[1],vlines[2]+bshift[2]);
  }
  Bool_t bEx=0;//turn on/off excitation plots

  loadfiles("_cut");
  Float_t pttext=0.06;
  TPaveText *pt;
  TText *text;  
  //1.) Time peaks-----------------------
  _filename0->cd();
  hET21->ProjectionX();
  if(!bOrig){
    slopex("hET21_px",1.045,-0.244);//corrects for B=2.00 field
    TH1F * hOutput1=(TH1F *) gROOT->FindObject("hET21_px_slope"); 
  }
  else
    TH1F * hOutput1=(TH1F *) gROOT->FindObject("hET21_px"); 
  hOutput1->SetDirectory(home);
  hOutput1->SetAxisRange(0,82,"X");
  hOutput1->SetYTitle("Counts per 535 ps");
  hOutput1->GetYaxis()->CenterTitle(1);
  hOutput1->GetYaxis()->SetTitleOffset(1.23);
  hOutput1->SetXTitle("T (ns)");
  hOutput1->GetXaxis()->CenterTitle(1);
  hOutput1->SetStats(kFALSE);
  hOutput1->SetTitle();
  mkCanvas2("cTOF","cTOF");
  cTOF->SetTopMargin(.02);
  cTOF->SetRightMargin(.02);
  cTOF->SetCanvasSize(width,(UInt_t)(width*ratio2));
  hOutput1->Draw();
  pt = new TPaveText(37,12677.06,45,14361.89,"br");
  text = pt->AddText("A/q=1, 1 turn");
  //text = pt->AddText("A/q=1 2 turns");
  pt->SetTextAlign(12);//Middle, Left
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();

  pt = new TPaveText(64,5018.749,72,6754.634,"br");
  text = pt->AddText("A/q=2, 1 turn");
  text = pt->AddText("A/q=1, 2 turns");
  //text = pt->AddText("A/q=2");
  //text = pt->AddText("2 turns");
  pt->SetTextAlign(20);//Middle, Middle
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();

  cTOF->SaveAs("cTOF_21_100.eps");
  
  merge("hEZg");
  cFit->Close();
  home->cd();
  /*
  //2.) Energy vs. Position (all)--------
  TH2F * hOutput2=(TH2F *) gROOT->FindObject("hEZg_all");
  hOutput2->SetAxisRange(plot_minZ,plot_maxZ,"X");
  hOutput2->SetAxisRange(minE,maxE,"Y");
  hOutput2->SetYTitle(ytitle);
  hOutput2->GetYaxis()->CenterTitle(1);
  hOutput2->SetXTitle(xtitle);
  hOutput2->GetXaxis()->CenterTitle(1);
  hOutput2->SetStats(kFALSE);
  hOutput2->SetTitle();
  hOutput2->SetMinimum(8);
  mkCanvas2("cComp","cComp");
  cComp->SetTopMargin(.02);
  cComp->SetRightMargin(.03);  
  cComp->SetCanvasSize(width,(UInt_t)(width*ratio1));
  cComp->SetLogz();  
  hOutput2->Draw("col");
  plotlabangle(ang,0);
  plotbore();

  //drawline("gs_line.txt","gs",0);
  //gs->SetLineWidth(1);
  //gs->SetLineStyle(1);
  //gs->SetLineColor(1);
  //  gs->Draw();
  plotvlines(vlines[0],vlines[1],vlines[2]);
  cComp->SaveAs("cEZg_all.eps");
  */
  //4.)Alternate e vs. z-------------------
  home->cd();
  if(bOrig){
    shiftx2("hEZg_100",8);
    shiftx2("hEZg_350",0);
    shiftx2("hEZg_500",0);
  }
  else{
    shiftx2("hEZg_100",bshift[0]);
    shiftx2("hEZg_350",bshift[1]);
    shiftx2("hEZg_500",bshift[2]);
  }

  add2("hEZg_100_shift","hEZg_500_shift","hEZg_shift");
  add2("hEZg_100_shift","hEZg_350_shift","hEZg_ang");
  add2("hEZg_350_shift","hEZg_shift","hEZg_shift");
  cFit->Close();

  TH2F * hOutput4=(TH2F *) gROOT->FindObject("hEZg_shift");
  hOutput4->SetAxisRange(plot_minZ,plot_maxZ,"X");
  hOutput4->SetAxisRange(minE,maxE,"Y");
  hOutput4->SetYTitle(ytitle);
  hOutput4->GetYaxis()->CenterTitle(1);
  hOutput4->SetXTitle(xtitle);
  hOutput4->GetXaxis()->CenterTitle(1);
  hOutput4->SetStats(kFALSE);
  hOutput4->SetTitle();
  hOutput4->SetMinimum(8);
  mkCanvas2("cShift","cShift");
  cShift->SetTopMargin(.02);
  cShift->SetRightMargin(.03);  
  cShift->SetCanvasSize(width,(UInt_t)(width*ratio1));
  cShift->SetLogz();  
  if(bw){gStyle->SetPalette(pal);hOutput4->SetMaximum(-1111);}
  //  hOutput4->Rebin2D(hOutput4->GetNbinsX()/512,hOutput4->GetNbinsY()/512);
  hOutput4->Draw("col");
  if(bOrig){
    plotlabangle(ang,0,1.007,1,B0);//ang=18.5,B=1.91
    plotbore();
    plotvlines(vlines[0]+8,vlines[1],vlines[2],!bOrig);  
    if(bslope)
      plotlines("_line.txt",7,7);
  }
  else{
    plotlabangle(ang,0,0,1.007,1,B0);
    plotbore(.462,1.007,1,B0);
    plotvlines(vlines[0]+bshift[0],vlines[1]+bshift[1],vlines[2]+bshift[2],1);  
    if(bslope)
      plotlines("_lineB.txt",7,7);
  }
  if(bThesis){
    pt = new TPaveText(-784,3,-784,3,"br");
    text = pt->AddText("A");
    pt->SetTextAlign(22);//Middle, Middle
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();
  }

  if(bw)
    cShift->SaveAs("cShift_bw.eps");
  else
    cShift->SaveAs("cShift.eps");

  if(bEx){
    //3.) Excitation energy(at one position)--------------------
    _filename0->cd();
    TH2F * hInput3=(TH2F *) gROOT->FindObject("hQTheta");
    hInput3->ProjectionY();
    bkgfit2("hQTheta_py");
    slopex("hQTheta_py_out",1,-0.141015);
    dump("hQTheta_py_out_slope","hEx_one_position.dat");
    cFit->Close();
    TH1F * hOutput3=(TH1F *) gROOT->FindObject("hQTheta_py_out_slope"); 
    hOutput3->SetDirectory(home);
    hOutput3->SetAxisRange(-1,10,"X");
    hOutput3->SetXTitle("Excitation Energy (MeV)");
    hOutput3->GetXaxis()->CenterTitle(1);
    hOutput3->SetYTitle("Counts");
    hOutput3->GetYaxis()->CenterTitle(1);
    hOutput3->GetYaxis()->SetTitleOffset(1.25);
    hOutput3->SetStats(kFALSE);
    hOutput3->SetTitle();
    mkCanvas2("cQy","cQy");
    cQy->SetTopMargin(.02);
    cQy->SetRightMargin(.02);
    cQy->SetCanvasSize(width,(UInt_t)(width*ratio2));
    hOutput3->Draw();

    cQy->SaveAs("cExcite_100.eps");
  
  }
  //5.) Individual position E vs. Z-------
  home->cd();
  if(bOrig){
    TH2F * hOutput5=(TH2F *) gROOT->FindObject("hEZg_350");
    hOutput5->SetAxisRange(vlines[1]-352.4,vlines[1]+10,"X");
  }
  else{
    shiftx2("hEZg_350",bshift[1],0);
    TH2F * hOutput5=(TH2F *) gROOT->FindObject("hEZg_350_shift");
    hOutput5->SetAxisRange(vlines[1]-352.4+bshift[1],vlines[1]+12+bshift[1],"X");
  }
  hOutput5->SetAxisRange(minE,maxE,"Y");
  hOutput5->SetYTitle(ytitle);
  hOutput5->GetYaxis()->CenterTitle(1);
  hOutput5->SetXTitle(xtitle);
  hOutput5->GetXaxis()->CenterTitle(1);
  hOutput5->SetStats(kFALSE);
  hOutput5->SetTitle();
  hOutput5->SetMinimum(8);
  mkCanvas2("cIndiv","cIndiv");
  cIndiv->SetTopMargin(.02);
  if(bOrig)
    cIndiv->SetRightMargin(.2);  
  else
    cIndiv->SetRightMargin(.03);  
  cIndiv->SetCanvasSize(width,(UInt_t)(width*ratio1));
  cIndiv->SetLogz();  
  if(bw){gStyle->SetPalette(pal);hOutput5->SetMaximum(-1111);}
  hOutput5->Draw("col");
  if(bOrig){
    plotvlines(0,vlines[1]);
    plotlabangle(ang,0,0,1.007,1,B0);
    plotlines("_line.txt",6,1);
    plotlabangle(ang_fan);
  }
  else{
    plotvlines(0,vlines[1]+bshift[1],0,0);
    plotlabangle(ang,0,0,1.007,1,B0);
    plotlines("_lineB.txt",6,7);
    plotlabangle(ang_fan,0,0,1.007,1,B0);
  }

  if(!bOrig){
    // f3->SetLineStyle(9);//set line style to match plotlabangle (upper bound)
    f3->Draw("same");
    Float_t left=vlines[1]+bshift[1]+1.5;
    Float_t right=left+7.5;
    pt = new TPaveText(left,7.8,right,8.6,"br");
    text = pt->AddText("a");
    pt->SetTextAlign(12);//Middle, Left
    pt->SetFillColor(0); 
    //    pt->SetFillStyle(0);
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    pt = new TPaveText(left,6.6,right,7.4,"br");
    text = pt->AddText("b");
    pt->SetTextAlign(12);//Middle, Left
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();
    
    pt = new TPaveText(left,5.9,right,6.7,"br");
    text = pt->AddText("c");
    pt->SetTextAlign(12);//Middle, Left
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();
    
    pt = new TPaveText(left,4.9,right,5.7,"br");
    text = pt->AddText("d");
    pt->SetTextAlign(12);//Middle, Left
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();
    
    pt = new TPaveText(left,4.3,right,5.0,"br");
    text = pt->AddText("e");
    pt->SetTextAlign(12);//Middle, Left
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();
    
    pt = new TPaveText(left,3.0,right,3.8,"br");
    text = pt->AddText("f");
    pt->SetTextAlign(22);//Middle, Middle
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    pt = new TPaveText(-525,8.8,-525,9.8,"br");
    text = pt->AddText("B");
    pt->SetTextAlign(22);//Middle, Middle
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    pt = new TPaveText(-625,1.5,-625,1.5,"br");
    text = pt->AddText("A");
    pt->SetTextAlign(22);//Middle, Middle
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

  }
  else{
    pt = new TPaveText(-345,7.8,-314,8.6,"br");
    text = pt->AddText("g.s.");
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    pt = new TPaveText(-350.8,6.6,-253,7.4,"br");
    text = pt->AddText("1.27 MeV");
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    pt = new TPaveText(-350.8,5.9,-253,6.7,"br");
    text = pt->AddText("2.03 MeV");
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    pt = new TPaveText(-350.8,4.9,-253,5.7,"br");
    text = pt->AddText("3.07 MeV");
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    pt = new TPaveText(-350.8,4.3,-253,5.0,"br");
    text = pt->AddText("3.62 MeV");
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    pt = new TPaveText(-350.8,3.0,-253,3.8,"br");
    text = pt->AddText("4.90 MeV");
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();
  
    pt = new TPaveText(-565.9,8.8,-517.2,9.8,"br");
    text = pt->AddText("#theta_{lab}=66^{#circ}");
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    pt = new TPaveText(-700.7022,0.8087121,-614.6112,1.667969,"br");
    text = pt->AddText("#theta_{lab}=18^{#circ}");
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();
  }
  if(bw)
    cIndiv->SaveAs("cEZg_350_bw.eps");
  else
    cIndiv->SaveAs("cEZg_350.eps");

  //6. double orbit ---------------------
  _filename3=new TFile("500_new.root");
  hEZg2->Clone();
  hEZg2->SetDirectory(home);
  home->cd();
  if(bOrig){
    TH2F * hOutput6=(TH2F *) gROOT->FindObject("hEZg2");
    hOutput6->SetAxisRange(vlines[2]-352.4,vlines[2]+10,"X");
  }
  else{
    shiftx2("hEZg2",bshift[2],0);
    TH2F * hOutput6=(TH2F *) gROOT->FindObject("hEZg2_shift");
    hOutput6->SetAxisRange(vlines[2]-352.4+bshift[2],vlines[2]+10+bshift[2],"X");
  }
  hOutput6->Rebin2D(hOutput6->GetNbinsX()/512,hOutput6->GetNbinsY()/512);
  hOutput6->SetAxisRange(minE,maxE,"Y");
  hOutput6->SetYTitle(ytitle);
  hOutput6->GetYaxis()->CenterTitle(1);
  hOutput6->SetXTitle(xtitle);
  hOutput6->GetXaxis()->CenterTitle(1);
  hOutput6->SetStats(kFALSE);
  hOutput6->SetTitle();
  if(bOrig){
    hOutput6->SetMinimum(8);
    hOutput6->SetMaximum(50);
  }
  else
    {
    }
  mkCanvas2("cDouble","cDouble");
  cDouble->SetTopMargin(.02);
  cDouble->SetRightMargin(.03);  
  cDouble->SetCanvasSize(width,(UInt_t)(width*ratio1));
  cDouble->SetLogz();  

  if(bw){gStyle->SetPalette(pal);hOutput6->SetMaximum(-1111);}
  hOutput6->Draw("col");

  if(bOrig){
    plotlabangle(ang_fan,0,0,1.007,1,B0,2);
    plotvlines(0,0,vlines[2]);
  }
  else{
    plotvlines(0,0,vlines[2]+bshift[2],0);
    plotlines("_line2B.txt",7,7);

    Float_t left=vlines[2]+bshift[2]+2;
    Float_t right=left;
    pt = new TPaveText(left,5.4,right,5.9,"br");
    text = pt->AddText("e");
    pt->SetTextAlign(12);//Middle, Left
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();
    
    pt = new TPaveText(left+3,4.1,right,4.8,"br");
    text = pt->AddText("f");
    pt->SetTextAlign(12);//Middle, Left
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    pt = new TPaveText(left-1,2.6,right,3.4,"br");
    text = pt->AddText("g");
    pt->SetTextAlign(12);//Middle, Left
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    plotlabangle(ang,0,0,1.007,1,B0,2);
    plotlabangle(ang_fan,0,0,1.007,1,B0,2);

    pt = new TPaveText(-818,5.8,-818,5.8,"br");
    text = pt->AddText("B");
    pt->SetTextAlign(22);//Middle, Middle
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    pt = new TPaveText(-580,1,-580,1,"br");
    text = pt->AddText("A");
    pt->SetTextAlign(22);//Middle, Middle
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    pt->Draw();

    
  }
  if(bw)
    cDouble->SaveAs("cDouble_bw.eps");
  else
    cDouble->SaveAs("cDouble.eps");

  //10.) Efficiency Plot----------------------------------
  home->cd();
  //  shiftx2("hEZg_100",bshift[0],0);
  TH2F * hOutput10=(TH2F *) gROOT->FindObject("hEZg_100_shift");

  hOutput10->SetAxisRange(vlines[0]-352.4+bshift[0],vlines[0]+12+bshift[0],"X");

  if((hOutput10->GetNbinsY())>512)
    hOutput10->Rebin2D(hOutput10->GetNbinsX()/512,hOutput10->GetNbinsY()/512);
  
  hOutput10->SetAxisRange(minE,maxE,"Y");
  hOutput10->SetYTitle(ytitle);
  hOutput10->GetYaxis()->CenterTitle(1);
  hOutput10->SetXTitle(xtitle);
  hOutput10->GetXaxis()->CenterTitle(1);
  hOutput10->SetStats(kFALSE);
  hOutput10->SetTitle();
  hOutput10->SetMinimum(10);
  hOutput10->SetMaximum(350);
  mkCanvas2("cEffic","cEffic");
  cEffic->SetTopMargin(.02);
  if(bOrig)
    cEffic->SetRightMargin(.2);  
  else
    cEffic->SetRightMargin(.03);  
  cEffic->SetCanvasSize(width,(UInt_t)(width*ratio1));
  cEffic->SetLogz();  
  if(bw){gStyle->SetPalette(pal);hOutput10->SetMaximum(-1111);}
  hOutput10->Draw("col");

  //plotvlines(0,vlines[0]+bshift[0],0,0);
  plotlabangle(ang,0,0,1.007,1,B0);
  //  plotlines("_lineB.txt",6,7);
  plotcuts(0,0.21,6,30);
  cEffic->SaveAs("cEZg_100_effic.eps");
  
  // Excitation Energy Plots------------------------------
  if(bEx){ 
    _filename0->cd();
    plotallpjy("hEcX");

    //7.) 1 detector per position excitation energy spectrum
    add2("hEcX19_py","hEcX2_py","hEc_100");
    add2("hEcX21_py","hEc_100","hEc_100");
    add2("hEcX22_py","hEc_100","hEc_100");
    add2("hEcX11_py","hEc_100","hEc_100");
    add2("hEcX24_py","hEc_100","hEc_100");
    hEc_100->SetDirectory(home);
    hEcX19_py->SetDirectory(home);
    home->cd();
    slopex("hEc_100",-0.969865,11.8159);
    bkgfit2("hEc_100_slope");
    //  cFit->Close();
    TH1F * hOutput7=(TH1F *) gROOT->FindObject("hEc_100_slope_out"); 
    dump("hEc_100_slope_out","hEx_six_detectors.dat");
    //  hOutput7->SetDirectory(home);
    hOutput7->SetAxisRange(-1,10,"X");
    hOutput7->SetXTitle("Excitation Energy (MeV)");
    hOutput7->GetXaxis()->CenterTitle(1);
    hOutput7->SetYTitle("Counts");
    hOutput7->GetYaxis()->CenterTitle(1);
    hOutput7->GetYaxis()->SetTitleOffset(1.25);
    hOutput7->SetStats(kFALSE);
    hOutput7->SetTitle();
    mkCanvas2("cEx_100","cEx_100");
    cEx_100->SetTopMargin(.02);
    cEx_100->SetRightMargin(.02);
    cEx_100->SetCanvasSize(width,(UInt_t)(width*ratio2));
    hOutput7->Draw();
    cEx_100->SaveAs("cExcite_100_6det.eps");

    //9). Single-detector excitation spectrum
    slopex("hEcX19_py",-0.978476,11.7848);
    bkgfit2("hEcX19_py_slope");
    TH1F * hOutput9=(TH1F *) gROOT->FindObject("hEcX19_py_slope_out"); 
    hOutput9->SetAxisRange(-1,10,"X");
    hOutput9->SetXTitle("Excitation Energy (MeV)");
    hOutput9->GetXaxis()->CenterTitle(1);
    hOutput9->SetYTitle("Counts");
    hOutput9->GetYaxis()->CenterTitle(1);
    hOutput9->GetYaxis()->SetTitleOffset(1.25);
    hOutput9->SetStats(kFALSE);
    hOutput9->SetTitle();
    hOutput9->Clone("hEx_19");
    dump("hEx_19","hEx_one_detector.dat");
    hEx_19->SetDirectory(home);
    mkCanvas2("cEcX19","cEcX19");
    cEcX19->SetTopMargin(.02);
    cEcX19->SetRightMargin(.02);
    cEcX19->SetCanvasSize(width,(UInt_t)(width*ratio2));
    hOutput9->Draw();

    pt = new TPaveText(.877,.772,.965,.934,"NDC");
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    text = pt->AddText("a");
    pt->Draw();

    cEcX19->SaveAs("cEcX19.eps");

    //8.)global excitation spectrum--------
    buildexcite();
    bkgfit2("hTilt_all_slope",2,4,.005);
    hTilt_all_slope_out->Smooth(1);
    dump("hTilt_all_slope_out","hEx_all_positions.dat");
    cFit->Close();
    TH1F * hOutput8=(TH1F *) gROOT->FindObject("hTilt_all_slope_out"); 
    hOutput8->SetAxisRange(-1,10,"X");
    hOutput8->SetXTitle("Excitation Energy (MeV)");
    hOutput8->GetXaxis()->CenterTitle(1);
    hOutput8->SetYTitle("Counts");
    hOutput8->GetYaxis()->CenterTitle(1);
    hOutput8->GetYaxis()->SetTitleOffset(1.25);
    hOutput8->SetStats(kFALSE);
    hOutput8->SetTitle();
    mkCanvas2("cEx_all","cEx_all");
    cEx_all->SetTopMargin(.02);
    cEx_all->SetRightMargin(.02);
    cEx_all->SetCanvasSize(width,(UInt_t)(width*ratio2));
    hOutput8->Draw();

    pt = new TPaveText(.877,.772,.965,.934,"NDC");
    pt->SetFillColor(0); 
    if(bModern)pt->SetShadowColor(0);
    pt->SetLineColor(0);  
    pt->SetTextSize(pttext); 
    text = pt->AddText("b");  
    pt->Draw();

    cEx_all->SaveAs("cExcite_all.eps");
    dr("hTilt_all_slope_out");
    peakfit("hTilt_all_slope_out","excite11.lst",2,1.2,.006);
    hTilt_all_slope_out->SetAxisRange(-0.35,8.7,"X");
    peakfit("hTilt_all_slope_out","excite11.lst",2,1.2,.006);
  }//end if(bEx)
}//end paperplots()

void buildexcite(Bool_t bin1=1,Bool_t limit1=0,Bool_t bin2=1,Bool_t limit2=0)
{//called in paperplots()
  Float_t min=-23+12,max=0+12;
  Float_t min1=min,max1=max;
  Float_t min2=min,max2=max;

  if(!limit1){
    min1=0;
    max1=0;
  }

  if(!limit2){
    min2=0;
    max2=0;
  }
  //loadfiles and perform first calibration
  loadfiles("_cut");
  _filename0->cd();
  tilt("hEZg",0.000413,0);
  slopex("hEZg_tilt_py",-0.973357,11.9907,bin1,min,max);   // 996 bins
  hEZg_tilt_py_slope->Clone("hTilt_100");
  hTilt_100->SetDirectory(home);
  _filename1->cd();
  tilt("hEZg",0.000780,0);
  slopex("hEZg_tilt_py",-0.97673,12.0245,bin1,min,max);    //1000 bins
  hEZg_tilt_py_slope->Clone("hTilt_350");
  hTilt_350->SetDirectory(home);
  _filename2->cd();
  tilt("hEZg",0.001193,0);
  slopex("hEZg_tilt_py",-1.017,12.3081,bin1,min,max);      //1041 bins
  hEZg_tilt_py_slope->Clone("hTilt_500");
  hTilt_500->SetDirectory(home);

  home->cd();  
  //perform second calibratoion depending on first
  if(bin1){
    if(limit1){//(1,1)
      slopex("hTilt_100",0.997653,0.00385078,bin2,min,max);
      slopex("hTilt_350",1.00166,-0.00108788,bin2,min,max);
      slopex("hTilt_500",0.977138,0.0546509,bin2,min,max);
    }
    else{//(1,0)
      slopex("hTilt_100",0.997555,0.00504853,bin2,min,max);// 993 bins
      slopex("hTilt_350",1.001,-0.00673543,bin2,min,max);  //1001 bins
      slopex("hTilt_500",0.976941,0.0574323,bin2,min,max); //1016 bins
    }
  }
  else{//(0,0)or(0,1)
    slopex("hTilt_100",0.99435,0.0188845,bin2,min,max);
    slopex("hTilt_350",1.00084,-0.00219133,bin2,min,max);
    slopex("hTilt_500",0.978779,0.0490576,bin2,min,max);
  }

  //dump to file and load into one histogram
  hTilt_100->Clone("hTilt_all");
  hTilt_all->Reset();
  dump("hTilt_100","hTilt_100.txt");
  dump("hTilt_350","hTilt_350.txt");
  dump("hTilt_500","hTilt_500.txt");
  cFit->cd(1);
  fill1("hTilt_100.txt","hTilt_all",0);
  fill1("hTilt_350.txt","hTilt_all",0);
  fill1("hTilt_500.txt","hTilt_all",0);
  dump("hTilt_all","hTilt_all.txt"); 
  hTilt_100_slope->Clone("hTilt_all_slope");
  hTilt_all_slope->Reset();
  dump("hTilt_100_slope","hTilt_100_slope.txt");
  dump("hTilt_350_slope","hTilt_350_slope.txt");
  dump("hTilt_500_slope","hTilt_500_slope.txt");
  cFit->cd(2);
  fill1("hTilt_100_slope.txt","hTilt_all_slope",0);
  fill1("hTilt_350_slope.txt","hTilt_all_slope",0);
  fill1("hTilt_500_slope.txt","hTilt_all_slope",0);

  dump("hTilt_all_slope","hTilt_all_slope.txt");
}

void test3(Float_t min=-2, Float_t max=10,Int_t bins=1024,Char_t * filename="hTilt_all.txt")
{
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();    
  hname="hnew";
  if ((TH1F *) gROOT->FindObject(hname.Data())) {
    gROOT->FindObject(hname)->Delete();  
    printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
  } 
  
  hnew = new TH1F(hname.Data(),"Small Histogram",bins,min,max);
  fill1(filename,"hnew",1);

}

void plotsim()
{
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  if(!((TH1F *) gROOT->FindObject("hError1"))){
    if(gROOT->GetVersionInt()<52400){
      printf("Deprecated (online) ROOT version %s.\n",gROOT->GetVersion());
      gROOT->LoadMacro("/net/helios/Si28/offline/lighthall/simulations/PlotSim.cxx"); 
    }
    else{
      gROOT->LoadMacro("C:/root/macros/PlotSim.cxx");
      gROOT->LoadMacro("PlotSim.cxx");
      gROOT->LoadMacro("simulations/PlotSim.cxx");
    }
    printf("PlotSims.cxx loaded.\n");
    loadparameters();
    makehists();
  }
  else{
    printf("PlotSims.cxx already loaded.\n");    
  }
}
void loadsim(Char_t *set="")
{
  plotsim();
  if(!bsimQ)  printf("Offset is set to %4.0f mm.\n",hp[14]);
  TH2F * hInput=(TH2F *) gROOT->FindObject("hEZg");
  //if(!(Bool_t)((TTree *)gROOT->FindObject("tree1"))){
  if(!(gROOT->FindObject("tree1")))
    printf("\"tree1\" not found.  Calling filltree...");
  else
    printf("\"tree1\" found\n");
  hname="simulations/094";
  hname+=set;
  fname=hname;
  hname+=".dat"; 
  fname+="_param.dat"; 

  if(bsimQ){
    printf("bsimQ is true!\n");
    filltree(hname,"tree1",fname);
    printf("Offset is set to %4.0f mm.\n",hp[14]);
    TH2F * hOutput1=(TH2F *) gROOT->FindObject("hEZg");
  }
  else{
    printf("bsimQ is false!\n");
    filltree(hname,"tree1");
    shiftx2("hEZg",-hp[14]+vlines[0]+bshift[0]);
    TH2F * hOutput1=(TH2F *) gROOT->FindObject("hEZg_shift");
  }

  hOutput1->Clone("hEZg_100");
  hInput->Reset();
  
  //(if!((TTree *)gROOT->FindObject("tree2")))
  hname="simulations/342";
  hname+=set;
  fname=hname;
  hname+=".dat"; 
  fname+="_param.dat"; 
    
  if(bsimQ){
    filltree(hname,"tree2",fname);
    printf("Offset is set to %4.0f mm.\n",hp[14]);
    TH2F * hOutput2=(TH2F *) gROOT->FindObject("hEZg");
  }
  else{
    filltree(hname,"tree2");
    shiftx2("hEZg",-hp[14]+vlines[1]+bshift[1]);
    TH2F * hOutput2=(TH2F *) gROOT->FindObject("hEZg_shift");
  }

  hOutput2->Clone("hEZg_350");
  hInput->Reset();
  
  //  (if!((TTree *)gROOT->FindObject("tree3")))
  hname="simulations/492";
  hname+=set;
  fname=hname;
  hname+=".dat"; 
  fname+="_param.dat"; 
  if(bsimQ){
    filltree(hname,"tree3",fname);
    printf("Offset is set to %4.0f mm.\n",hp[14]);
    TH2F * hOutput3=(TH2F *) gROOT->FindObject("hEZg");
  }
  else{
    filltree(hname,"tree3");
    shiftx2("hEZg",-hp[14]+vlines[2]+bshift[2]);
    TH2F * hOutput3=(TH2F *) gROOT->FindObject("hEZg_shift");
  }
  
  hOutput3->Clone("hEZg_500");
  hInput->Reset();
  
  add2("hEZg_100","hEZg_350");
  add2("hEZg_100_copy","hEZg_500","hEZg");
}

void plotsimdata(Char_t* hname="test.dat")
{//used for PlotSim.cxx
  plotsim();
  filltree(hname,"tree");
  analyze("tree");
}

void simpaperplots()
{//To be run from the /simulations directory (on Mac, parent of).
  gStyle->SetOptDate(0);
  if(bOrig){
    B0=1.91585;   
    printf("Detector positions are: %.2f, %.2f, %.2f\n",vlines[0],vlines[1],vlines[2]);
  }
  else{
    B0=2.0020;   
    printf("Detector positions are: %.2f, %.2f, %.2f\n",vlines[0]+bshift[0],vlines[1]+bshift[1],vlines[2]+bshift[2]);
  }
  if(gROOT->GetVersionInt()<52400){
    printf("Deprecated (online) ROOT version.\n");
    bModern=0;
  }
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  if(!((TTree *)gROOT->FindObject("tree1")))
    loadsim();
  else
    printf("Simulations already loaded.\n");
 
  TH2F * hInput=(TH2F *) gROOT->FindObject("hEZg");
  dr(hInput->GetName());
  hInput->SetAxisRange(plot_minZ,plot_maxZ,"X");
  hInput->SetAxisRange(minE,maxE,"Y");
  hInput->SetYTitle(ytitle);
  hInput->GetYaxis()->CenterTitle(1);
  hInput->SetXTitle(xtitle);
  hInput->GetXaxis()->CenterTitle(1);
  hInput->SetStats(kFALSE);
  hInput->SetTitle();

  mkCanvas2("cSim","cSim");
  cSim->SetCanvasSize(width,(UInt_t)(width*ratio1));
  cSim->SetTopMargin(.02);
  cSim->SetRightMargin(.03);
  if(bw){gStyle->SetPalette(pal);hInput->SetMaximum(7);}
  hInput->Draw("col");
  if(bOrig){
    plotlabangle(ang,0,1.007,1,B0);//ang=18.5,B=1.91
    plotbore();
    plotvlines(vlines[0]+8,vlines[1],vlines[2],!bOrig);  
  }
  else{
    plotlabangle(ang,0,0,1.007,1,B0);
    plotbore(.462,1.007,1,B0);
    plotvlines(vlines[0]+bshift[0],vlines[1]+bshift[1],vlines[2]+bshift[2],!bOrig);  
  }  
  
  /*
    drawline("gs_line.txt","gline",0);
    gline->SetLineWidth(1);
    gline->SetLineStyle(1);
    gline->SetLineColor(1);
    gline->Draw();
  */ 
  if(bw)
    cSim->SaveAs("cSim_bw.eps");
  else
    cSim->SaveAs("cSim.eps");
  cFit->Close();
}

void setvlines(Float_t ymin=-999, Float_t ymax=999)
{
  minE=ymin;
  maxE=ymax;
  printf("Vertical extent of plotted lines has been set to (%f,%f).\n",minE,maxE);
}

void plotvlines(Float_t edge1=0., Float_t edge2=0., Float_t edge3=0.,Bool_t span=0, Int_t color=1)
{//shows detector coverage area(s)
  Float_t z;
  Float_t gap=.625;
  if(span)
    printf("Sensor span is %f\n",aspan);
  if(edge1!=0){
    z=edge1;
    TLine *line = new TLine(z,minE,z,maxE);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(color);
    line->Draw();
    if(span){    
      Float_t pttext=0.06;
      TPaveText *pt = new TPaveText(z-aspan/2,9.125,z-aspan/2,9.125,"br");
      if(roman){
	TText *text = pt->AddText("I");
	//	pt->SetTextFont(132);//22times bold, 132times
	pt->SetTextSize(0.05); 
      }
      else{
	TText *text = pt->AddText("a");
	pt->SetTextSize(pttext); 
      }
      pt->SetFillColor(0); 
      if(bModern)pt->SetShadowColor(0);
      pt->SetLineColor(0);  

      pt->SetTextAlign(21);//Bottom, Middle
      pt->SetFillStyle(0); 
      pt->Draw();

      line = new TLine(z,9,z-aspan,9);
      line->SetLineWidth(3);
      line->Draw();     
    }

    z-=aspan;
    line = new TLine(z,minE,z,maxE);
    line->SetLineStyle(2);//dashed
    line->SetLineWidth(2);
    line->Draw();
  }

  if(edge2!=0){
    z=edge2;
    line = new TLine(z,minE,z,maxE);
    line->SetLineStyle(1);
    line->SetLineWidth(2);
    line->SetLineColor(color);
    line->Draw();
    if(span){
      pt = new TPaveText(z-aspan/2,9-gap+0.125,z-aspan/2,9-gap+0.125,"br");
      if(roman){
	text = pt->AddText("II");
	//pt->SetTextFont(132);
	pt->SetTextSize(0.05); 
      }
      else{
	TText *text = pt->AddText("b");      
	pt->SetTextSize(pttext); 
      }
      pt->SetFillColor(0); 
      if(bModern)pt->SetShadowColor(0);
      pt->SetLineColor(0);  

      pt->SetTextAlign(21);//Bottom, Middle
      pt->SetFillStyle(0); 
      pt->Draw();

      line = new TLine(z,9-gap,z-aspan,9-gap);
      line->SetLineWidth(3);
      line->Draw();
    }
    z-=aspan;
    line = new TLine(z,minE,z,maxE);
    line->SetLineStyle(1);//solid
    line->SetLineWidth(2);
    line->Draw();
  }

  if(edge3!=0){
    z=edge3;
    line = new TLine(z,minE,z,maxE);
    line->SetLineStyle(4);//dot-dash
    line->SetLineWidth(2);
    line->SetLineColor(color);
    line->Draw();
    if(span){   
      pt = new TPaveText(z-aspan/2,9-2*gap+0.125,z-aspan/2,9-2*gap+0.125,"br");
      if(roman){
	text = pt->AddText("III");
	//	pt->SetTextFont(132);
	pt->SetTextSize(0.05); 
      }
      else{
	text = pt->AddText("c");
	pt->SetTextSize(pttext); 
      }
      pt->SetFillColor(0); 
      if(bModern)pt->SetShadowColor(0);
      pt->SetLineColor(0);  
      pt->SetTextAlign(21);//Bottom, Middle
      pt->SetFillStyle(0); 
      pt->Draw();

      line = new TLine(z,9-2*gap,z-aspan,9-2*gap);
      line->SetLineWidth(3);
      line->Draw();
    }
    z-=aspan;
    line = new TLine(z,minE,z,maxE);
    line->SetLineStyle(4);
    line->SetLineWidth(2);
    line->Draw();
  }
}

void plothlines(Float_t edge1=0., Float_t edge2=0., Float_t edge3=0.,
		Float_t xmin=-1e5,Float_t xmax=1e5,Int_t color=1)
{
  Float_t z;
  if(edge1!=0){
    z=edge1;
    TLine *line = new TLine(xmin,z,xmax,z);
    line->SetLineStyle(2);//dashed
    line->SetLineWidth(2);
    line->SetLineColor(color);
    line->Draw();
  }

  if(edge2!=0){
    z=edge2;
    line = new TLine(xmin,z,xmax,z);
    line->SetLineStyle(1);//solid
    line->SetLineWidth(2);
    line->SetLineColor(color);
    line->Draw();
  }

  if(edge3!=0){
    z=edge3;
    line = new TLine(xmin,z,xmax,z);
    line->SetLineStyle(4);//dot-dash
    line->SetLineWidth(2);
    line->SetLineColor(color);
    line->Draw();
  }
}

void plotline(Float_t xmin=0,Float_t ymin=0,Float_t xmax=10,Float_t ymax=10,Int_t color=1)
{
  TLine *line = new TLine(xmin,ymin,xmax,ymax);
  line->SetLineStyle(2);//dashed
  line->SetLineWidth(2);
  line->SetLineColor(color);
  line->Draw();
}

void loadPSDsource()
{//works with psd_analyze.cc
  _filename0=new TFile("PSD/run10.root");
  _filename1=new TFile("PSD/run12.root");
  //The TADC class may not be defined before loading .root files!?
  if(gROOT->GetVersionInt()<52400){
    printf("Deprecated (online) ROOT version.\n");
    gROOT->LoadMacro("/net/helios/Si28/offline/lighthall/macros/psd_analyze.cc");  
  }
  else
    gROOT->LoadMacro("C:/root/macros/psd_analyze.cc");  
  
  _filename0->cd();
  analyze();
  hEX->Clone("hEX_10");
  hEX_10->SetDirectory(home);
  hXFXN->Clone("hXFXN_10");
  hXFXN_10->SetDirectory(home);
  hESum->Clone("hESum_10");
  hESum_10->SetDirectory(home); 
  _filename1->cd();
  analyze();
  hEX->Clone("hEX_12");
  hEX_12->SetDirectory(home);
  hXFXN->Clone("hXFXN_12");
  hXFXN_12->SetDirectory(home);
  hESum->Clone("hESum_12");
  hESum_12->SetDirectory(home); 
  home->cd();
  add2("hEX_10","hEX_12","hEX_src");
  //add2("hXFXN_10","hXFXN_12","hXFXN_src");
  //add2("hESum_10","hESum_12","hESum_src");
  dr("hEX_src");
}

void loadPSD()
{//works with psd_analyze.cc - To be run from parent directory (above /PSD)
  _filename0=new TFile("PSD/run10.root"); //masked source             1/23/07
  _filename1=new TFile("PSD/run11.root"); //target at 45 deg (masked) 1/24/07
  _filename2=new TFile("PSD/run12.root"); //masked source (no beam)
  _filename3=new TFile("PSD/run13.root"); //5MeV
  _filename4=new TFile("PSD/run14.root"); //5MeV                      1/25/07
  _filename5=new TFile("PSD/run15.root"); //4MeV
  _filename6=new TFile("PSD/run16.root"); //3MeV
  _filename7=new TFile("PSD/run17.root"); //2MeV
  //The TADC class may not be defined before loading .root files!?
  if(gROOT->GetVersionInt()<52400){//load online/offline analysis file
    printf("Deprecated (online) ROOT version.\n");
    gROOT->LoadMacro("/net/helios/Si28/offline/lighthall/macros/psd_analyze.cc");  
    bModern=0;
  }
  else
    gROOT->LoadMacro("C:/root/macros/psd_analyze.cc");  
  _filename0->cd();
  analyze();
  hEX->Clone("hEX_10");
  hEX_10->SetDirectory(home);
  hXFXN->Clone("hXFXN_10");
  hXFXN_10->SetDirectory(home);
  hESum->Clone("hESum_10");
  hESum_10->SetDirectory(home); 
  _filename1->cd();
  //makehists();
  analyze();
  hEX->Clone("hEX_11");
  hEX_11->SetDirectory(home);
  hXFXN->Clone("hXFXN_11");
  hXFXN_11->SetDirectory(home);
  hESum->Clone("hESum_11");
  hESum_11->SetDirectory(home); 
  _filename2->cd();
  analyze();
  hEX->Clone("hEX_12");
  hEX_12->SetDirectory(home);
  hXFXN->Clone("hXFXN_12");
  hXFXN_12->SetDirectory(home);
  hESum->Clone("hESum_12");
  hESum_12->SetDirectory(home);  
  _filename3->cd();
  analyze();
  hEX->Clone("hEX_13");
  hEX_13->SetDirectory(home);
  hXFXN->Clone("hXFXN_13");
  hXFXN_13->SetDirectory(home);
  hESum->Clone("hESum_13");
  hESum_13->SetDirectory(home); 
  _filename4->cd();
  analyze();
  hEX->Clone("hEX_14");
  hEX_14->SetDirectory(home);
  hXFXN->Clone("hXFXN_14");
  hXFXN_14->SetDirectory(home);
  hESum->Clone("hESum_14");
  hESum_14->SetDirectory(home); 
  _filename5->cd();
  analyze();
  hEX->Clone("hEX_15");
  hEX_15->SetDirectory(home);
  hXFXN->Clone("hXFXN_15");
  hXFXN_15->SetDirectory(home);
  hESum->Clone("hESum_15");
  hESum_15->SetDirectory(home); 
  _filename6->cd();
  analyze();
  hEX->Clone("hEX_16");
  hEX_16->SetDirectory(home);
  hXFXN->Clone("hXFXN_16");
  hXFXN_16->SetDirectory(home);
  hESum->Clone("hESum_16");
  hESum_16->SetDirectory(home); 
  _filename7->cd();
  analyze();
  hEX->Clone("hEX_17");
  hEX_17->SetDirectory(home);
  hXFXN->Clone("hXFXN_17");
  hXFXN_17->SetDirectory(home);
  hESum->Clone("hESum_17");
  hESum_17->SetDirectory(home); 
  home->cd();
  add2("hEX_10","hEX_12","hEX_src");
  add2("hEX_13","hEX_14","hEX_BE5");
  add2("hEX_15","hEX_BE5","hEX_all");
  add2("hEX_16","hEX_all","hEX_all");
  add2("hEX_17","hEX_all","hEX_all");

  add2("hXFXN_10","hXFXN_12","hXFXN_src"); //"source"
  add2("hXFXN_13","hXFXN_14","hXFXN_BE5"); //"beam energy 5MeV"
  add2("hXFXN_src","hXFXN_BE5","hXFXN_all");
  add2("hXFXN_15","hXFXN_all","hXFXN_all"); //4MeV
  add2("hXFXN_16","hXFXN_all","hXFXN_all"); //3MeV
  add2("hXFXN_17","hXFXN_all","hXFXN_all"); //2MeV

  add2("hESum_10","hESum_12","hESum_src");
  add2("hESum_13","hESum_14","hESum_all");
  add2("hESum_15","hESum_all","hESum_all");
  add2("hESum_16","hESum_all","hESum_all");
  add2("hESum_17","hESum_all","hESum_all");
}

void PSDpaperplots()
{
  if(!((TH2F *) gROOT->FindObject("hEX_all")))
    loadPSD();
  TString ytitle2="Energy (MeV)";
  TString xtitle2="Position (mm)"; 
  copy2("hEX_all",2);
  TH2F * hOutput1=(TH2F *) gROOT->FindObject("hEX_all_copy");
  pjy(hOutput1->GetName(),27.5,32.5);
  TH1F * hOutput2=(TH1F *) gROOT->FindObject("yproj");
  
  pjx(hOutput1->GetName(),4.3,5);
  if(gROOT->FindObject("hHigh"))hHigh->Delete();
  xproj->Clone("hHigh");
  TH1F * hOutput3=(TH1F *) gROOT->FindObject("hHigh");
  pjx(hOutput1->GetName(),1.7,2.2);
  if(gROOT->FindObject("hLow"))hLow->Delete();
  xproj->Clone("hLow");
  TH1F * hOutput4=(TH1F *) gROOT->FindObject("hLow");

  hOutput1->SetStats(kFALSE); //turn off statistics box
  hOutput2->SetStats(kFALSE);
  hOutput3->SetStats(kFALSE);
  hOutput4->SetStats(kFALSE);
  hOutput1->SetTitle();//turn off title box
  hOutput2->SetTitle();
  hOutput3->SetTitle();
  hOutput4->SetTitle();
 
  hOutput1->SetAxisRange(6,56.3,"X");
  hOutput2->SetAxisRange(6,54,"X"); 
  hOutput1->SetYTitle(ytitle2);
  hOutput2->SetYTitle("Counts");
  hOutput3->SetYTitle("Counts");
  hOutput4->SetYTitle("Counts");
  hOutput1->GetYaxis()->CenterTitle(1);
  hOutput2->GetYaxis()->CenterTitle(1);
  hOutput3->GetYaxis()->CenterTitle(1);
  hOutput4->GetYaxis()->CenterTitle(1);
  hOutput1->GetYaxis()->SetTitleOffset(0.8); 
  hOutput2->GetYaxis()->SetTitleOffset(1.125);
  hOutput3->GetYaxis()->SetTitleOffset(1.125);
  hOutput4->GetYaxis()->SetTitleOffset(1.125);
  hOutput1->SetXTitle(xtitle2);
  hOutput2->SetXTitle(ytitle2);
  hOutput3->SetXTitle(xtitle2);
  hOutput4->SetXTitle(xtitle2);
  hOutput1->GetXaxis()->CenterTitle(1);
  hOutput2->GetXaxis()->CenterTitle(1);
  hOutput3->GetXaxis()->CenterTitle(1);
  hOutput4->GetXaxis()->CenterTitle(1);
  yproj->Clone("hLog");
  TH1F * hOutput5=(TH1F *) gROOT->FindObject("hLog");  
 
  mkCanvas2("cOut1","cOut1");
  cOut1->SetTopMargin(.02);
  cOut1->SetRightMargin(.02);
  cOut1->SetCanvasSize(width,(UInt_t)(width*ratio1));
  if(bw){
    gStyle->SetPalette(pal);
    cOut1->SetLogz();  
    hOutput1->SetMaximum(-1111);
  }

  hOutput1->Draw("col");
  if(bw)
    linebox(9,0.5,48,1.72,13,9,2);
  else
    linebox(9,0.5,48,1.72,2,9,2);

  TArrow *arrow = new TArrow(7,5.55,9,5.55,0.02,"|>");

  arrow->SetFillStyle(1001);
  if(bw){
    arrow->SetLineColor(13);
    arrow->SetFillColor(13);
  }
  else{
    arrow->SetLineColor(2);
    arrow->SetFillColor(2);
  }
  arrow->SetLineWidth(2);
  arrow->SetAngle(48);
  arrow->Draw(); 

  pt = new TPaveText(8.5,5.8,8.5,5.8,"br");
  text = pt->AddText("#alpha");
  pt->SetTextAlign(22);//Middle, Middle
  pt->SetFillStyle(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(0.045); 
  pt->SetBorderSize(0);
  pt->Draw();

  pt = new TPaveText(40,1.2,43,1.5,"br");
  text = pt->AddText("^{1}H(p,p)");
  pt->SetTextAlign(21);//Bottom, Middle
  pt->SetFillStyle(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(0.045); 
  pt->SetBorderSize(0);
  pt->Draw();

  plotthresh(49.75,6.4,0.112,0,1);

  if(bw)
    cOut1->SaveAs("out1_bw.eps");
  else
    cOut1->SaveAs("out1.eps");

 
  mkCanvas2("cOut2","cOut2");
  cOut2->SetTopMargin(.02);
  cOut2->SetRightMargin(.02);
  cOut2->SetCanvasSize(width,(UInt_t)(width*ratio2));
  hOutput2->Draw(     );
  Float_t pttext=0.06;
  TPaveText *pt = new TPaveText(0.8068182,472.3586,1.636364,754.25,"br");
  TText *text = pt->AddText("^{1}H(p,p)");
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();
  pt = new TPaveText(1.488636,1782.773,2.363636,2057.045,"br");
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  text = pt->AddText("^{12}C(p,p)");
  pt->Draw();
  pt = new TPaveText(5.170455,182.8485,5.931818,479.9773,"br");
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  text = pt->AddText("#alpha");
  pt->Draw();
  cOut2->SaveAs("out2.eps");
  //High-energy position resolution

  mkCanvas2("cOut3","cOut3");
  cOut3->SetTopMargin(.02);
  cOut3->SetRightMargin(.02);
  cOut3->SetCanvasSize(width,(UInt_t)(width*ratio2));
  hOutput3->Draw(     );
  pt = new TPaveText(.877,.772,.965,.934,"NDC");
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  text = pt->AddText("a");
  pt->Draw();
  cOut3->SaveAs("out3.eps");
  mkCanvas2("cOut4","cOut4");
  cOut4->SetTopMargin(.02);
  cOut4->SetRightMargin(.02);
  cOut4->SetCanvasSize(width,(UInt_t)(width*ratio2));
  hOutput4->Draw(     );
  pt = new TPaveText(.877,.772,.965,.934,"NDC");
  text = pt->AddText("b");
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();
  cOut4->SaveAs("out4.eps");
  //Log-scale histogram----------------------
  mkCanvas2("cOut5","cOut5");
  cOut5->SetTopMargin(.02);
  cOut5->SetRightMargin(.02);
  cOut5->SetCanvasSize(width,(UInt_t)(width*ratio2));
  cOut5->SetLogy(1);
  hOutput5->SetMinimum(9);
  hOutput5->Draw(     );
  pt = new TPaveText(0.1933333,0.6735905,0.2816667,0.768546,"brNDC");
  text = pt->AddText("^{1}H(p,p)");
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();
  pt = new TPaveText(2.897727,1177.476,3.772727,3508.627,"br");
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  text = pt->AddText("^{12}C(p,p)");
  pt->Draw();
  pt = new TPaveText(5.181818,129.5942,5.943182,278.0051,"br");
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  text = pt->AddText("#alpha");
  pt->Draw();
  pt = new TPaveText(4.75,32.86088,5.477273,107.2456,"br");
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  text = pt->AddText("^{16}O(p,p)");
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();
  TArrow *arrow = new TArrow(4.772727,25.01111,5.181818,38.24178,0.02,"<|");
  if(bw){
    arrow->SetFillColor(13);
    arrow->SetLineColor(13);  
  }
  else{
    arrow->SetFillColor(2);
    arrow->SetLineColor(2);  
  }
  arrow->SetFillStyle(1001);
  arrow->SetLineWidth(2);
  arrow->SetAngle(48);
  arrow->Draw();
  
  if(bw)
    cOut5->SaveAs("out5_bw.eps");
  else
    cOut5->SaveAs("out5.eps");

  cFit->Close();
}

void loadbigPSD()
{//works with psd_analyze.cc
  _filename0=new TFile("PSD/run04.root");//50V bias, unmasked       10/04/06
  _filename1=new TFile("PSD/run05.root");//60V bias, unmasked
  _filename2=new TFile("PSD/run06.root");//70V bias (no difference)
  _filename3=new TFile("PSD/run07.root");//60V bias, masked, 18hrs
  _filename4=new TFile("PSD/run08.root");//masked, short
  _filename5=new TFile("PSD/run18.root");//50V 1MeV + 241Am         01/26/07
  _filename6=new TFile("PSD/run19.root");//60V
  _filename7=new TFile("PSD/run20.root");//70V
  _filename8=new TFile("PSD/run21.root");//80V
  _filename9=new TFile("PSD/run22.root");//70V 5MeV
  _filename10=new TFile("PSD/run23.root");//70V 4MeV

  //The TADC class may not be defined before loading .root files!?
  if(gROOT->GetVersionInt()<52400){
    printf("Deprecated (online) ROOT version.\n");
    gROOT->LoadMacro("/net/helios/Si28/offline/lighthall/macros/psd_analyze.cc");  
  }
  else
    gROOT->LoadMacro("C:/root/macros/psd_analyze.cc");  
  
  _filename0->cd();
  analyze("Event","temp5.lst",kTRUE);
  hEX->Clone("hEX_04");
  hEX_04->SetDirectory(home);
  hXFXN->Clone("hXFXN_04");
  hXFXN_04->SetDirectory(home);
  hESum->Clone("hESum_04");
  hESum_04->SetDirectory(home); 
  _filename1->cd();
  analyze("Event","temp5.lst",kTRUE);
  hEX->Clone("hEX_05");
  hEX_05->SetDirectory(home);
  hXFXN->Clone("hXFXN_05");
  hXFXN_05->SetDirectory(home);
  hESum->Clone("hESum_05");
  hESum_05->SetDirectory(home); 
  _filename2->cd();
  analyze("Event","temp5.lst",kTRUE);
  hEX->Clone("hEX_06");
  hEX_06->SetDirectory(home);
  hXFXN->Clone("hXFXN_06");
  hXFXN_06->SetDirectory(home);
  hESum->Clone("hESum_06");
  hESum_06->SetDirectory(home); 
  _filename3->cd();
  analyze("Event","temp5.lst",kTRUE);
  hEX->Clone("hEX_07");
  hEX_07->SetDirectory(home);
  hXFXN->Clone("hXFXN_07");
  hXFXN_07->SetDirectory(home);
  hESum->Clone("hESum_07");
  hESum_07->SetDirectory(home); 
  _filename4->cd();
  analyze("Event","temp5.lst",kTRUE);
  hEX->Clone("hEX_08");
  hEX_08->SetDirectory(home);
  hXFXN->Clone("hXFXN_08");
  hXFXN_08->SetDirectory(home);
  hESum->Clone("hESum_08");
  hESum_08->SetDirectory(home); 
  _filename5->cd();
  analyze("Event","temp5.lst",kTRUE,2);
  hEX->Clone("hEX_18");
  hEX_18->SetDirectory(home);
  hXFXN->Clone("hXFXN_18");
  hXFXN_18->SetDirectory(home);
  hESum->Clone("hESum_18");
  hESum_18->SetDirectory(home); 
  hETheta->Clone("hETheta_18");
  hETheta_18->SetDirectory(home); 
  _filename6->cd();
  analyze("Event","temp5.lst",kTRUE,2);
  hEX->Clone("hEX_19");
  hEX_19->SetDirectory(home);
  hXFXN->Clone("hXFXN_19");
  hXFXN_19->SetDirectory(home);
  hESum->Clone("hESum_19");
  hESum_19->SetDirectory(home); 
  hETheta->Clone("hETheta_19");
  hETheta_19->SetDirectory(home); 
  _filename7->cd();
  analyze("Event","temp5.lst",kTRUE,2);
  hEX->Clone("hEX_20");
  hEX_20->SetDirectory(home);
  hXFXN->Clone("hXFXN_20");
  hXFXN_20->SetDirectory(home);
  hESum->Clone("hESum_20");
  hESum_20->SetDirectory(home); 
  hETheta->Clone("hETheta_20");
  hETheta_20->SetDirectory(home); 
  _filename8->cd();
  analyze("Event","temp5.lst",kTRUE);
  hEX->Clone("hEX_21");
  hEX_21->SetDirectory(home);
  hXFXN->Clone("hXFXN_21");
  hXFXN_21->SetDirectory(home);
  hESum->Clone("hESum_21");
  hESum_21->SetDirectory(home); 
  _filename9->cd();
  analyze("Event","temp22.lst",kTRUE,1);
  hEX->Clone("hEX_22");
  hEX_22->SetDirectory(home);
  hXFXN->Clone("hXFXN_22");
  hXFXN_22->SetDirectory(home);
  hESum->Clone("hESum_22");
  hESum_22->SetDirectory(home); 
  hETheta->Clone("hETheta_22");
  hETheta_22->SetDirectory(home); 
  _filename10->cd();
  analyze("Event","temp23.lst",kTRUE,1);
  hEX->Clone("hEX_23");
  hEX_23->SetDirectory(home);
  hXFXN->Clone("hXFXN_23");
  hXFXN_23->SetDirectory(home);
  hESum->Clone("hESum_23");
  hESum_23->SetDirectory(home); 
  hETheta->Clone("hETheta_23");
  hETheta_23->SetDirectory(home); 
  home->cd();
  add2("hEX_05","hEX_06","hEX_unmasked");
  add2("hEX_07","hEX_08","hEX_masked");
  add2("hXFXN_05","hXFXN_06","hXFXN_unmasked");
  add2("hXFXN_07","hXFXN_08","hXFXN_masked");
  add2("hESum_05","hESum_06","hESum_unmasked");
  add2("hESum_07","hESum_08","hESum_masked");
  add2("hEX_18","hEX_19","hEX_BE1");
  add2("hEX_20","hEX_BE1","hEX_BE1");
  add2("hEX_22","hEX_23","hEX_BE45");
  add2("hETheta_18","hETheta_19","hETheta_BE1");
  add2("hETheta_20","hETheta_BE1","hETheta_masked2");
  add2("hETheta_22","hETheta_23","hETheta_BE45");
  add2("hETheta_BE1","hETheta_BE45","hETheta_all");
  add2("hEX_BE1","hEX_BE45","hEX_all");
}

void loadbigPSDsource()
{//works with psd_analyze.cc
  _filename0=new TFile("PSD/run04.root");//50V bias, unmasked       10/04/06
  _filename1=new TFile("PSD/run05.root");//60V bias, unmasked
  _filename2=new TFile("PSD/run06.root");//70V bias (no difference)
  _filename3=new TFile("PSD/run07.root");//60V bias, masked, 18hrs
  _filename4=new TFile("PSD/run08.root");//masked, short
  _filename5=new TFile("PSD/run18.root");//50V 1MeV + 241Am         01/26/07
  _filename6=new TFile("PSD/run19.root");//60V
  _filename7=new TFile("PSD/run20.root");//70V
  _filename8=new TFile("PSD/run21.root");//80V
  _filename9=new TFile("PSD/run22.root");//70V 5MeV
  _filename10=new TFile("PSD/run23.root");//70V 5MeV

  //The TADC class may not be defined before loading .root files!?
  if(gROOT->GetVersionInt()<52400){
    printf("Deprecated (online) ROOT version.\n");
    gROOT->LoadMacro("/net/helios/Si28/offline/lighthall/macros/psd_analyze.cc");  
  }
  else
    gROOT->LoadMacro("C:/root/macros/psd_analyze.cc");  
  
  _filename1->cd();
  analyze("Event","temp5.lst",kTRUE);
  hEX->Clone("hEX_05");
  hEX_05->SetDirectory(home);
  hXFXN->Clone("hXFXN_05");
  hXFXN_05->SetDirectory(home);
  hESum->Clone("hESum_05");
  hESum_05->SetDirectory(home); 
  _filename2->cd();
  analyze("Event","temp5.lst",kTRUE);
  hEX->Clone("hEX_06");
  hEX_06->SetDirectory(home);
  hXFXN->Clone("hXFXN_06");
  hXFXN_06->SetDirectory(home);
  hESum->Clone("hESum_06");
  hESum_06->SetDirectory(home); 

  home->cd();
  add2("hEX_05","hEX_06","hEX_unmasked");
  add2("hXFXN_05","hXFXN_06","hXFXN_unmasked");
  add2("hESum_05","hESum_06","hESum_unmasked");
} 

void bigPSDpaperplots()
{
  if(!((TH2F *) gROOT->FindObject("hEX_masked")))
    loadbigPSD();
  TString ytitle2="Energy (MeV)";
  TString xtitle2="Lab Angle (deg)"; 
  copy2("hETheta_all",1);
  TH2F * hOutput1=(TH2F *) gROOT->FindObject("hETheta_all_copy");
  
  hOutput1->SetStats(kFALSE); //turn off statistics box
  hOutput1->SetTitle();//turn off title box
  hOutput1->SetAxisRange(38,72,"X");
  hOutput1->SetYTitle(ytitle2);
  hOutput1->GetYaxis()->CenterTitle(1);
  hOutput1->GetYaxis()->SetTitleOffset(0.8); 
  hOutput1->SetXTitle(xtitle2);
  hOutput1->GetXaxis()->CenterTitle(1);
   
  mkCanvas2("cOut1","cOut1");
  cOut1->SetTopMargin(.02);
  cOut1->SetRightMargin(.02);
  cOut1->SetCanvasSize(width,(UInt_t)(width*ratio1));
  if(bw){
    gStyle->SetPalette(pal);
    cOut1->SetLogz();  
    hOutput1->SetMaximum(-1111);
  }

  hOutput1->Draw("col");
  if(bw)
    linebox(39.5,1.05,65,2.8,13,9,2);
  else
    linebox(39.5,1.05,65,2.8, 2,9,2);

  TArrow *arrow = new TArrow(38.3,5.55,40,5.55,0.02,"|>");

  arrow->SetFillStyle(1001);
  if(bw){
    arrow->SetLineColor(13);
    arrow->SetFillColor(13);
  }
  else{
    arrow->SetLineColor(2);
    arrow->SetFillColor(2);
  }
  arrow->SetLineWidth(2);
  arrow->SetAngle(48);
  arrow->Draw(); 

  pt = new TPaveText(39.2,5.7,39.2,5.7,"br");
  text = pt->AddText("#alpha");
  pt->SetTextAlign(21);//Bottom, Middle
  pt->SetFillStyle(0); 
  //if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(0.045); 
  pt->SetBorderSize(0);
  pt->Draw();

  pt = new TPaveText(61,2.4,61,2.7,"br");
  text = pt->AddText("^{1}H(p,p)");
  pt->SetTextAlign(21);//Bottom, Middle
  pt->SetFillStyle(0); 
  //if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(0.045); 
  pt->SetBorderSize(0);
  pt->Draw();

  if(bw)
    cOut1->SaveAs("bigout1_bw.eps");
  else
    cOut1->SaveAs("bigout1.eps");
 
  cFit->Close();
}

void loadtimePSD()
{//works with psd_analyze.cc
  _filename1=new TFile("PSD/time/run01.root");//garbage
  _filename2=new TFile("PSD/time/run02.root");//garbage
  _filename3=new TFile("PSD/time/run03.root");//5MeV, masked, 3/26/07
  _filename4=new TFile("PSD/time/run04.root");//5MeV, no mask, 3/27/07
  _filename5=new TFile("PSD/time/run05.root");//5MeV
  _filename6=new TFile("PSD/time/run06.root");//7.5 MeV, coinc
  _filename7=new TFile("PSD/time/run07.root");//7.5 MeV, no coinc
  _filename8=new TFile("PSD/time/run08.root");//9.0 MeV, coinc
  _filename9=new TFile("PSD/time/run09.root");//9.0 MeV, no coinc
  _filename10=new TFile("PSD/time/run10.root");//10 MeV, no coinc
  _filename11=new TFile("PSD/time/run11.root");//10MeV, coinc

  //The TADC class may not be defined before loading .root files!?
  if(gROOT->GetVersionInt()<52400){
    printf("Deprecated (online) ROOT version.\n");
    gROOT->LoadMacro("/net/helios/Si28/offline/lighthall/macros/jcl_psdtime.cc");  
  }
  else
    gROOT->LoadMacro("C:/root/macros/jcl_psdtime.cc");  

  _filename3->cd();
  analyze();
  hEXc->Clone("hEXc_03");
  hEXc_03->SetDirectory(home);
  hESiT->Clone("hESiT_03");
  hESiT_03->SetDirectory(home);
  hTX->Clone("hTX_03");
  hTX_03->SetDirectory(home);
  _filename4->cd();
  analyze();
  hEXc->Clone("hEXc_04");
  hEXc_04->SetDirectory(home);
  hESiT->Clone("hESiT_04");
  hESiT_04->SetDirectory(home);
  hTX->Clone("hTX_04");
  hTX_04->SetDirectory(home);
  _filename5->cd();
  analyze();
  hEXc->Clone("hEXc_05");
  hEXc_05->SetDirectory(home);
  hESiT->Clone("hESiT_05");
  hESiT_05->SetDirectory(home);
  hTX->Clone("hTX_05");
  hTX_05->SetDirectory(home);
  _filename6->cd();
  analyze();
  hEXc->Clone("hEXc_06");
  hEXc_06->SetDirectory(home);
  hESiT->Clone("hESiT_06");
  hESiT_06->SetDirectory(home);
  hTX->Clone("hTX_06");
  hTX_06->SetDirectory(home);
  _filename7->cd();
  analyze();
  hEXc->Clone("hEXc_07");
  hEXc_07->SetDirectory(home);
  hESiT->Clone("hESiT_07");
  hESiT_07->SetDirectory(home);
  hTX->Clone("hTX_07");
  hTX_07->SetDirectory(home);
  _filename8->cd();
  analyze();
  hEXc->Clone("hEXc_08");
  hEXc_08->SetDirectory(home);
  hESiT->Clone("hESiT_08");
  hESiT_08->SetDirectory(home);
  hTX->Clone("hTX_08");
  hTX_08->SetDirectory(home);
  _filename9->cd();
  analyze();
  hEXc->Clone("hEXc_09");
  hEXc_09->SetDirectory(home);
  hESiT->Clone("hESiT_09");
  hESiT_09->SetDirectory(home);
  hTX->Clone("hTX_09");
  hTX_09->SetDirectory(home);
  _filename10->cd();
  analyze();
  hEXc->Clone("hEXc_10");
  hEXc_10->SetDirectory(home);
  hESiT->Clone("hESiT_10");
  hESiT_10->SetDirectory(home);
  hTX->Clone("hTX_10");
  hTX_10->SetDirectory(home);
  _filename11->cd();
  analyze();
  hEXc->Clone("hEXc_11");
  hEXc_11->SetDirectory(home);
  hESiT->Clone("hESiT_11");
  hESiT_11->SetDirectory(home);
  hTX->Clone("hTX_11");
  hTX_11->SetDirectory(home);
  home->cd();
  add2("hEXc_04","hEXc_05","hEXc_BE5");
  add2("hEXc_06","hEXc_07","hEXc_BE75");
  add2("hEXc_08","hEXc_09","hEXc_BE9");
  add2("hEXc_10","hEXc_11","hEXc_BE10");
  add2("hEXc_BE5","hEXc_BE75","hEXc_all");
  add2("hEXc_BE75","hEXc_all","hEXc_all");
  add2("hEXc_BE9","hEXc_all","hEXc_all");
  add2("hEXc_BE10","hEXc_all","hEXc_all");
  
  shiftadd2("hESiT_04",327.68,"hESiT_05");
  hESiT_04_shift_copy->Clone("hESiT_BE5");
  add2("hESiT_06","hESiT_07","hESiT_BE75");
  add2("hESiT_08","hESiT_09","hESiT_BE9");
  add2("hESiT_10","hESiT_11","hESiT_BE10");
  add2("hESiT_BE5","hESiT_BE75","hESiT_all");
  add2("hESiT_BE75","hESiT_all","hESiT_all");
  add2("hESiT_BE9","hESiT_all","hESiT_all");
  add2("hESiT_BE10","hESiT_all","hESiT_all");
} 

void timePSDpaperplots()
{
  if(!((TH2F *) gROOT->FindObject("hESiT_all")))
    loadtimePSD();
  TString ytitle2="Energy (MeV)";
  TString xtitle2="Time (ns)"; 
  slopexy("hESiT_BE5",20.48,0,347.891,-48.328);
  slopexy("hESiT_BE75",20.48,0,347.891,-48.328);
  slopexy("hESiT_BE9",20.48,0,347.891,-48.328);
  slopexy("hESiT_BE10",20.48,0,347.891,-48.328);
  slopexy("hESiT_all",20.48,0,347.891,-48.328);
  TH2F * hOutput1=(TH2F *) gROOT->FindObject("hESiT_all_slope");
  
  hOutput1->SetStats(kFALSE); //turn off statistics box
  hOutput1->SetTitle();//turn off title box
  hOutput1->SetAxisRange(15,140,"X");
  hOutput1->SetAxisRange(0,8.5,"Y");
  hOutput1->SetYTitle(ytitle2);
  hOutput1->GetYaxis()->CenterTitle(1);
  hOutput1->GetYaxis()->SetTitleOffset(0.8); 
  hOutput1->SetXTitle(xtitle2);
  hOutput1->GetXaxis()->CenterTitle(1);
  
  mkCanvas2("cOut1","cOut1");
  cOut1->SetTopMargin(.02);
  cOut1->SetRightMargin(.02);
  cOut1->SetCanvasSize(width,(UInt_t)(width*ratio1));
  if(bw){
    gStyle->SetPalette(pal);
    cOut1->SetLogz();  
    hOutput1->SetMaximum(-1111);
    hOutput1->SetMinimum(20);
  }
  else{
    hOutput1->SetMinimum(3);
    hOutput1->SetMaximum(800);
  }
  

  hOutput1->Draw("col");
 
  TArrow *arrow = new TArrow(64,6.9,64,6,0.02,"|>");

  arrow->SetFillStyle(1001);
  if(bw){
    arrow->SetLineColor(13);
    arrow->SetFillColor(13);
  }
  else{
    arrow->SetLineColor(2);
    arrow->SetFillColor(2);
  }
  arrow->SetLineWidth(2);
  arrow->SetAngle(48);
  arrow->Draw(); 

  if(bw)
    cOut1->SaveAs("timetest_bw.eps");
  else
    cOut1->SaveAs("timetest.eps");
 
  cFit->Close();
}

void lable(Char_t *histin="hEZ", TString ytitle2="Energy (MeV)",  TString xtitle2="Time (ns)")
{
  hname=histin;
  if(!((TH2F *) gROOT->FindObject(hname.Data())))
    printf("Histogram \"%s\" not recognized.\n",hname);
    
  
  TH2F * hOutput1=(TH2F *) gROOT->FindObject(hname.Data());  
  
  hOutput1->SetStats(kFALSE); //turn off statistics box
  hOutput1->SetTitle();//turn off title box
  //  hOutput1->SetAxisRange(15,140,"X");
  //hOutput1->SetAxisRange(0,8.5,"Y");
  hOutput1->SetYTitle(ytitle2);
  hOutput1->GetYaxis()->CenterTitle(1);
  hOutput1->GetYaxis()->SetTitleOffset(0.8); 
  hOutput1->SetXTitle(xtitle2);
  hOutput1->GetXaxis()->CenterTitle(1);
  hOutput1->Draw("col");

}

void algor1()
{//runs 2000,3037,3038
  _filename0=new TFile("900.root");
  copy2("hXFXN10",1);
  
  fitpfx("hXFXN10_copy",0,4000,0,4000,2);
  hname="hXFXN10_copy_pfx";
  if(!((TH2F *) gROOT->FindObject(hname.Data())))
    printf("Histogram \"%s\" not recognized.\n",hname);

  TH2F * hOutput1=(TH2F *) gROOT->FindObject(hname.Data());  
  TString ytitle2="X_{far} (channel no.)";
  TString xtitle2="X_{near} (channel no.)";
  ytitle2="XF (channel no.)";
  xtitle2="XN (channel no.)";

  hOutput1->SetStats(kFALSE); //turn off statistics box
  hOutput1->SetTitle();//turn off title box
  hOutput1->SetAxisRange(0,3500,"X");
  hOutput1->SetAxisRange(0,3500,"Y");
  hOutput1->SetYTitle(ytitle2);
  hOutput1->GetYaxis()->CenterTitle(1);
  hOutput1->GetYaxis()->SetTitleOffset(1.10); 
  hOutput1->SetXTitle(xtitle2);
  hOutput1->GetXaxis()->CenterTitle(1);
  hOutput1->GetXaxis()->SetTitleOffset(1.24); 

  mkCanvas2("cXFXN","cXFXN");
  cXFXN->SetTopMargin(.02);
  cXFXN->SetRightMargin(.02);
  cXFXN->SetCanvasSize(width,(UInt_t)(width*ratio2));

  Float_t  p0=hProf->GetFunction("pol2")->GetParameter(0);//"c"
  Float_t  p1=hProf->GetFunction("pol2")->GetParameter(1);//"b"
  Float_t  p2=hProf->GetFunction("pol2")->GetParameter(2);//"a"
  Float_t zero1=(-p1-sqrt(p1*p1-(4*p2*p0)))/(2*p2);

  TLine *line = new TLine(0,p0,zero1,0);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  hXFXN10_copy_pfx->SetLineColor(1);  
  hXFXN10_copy_pfx->SetMarkerColor(2);
  hXFXN10_copy_pfx->SetMarkerStyle(20);
  hXFXN10_copy_pfx->SetMarkerSize(0.8);

  pol2->SetLineColor(1);
 
  hOutput1->Draw("col");
  line->Draw();

  cXFXN->SaveAs("cXFXN.eps");

  ytitle2="E (channel no.)";
  xtitle2="X (relative position)";
 
  copy2("hEX10",0);
  
  fitpfx("hEX10_copy",2000,2800,0,1,6);

  TH2F * hOutput2=(TH2F *) gROOT->FindObject("hEX10_copy_pfx");
  
  hOutput2->SetStats(kFALSE); //turn off statistics box
  hOutput2->SetTitle();//turn off title box
  hOutput2->SetAxisRange(-.05,1.05,"X");

  hOutput2->SetAxisRange(0,3500,"Y");
  hOutput2->SetYTitle(ytitle2);
  hOutput2->GetYaxis()->CenterTitle(1);
  hOutput2->GetYaxis()->SetTitleOffset(1.10); 
  hOutput2->SetXTitle(xtitle2);
  hOutput2->GetXaxis()->CenterTitle(1);
  hOutput2->GetXaxis()->SetTitleOffset(1.24); 
  
  hEX10_copy_pfx->SetLineColor(1);  
  hEX10_copy_pfx->SetMarkerColor(2);
  hEX10_copy_pfx->SetMarkerStyle(20);
  hEX10_copy_pfx->SetMarkerSize(0.8);

  p0=hProf->GetFunction("pol6")->GetParameter(0);//"c"

  TLine *line = new TLine(0,p0,1,p0);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  mkCanvas2("cEX","cEX");
  cEX->SetTopMargin(.02);
  cEX->SetRightMargin(.02);
  cEX->SetCanvasSize(width,(UInt_t)(width*ratio2));
  hOutput2->Draw("col");
  line->Draw("same");

  cEX->SaveAs("cEX.eps");
  
  cFit->Close();
}

void algor2(Bool_t tall=0)
{//uses data from run 27
  Float_t pttext=0.06;
  TPaveText *pt;
  Float_t aspect=0;
  Float_t yoffset=0;
  Float_t xoffset=0;
  Float_t margin=0;
  if(tall){
    aspect=ratio3;
    yoffset=1.80;
    margin=0.15;
    xoffset=1.00;
  }
  else{
    aspect=ratio1;
    yoffset=1.40;
    margin=0.11;
    xoffset=1.24;
  }

  _filename0=new TFile("500_sum.root");
 
  //1.) hEX gated------------------------
  TH2F * hOutput1=(TH2F *) gROOT->FindObject("hEX19");  

  TString ytitle2="X_{far} (channel no.)";
  TString xtitle2="X_{near} (channel no.)";

  ytitle2="E (channel no.)";
  xtitle2="X (relative position)";

  hOutput1->SetStats(kFALSE); //turn off statistics box
  hOutput1->SetTitle();//turn off title box
  hOutput1->SetAxisRange(-.05,1.05,"X");
  hOutput1->SetAxisRange(0,3500,"Y");

  hOutput1->SetYTitle(ytitle2);
  hOutput1->GetYaxis()->CenterTitle(1);
  hOutput1->GetYaxis()->SetTitleOffset(yoffset); 
  hOutput1->SetXTitle(xtitle2);
  hOutput1->GetXaxis()->CenterTitle(1);
  hOutput1->GetXaxis()->SetTitleOffset(xoffset); 

  mkCanvas2("cEX_sum","cEX_sum");
  cEX_sum->SetTopMargin(.02);
  cEX_sum->SetRightMargin(.02);
  cEX_sum->SetLeftMargin(margin);
  cEX_sum->SetCanvasSize(width,(UInt_t)(width*aspect));

  hOutput1->Draw("col");

  pt = new TPaveText(0.95,3250,0.95,3250,"br");
  text = pt->AddText("d");
  pt->SetTextAlign(12);//Middle, Left
  pt->SetFillColor(0); 
  //    pt->SetFillStyle(0);
  //  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();

  cEX_sum->SaveAs("cEX_sum.eps");

  //2.) hESum projection-----------------
  mkCanvas2("cESum_pjy","cESum_pjy");
  cESum_pjy->cd();
  pjy("hESum19");

  TH1F * hOutput2=(TH1F *) gROOT->FindObject("yproj"); 

  hOutput2->SetStats(kFALSE); //turn off statistics box
  hOutput2->SetTitle();//turn off title box

  ytitle2="Counts";
  xtitle2="E-(XF+XN) (channel no.)";

  hOutput2->SetAxisRange(-100,100,"X");
  //  hOutput2->SetAxisRange(0,500,"Y");
  hOutput2->SetYTitle(ytitle2);
  hOutput2->GetYaxis()->CenterTitle(1);
  hOutput2->GetYaxis()->SetTitleOffset(yoffset); 
  hOutput2->SetXTitle(xtitle2);
  hOutput2->GetXaxis()->CenterTitle(1);
  hOutput2->GetXaxis()->SetTitleOffset(xoffset); 

  cESum_pjy->SetTopMargin(.02);
  cESum_pjy->SetRightMargin(.03);
  cESum_pjy->SetLeftMargin(margin);
  cESum_pjy->SetCanvasSize(width,(UInt_t)(width*aspect));

  hOutput2->Draw("col");
  yproj->Draw();

  pt = new TPaveText(80,450,80,450,"br");
  text = pt->AddText("b");
  pt->SetTextAlign(12);//Middle, Left
  pt->SetFillColor(0); 
  //    pt->SetFillStyle(0);
  //  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();

  cESum_pjy->SaveAs("cESum_pjy.eps");

  //5.) peakfit spectrum-------------------------------------
  peakfity("hEX19","7alphas.lst");
  TH1F * hOutput5=(TH1F *) gROOT->FindObject("hEX19_py"); 
  hOutput5->SetAxisRange(1200,2500,"X");
  peakfity("hEX19","7alphas.lst",2,1.1,.2);

  hOutput5->SetStats(kFALSE); //turn off statistics box
  hOutput5->SetTitle();//turn off title box

  ytitle2="Counts";
  xtitle2="E (channel no.)";

  //  hOutput5->SetAxisRange(-100,100,"X");
  //  hOutput5->SetAxisRange(0,500,"Y");
  hOutput5->SetYTitle(ytitle2);
  hOutput5->GetYaxis()->CenterTitle(1);
  hOutput5->GetYaxis()->SetTitleOffset(1.40); 
  hOutput5->SetXTitle(xtitle2);
  hOutput5->GetXaxis()->CenterTitle(1);
  hOutput5->GetXaxis()->SetTitleOffset(1.24); 

  mkCanvas2("cEpeaks","cEpeaks");
  cEpeaks->SetTopMargin(.02);
  cEpeaks->SetRightMargin(.03);
  cEpeaks->SetLeftMargin(margin);
  cEpeaks->SetCanvasSize(width,(UInt_t)(width*ratio2));

  hOutput5->Draw("col");
  hEX19_py->GetFunction("gaus")->Delete();
  hEX19_py->Draw();
  cFit->Close();

  cEpeaks->SaveAs("cEpeaks.eps");

  //3.) hEX ungated------------------------------------------
  _filename1=new TFile("500_nosum.root");
 
  TH2F * hOutput3=(TH2F *) gROOT->FindObject("hEX19");  

  ytitle2="E (channel no.)";
  xtitle2="X (relative position)";

  hOutput3->SetStats(kFALSE); //turn off statistics box
  hOutput3->SetTitle();//turn off title box
  hOutput3->SetAxisRange(-.05,1.05,"X");
  hOutput3->SetAxisRange(0,3500,"Y");

  hOutput3->SetYTitle(ytitle2);
  hOutput3->GetYaxis()->CenterTitle(1);
  hOutput3->GetYaxis()->SetTitleOffset(yoffset); 
  hOutput3->SetXTitle(xtitle2);
  hOutput3->GetXaxis()->CenterTitle(1);
  hOutput3->GetXaxis()->SetTitleOffset(xoffset); 

  mkCanvas2("cEX_nosum","cEX_nosum");
  cEX_nosum->SetTopMargin(.02);
  cEX_nosum->SetRightMargin(.02);
  cEX_nosum->SetLeftMargin(margin);
  cEX_nosum->SetCanvasSize(width,(UInt_t)(width*aspect));

  hOutput3->Draw("col");

  TArrow *arrow = new TArrow(0.035,2900,0.035,2500,0.02,"|>");
  
  arrow->SetFillStyle(1001);
  if(bw){
    arrow->SetLineColor(13);
    arrow->SetFillColor(13);
  }
  else{
    arrow->SetLineColor(2);
    arrow->SetFillColor(2);
  }
  arrow->SetLineWidth(2);
  arrow->SetAngle(48);
  arrow->Draw(); 

  pt = new TPaveText(0.95,3250,0.95,3250,"br");
  text = pt->AddText("c");
  pt->SetTextAlign(12);//Middle, Left
  pt->SetFillColor(0); 
  //    pt->SetFillStyle(0);
  //  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();

  cEX_nosum->SaveAs("cEX_nosum.eps");

  //4.) hEsum------------------------------------------------
  TH2F * hOutput4=(TH2F *) gROOT->FindObject("hESum19");  
 
  ytitle2="E (channel no.)";
  xtitle2="(XF+XN) (channel no.)";
  
  hOutput4->SetStats(kFALSE); //turn off statistics box
  hOutput4->SetTitle();//turn off title box
  hOutput4->SetAxisRange(0,3500,"X");
  hOutput4->SetAxisRange(0,3500,"Y");

  hOutput4->SetYTitle(ytitle2);
  hOutput4->GetYaxis()->CenterTitle(1);
  hOutput4->GetYaxis()->SetTitleOffset(yoffset); 
  hOutput4->SetXTitle(xtitle2);
  hOutput4->GetXaxis()->CenterTitle(1);
  hOutput4->GetXaxis()->SetTitleOffset(xoffset); 

  mkCanvas2("cESum","cESum");
  cESum->SetTopMargin(.02);
  cESum->SetRightMargin(.04);
  cESum->SetLeftMargin(margin);
  cESum->SetCanvasSize(width,(UInt_t)(width*aspect));

  hOutput4->Draw("col");

  TArrow *arrow = new TArrow(2750,2750,2450,2450,0.02,"|>");
  
  arrow->SetFillStyle(1001);
  if(bw){
    arrow->SetLineColor(13);
    arrow->SetFillColor(13);
  }
  else{
    arrow->SetLineColor(2);
    arrow->SetFillColor(2);
  }
  arrow->SetLineWidth(2);
  arrow->SetAngle(48);
  arrow->Draw(); 

  pt = new TPaveText(3200,3250,3200,3250,"br");
  text = pt->AddText("a");
  pt->SetTextAlign(12);//Middle, Left
  pt->SetFillColor(0); 
  //    pt->SetFillStyle(0);
  //  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();

  cESum->SaveAs("cESum.eps");
}

void algor3()
{//made from runs run0002[3,5,6] (thick target setting)
  _filename0=new TFile("500_time.root");
  fitpqpfy("hET1",400,1100,400,3000);

  TH2F * hOutput1=(TH2F *) gROOT->FindObject("hET1_pfy");  
  TString ytitle2="X_{far} (channel no.)";
  TString xtitle2="X_{near} (channel no.)";

  ytitle2="T (channel no.)";
  xtitle2="E (channel no.)";

  hOutput1->SetStats(kFALSE); //turn off statistics box
  hOutput1->SetTitle();//turn off title box
  hOutput1->SetAxisRange(0,2600,"X");
  hOutput1->SetAxisRange(400,1100,"Y");
  hOutput1->SetYTitle(ytitle2);
  hOutput1->GetYaxis()->CenterTitle(1);
  hOutput1->GetYaxis()->SetTitleOffset(1.10); 
  hOutput1->SetXTitle(xtitle2);
  hOutput1->GetXaxis()->CenterTitle(1);
  hOutput1->GetXaxis()->SetTitleOffset(1.24); 

  mkCanvas2("cET_pfy","cET_pfy");
  cET_pfy->SetTopMargin(.02);
  cET_pfy->SetRightMargin(.02);
  cET_pfy->SetCanvasSize(width,(UInt_t)(width*ratio2));

  Float_t  p0=hProf->GetFunction("pol2")->GetParameter(0);//"c"
  Float_t  p1=hProf->GetFunction("pol2")->GetParameter(1);//"b"
  Float_t  p2=hProf->GetFunction("pol2")->GetParameter(2);//"a"
  Float_t zero1=(-p1-sqrt(p1*p1-(4*p2*p0)))/(2*p2);
  Float_t cp=(-p1/(2*p2));
  
  TLine *line = new TLine(0,p0,zero1,0);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  hET1_pfy->SetLineColor(1);  
  hET1_pfy->SetMarkerColor(2);
  hET1_pfy->SetMarkerStyle(20);
  hET1_pfy->SetMarkerSize(0.8);

  pol2->SetLineColor(1);
 
  hOutput1->Draw("col");

  cET_pfy->SaveAs("cET_pfy.eps");

  ytitle2="E (channel no.)";
  xtitle2="T (channel no.)";
 
  copy2("hET1",0);
  
  TH2F * hOutput2=(TH2F *) gROOT->FindObject("hET1_copy");
  
  hOutput2->SetStats(kFALSE); //turn off statistics box
  hOutput2->SetTitle();//turn off title box
  hOutput2->SetAxisRange(0,1100,"X");

  hOutput2->SetAxisRange(0,2600,"Y");
  hOutput2->SetYTitle(ytitle2);
  hOutput2->GetYaxis()->CenterTitle(1);
  hOutput2->GetYaxis()->SetTitleOffset(1.40); 
  hOutput2->SetXTitle(xtitle2);
  hOutput2->GetXaxis()->CenterTitle(1);
  hOutput2->GetXaxis()->SetTitleOffset(1.24); 
  
  mkCanvas2("cET","cET");
  cET->SetTopMargin(.02);
  cET->SetRightMargin(.02);
  cET->SetLeftMargin(.11);
  cET->SetCanvasSize(width,(UInt_t)(width*ratio1));
  hOutput2->Draw("col");

  cET->SaveAs("cET.eps");
  
  cFit->Close();
}

void algor4(Bool_t tall=1)
{
  gStyle->SetOptDate(0);
  Float_t aspect=0;
  Float_t yoffset=0;
  Float_t xoffset=0;
  Float_t margin=0;
  if(tall){
    aspect=ratio3;
    yoffset=1.80;
    margin=0.15;
    xoffset=1.00;
  }
  else{
    aspect=ratio1;
    yoffset=1.40;
    margin=0.11;
    xoffset=1.24;
  }

  loadfiles();
  _filename3=new TFile("100_algor.root");
  _filename0->cd();
  hTX21->Clone("hTX_cal");
  hTX_cal->SetDirectory(home);
  //  plotall("hTX");
  _filename3->cd();
  hTX21->Clone("hTX_nocal");
  hTX_nocal->SetDirectory(home);
  //plotall("hTX",0,0,0,0,1,"cFit2");
  home->cd();
  
  TH2F * hOutput1=(TH2F *) gROOT->FindObject("hTX_cal");  
 
  TString xtitle2="X (relative position)";
  TString ytitle2="Time (ns)";
  
  hOutput1->SetStats(kFALSE); //turn off statistics box
  hOutput1->SetTitle();//turn off title box
  hOutput1->SetAxisRange(-.05,1.05,"X");
  hOutput1->SetAxisRange(0,82,"Y");

  hOutput1->SetYTitle(ytitle2);
  hOutput1->GetYaxis()->CenterTitle(1);
  hOutput1->GetYaxis()->SetTitleOffset(yoffset); 
  hOutput1->SetXTitle(xtitle2);
  hOutput1->GetXaxis()->CenterTitle(1);
  hOutput1->GetXaxis()->SetTitleOffset(xoffset); 
  hOutput1->SetMinimum(6);

  mkCanvas2("cTcal","cTcal");
  cTcal->SetTopMargin(.02);
  cTcal->SetRightMargin(.04);
  cTcal->SetLeftMargin(margin);
  cTcal->SetCanvasSize(width,(UInt_t)(width*aspect));

  hOutput1->Draw("col");
  cTcal->SaveAs("cTcal.eps");

  TH2F * hOutput2=(TH2F *) gROOT->FindObject("hTX_nocal");  
  ytitle2="Time (channel no.)"; 

  hOutput2->SetStats(kFALSE); //turn off statistics box
  hOutput2->SetTitle();//turn off title box
  hOutput2->SetAxisRange(0,1,"X");
  Float_t offset=0;
  hOutput2->SetAxisRange(0+offset,1476+offset,"Y");

  hOutput2->SetYTitle(ytitle2);
  hOutput2->GetYaxis()->CenterTitle(1);
  hOutput2->GetYaxis()->SetTitleOffset(yoffset); 
  hOutput2->SetXTitle(xtitle2);
  hOutput2->GetXaxis()->CenterTitle(1);
  hOutput2->GetXaxis()->SetTitleOffset(xoffset); 
  hOutput2->SetMinimum(6);
  hOutput2->SetMaximum(145);

  mkCanvas2("cTnocal","cTnocal");
  cTnocal->SetTopMargin(.02);
  cTnocal->SetRightMargin(.04);
  cTnocal->SetLeftMargin(margin);
  cTnocal->SetCanvasSize(width,(UInt_t)(width*aspect));

  hOutput2->Draw("col");
  cTnocal->SaveAs("cTnocal.eps");
}

void deffunc(Float_t b=2048)
{
  fEDiff1=new TF1("fEDiff1","x",0,4096);
  fEDiff2=new TF1("fEDiff2","-x",-4096,0);
  fXFXN=new TF1("fXFXN","[0]-x",0,4096);
  fXFXN->SetParameter(0,b);
  fXFXN->SetLineWidth(2);
  fXFXN->SetLineStyle(9);
  fXFXN->SetLineColor(3);
}

void linebox(Float_t X1, Float_t Y1, Float_t X2, Float_t Y2, Int_t color=1, Int_t style=1, Int_t width=1)
{//draws a box composed of lines
  TLine *line;
  line = new TLine(X1,Y1,X2,Y1);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(width);
  line->Draw();
  line = new TLine(X2,Y1,X2,Y2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(width);
  line->Draw();
  line = new TLine(X2,Y2,X1,Y2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(width);
  line->Draw();
  line = new TLine(X1,Y2,X1,Y1);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(width);
  line->Draw();
}

void plotthresh(Float_t L=1, Float_t start=0, Float_t thresh=.049, Bool_t bVerb=kTRUE,Bool_t bReplace=kFALSE)
{
  TString fname="thresh.dat";
  Float_t X=0,cutoff=0;
  Int_t steps=200;
  FILE * outfile;
  if(bReplace&&gROOT->FindObject("gline"))
    gROOT->FindObject("gline")->Delete();
  
  outfile=fopen(fname.Data(),"w");
  if(bVerb) printf("The contents of \"%s\" are: ",fname.Data());
  for(Int_t i=1;i<(steps-1);i++){
    X=(Float_t)i/steps;
    if(bVerb) printf("i = %2d, X = %.3f, ",i, X);
    if(X<=0.5)
      cutoff=thresh/X;
    else
      cutoff=thresh/(1.0-X);

    X=X*L+start;
    fprintf(outfile,"%g  %g\n",X,cutoff);
    if(bVerb) printf("%.3f  %.3f\n",X,cutoff);
  }
  fclose(outfile);
  drawline(fname.Data(),"gline",0,2,2);
}

void ploteffic(Float_t L=1, Float_t start=0, Float_t top=1, Float_t factor=10, Bool_t bVerb=kTRUE,Bool_t bReplace=kFALSE)
{
  TString fname="thresh.dat";
  Float_t X=0,cutoff=0;
  Int_t steps=200;
  FILE * outfile;
  if(bReplace)
    if(gROOT->FindObject("gline"))
      gROOT->FindObject("gline")->Delete();
  
  outfile=fopen(fname.Data(),"w");
  if(bVerb) printf("The contents of \"%s\" are: ",fname.Data());
  for(Int_t i=1;i<(steps-1);i++){
    X=(Float_t)i/steps;
    if(bVerb) printf("i = %2d, X = %.3f, ",i, X);
    
    cutoff=top-factor*fabs(X-0.5);

    X=X*L+start;
    fprintf(outfile,"%g  %g\n",X,cutoff);
    if(bVerb) printf("%.3f  %.3f\n",X,cutoff);
  }
  fclose(outfile);
  drawline(fname.Data(),"gline",0,2,2);
}

void plotcuts(Int_t set=0, Float_t thresh=0.21, Float_t top=5, Float_t factor=17,Float_t L=50.5 )
{
  for(Int_t i=0;i<6;i++){
    plotthresh(L,edges[set][i],thresh,0,0);
    ploteffic(L,edges[set][i],top,factor,0,0);
  }
}

void clearlines(TString name="gline")
{
  
  while (gROOT->FindObject(name.Data()))
    gROOT->FindObject(name.Data())->Delete();
}

void B12plots()
{
  gStyle->SetOptDate(0);
  bw=1;
  _filename2 = new TFile("helios_b11dp_off3_scl_ahw_hy.root");
  TH2F * hOutput4=(TH2F *) gROOT->FindObject("hEZg4");
  //hOutput4->SetAxisRange(plot_minZ,plot_maxZ,"X");
  maxE=3.5;
  hOutput4->SetAxisRange(0,maxE,"Y");
  hOutput4->SetYTitle(ytitle);
  hOutput4->GetYaxis()->CenterTitle(1);
  hOutput4->SetXTitle(xtitle);
  hOutput4->GetXaxis()->CenterTitle(1);
  hOutput4->SetStats(kFALSE);
  hOutput4->SetTitle();
  //  hOutput4->SetMinimum(8);
  mkCanvas2("cB11","cB11");
  cB11->SetTopMargin(.02);
  cB11->SetRightMargin(.03);  
  cB11->SetCanvasSize(width,(UInt_t)(width*ratio1));
  cB11->SetLogz();  
  //  if(bw){gStyle->SetPalette(pal);hOutput4->SetMaximum(-1111);}
  hOutput4->Rebin2D(hOutput4->GetNbinsX()/512,hOutput4->GetNbinsY()/512);
  hOutput4->Draw("col");
  B0=1.0424;
  plotlabangle(4.5,0,0,1.007,1,B0);
  reaction=1;//correct
  plotbore(.462,1.007,1,B0);
  aspan=345.30;

  plotvlines(-368.7);

  Float_t pttext=0.06;
  TPaveText *pt;
  TText *text;  

  plotlines("_B11line.txt",4,7);//missing!

  pt = new TPaveText(-604,3.23,-590.9,3.41,"br");
  text = pt->AddText("a");
  pt->SetTextAlign(12);//Middle, Left
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();
  
  pt = new TPaveText(-497.8,3.09,-484.5,3.27,"br");
  text = pt->AddText("b");
  pt->SetTextAlign(12);//Middle, Left
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();
  
  pt = new TPaveText(-397,2.72,-384,2.89,"br");
  text = pt->AddText("c");
  pt->SetTextAlign(12);//Middle, Left
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();
  
  pt = new TPaveText(-366,2.3,-360,2.3,"br");
  text = pt->AddText("d");
  pt->SetTextAlign(12);//Middle, Left
  pt->SetFillColor(0); 
  if(bModern)pt->SetShadowColor(0);
  pt->SetLineColor(0);  
  pt->SetTextSize(pttext); 
  pt->Draw();
  
  cB11->SaveAs("cB11.eps");
}

void B13plots()
{
  gStyle->SetOptDate(0);
  plotsim();
  filltree("simulations/b12_d_p_sample.dat","tree","simulations/B13_param.dat");
  fixsim();

  TH2F * hOutput4=(TH2F *) gROOT->FindObject("hEZ_gated_shift");
  hOutput4->SetAxisRange(-720,-360,"X");
  maxE=3.5;
  hOutput4->SetAxisRange(0,maxE,"Y");
  hOutput4->SetYTitle(ytitle);
  hOutput4->GetYaxis()->CenterTitle(1);
  hOutput4->SetXTitle(xtitle);
  hOutput4->GetXaxis()->CenterTitle(1);
  hOutput4->SetStats(kFALSE);
  hOutput4->SetTitle();
  //  hOutput4->SetMinimum(8);
  mkCanvas2("cB11","cB11");
  cB11->SetTopMargin(.02);
  cB11->SetRightMargin(.03);  
  cB11->SetCanvasSize(width,(UInt_t)(width*ratio1));
  cB11->SetLogz();  
  //  if(bw){gStyle->SetPalette(pal);hOutput4->SetMaximum(-1111);}
  hOutput4->Rebin2D(hOutput4->GetNbinsX()/512,hOutput4->GetNbinsY()/512);
  hOutput4->Draw("col");
  B0=1.0424;
  plotlabangle(4.5,0,0,1.007,1,B0);
  reaction=1;//wrong!
  plotbore(.462,1.007,1,B0);
  aspan=345.30;

  plotvlines(-368.7);

  cB11->SaveAs("cB12.eps");
}

void Oxplots()
{
  gStyle->SetOptDate(0);
  plotsim();
  cFit->Close();
  filltree("simulations/Ox20.dat","tree","simulations/Ox20_param.dat");
  TH2F * hOutput4=(TH2F *) gROOT->FindObject("hPhiCoverage2");
 
  hOutput4->SetAxisRange(0,82,"Y");
  hOutput4->SetYTitle("Counts per degree");
  hOutput4->GetYaxis()->CenterTitle(1);
  hOutput4->SetXTitle("Angle of Rotation (deg)");
  hOutput4->GetXaxis()->CenterTitle(1);
  hOutput4->SetStats(kFALSE);
  hOutput4->SetTitle();

  mkCanvas2("cPhiCov","cPhiCov");
  cPhiCov->SetTopMargin(.02);
  cPhiCov->SetRightMargin(.03);  
  cPhiCov->SetCanvasSize(width,(UInt_t)(width*ratio1));
 
  hPhiCoverage2->Draw();  
  hPhiCoverage1->Draw("same");

  cPhiCov->SaveAs("cPhiCov.eps");
}

void simplots2()
{
  gStyle->SetOptDate(0);
  plotsim();
  filltree("simulations/21res84tgt2map1bsp_94_ztgt0.dat","tree1","simulations/094_param.dat");
  pjy("hQ2Z");
  slopex("yproj",1.014,0.088);
  yproj_slope->Clone("hmap2");
  
  hQ2Z->Reset();
  filltree("simulations/21res84tgt0map1bsp_94.dat","tree2","simulations/094_param.dat");
  pjy("hQ2Z");
  slopex("yproj",1.013,0.106);
  yproj_slope->SetAxisRange(-0.25,0.25,"X");
  yproj_slope->SetYTitle("Counts per 6 keV");
  yproj_slope->GetYaxis()->CenterTitle(1);
  yproj_slope->SetXTitle("Excitation Energy");
  yproj_slope->GetXaxis()->CenterTitle(1);
  yproj_slope->SetStats(kFALSE);
  yproj_slope->SetTitle();

  cFit->Close();
  mkCanvas2("csimQ","csimQ");
  csimQ->SetTopMargin(.02);
  csimQ->SetRightMargin(.03);  
  csimQ->SetCanvasSize(width,(UInt_t)(width*ratio1));

  yproj_slope->SetLineWidth(2);
  yproj_slope->Draw();
  hmap2->SetLineWidth(2);
  hmap2->SetLineColor(2);
  hmap2->Draw("same");
  csimQ->SaveAs("csimQ.eps");
}

void simplots3()
{
  gStyle->SetOptDate(0);
  plotsim();
  filltree("simulations/21res84tgt2map_450.dat","tree1","simulations/492_param.dat");
  shiftx2("hEZg",-hp[14]-450);
  hEZg_shift->Clone("hEZg_450");
  hEZg->Reset();
  
  filltree("simulations/21res84tgt2map_100.dat","tree2","simulations/094_param.dat");
  shiftx2("hEZg",-hp[14]-100);
  add2("hEZg_shift","hEZg_450");
  TH2F * hOutput1=(TH2F *) gROOT->FindObject("hEZg_shift_copy");
  
  hOutput1->SetAxisRange(plot_minZ,plot_maxZ,"X");
  hOutput1->SetAxisRange(minE,maxE,"Y");
  hOutput1->SetYTitle(ytitle);
  hOutput1->GetYaxis()->CenterTitle(1);
  hOutput1->SetXTitle(xtitle);
  hOutput1->GetXaxis()->CenterTitle(1);
  hOutput1->SetStats(kFALSE);
  hOutput1->SetTitle();

  cFit->Close();
  mkCanvas2("csimB","csimB");
  csimB->SetTopMargin(.02);
  csimB->SetRightMargin(.03);  
  csimB->SetCanvasSize(width,(UInt_t)(width*ratio1));

  hOutput1->Draw("col");
  plotlines("_lineB.txt",7,2);
  plotvlines(-100,0,-450);

  csimB->SaveAs("csimB.eps");
}

void timeplot3()
{
  gStyle->SetOptDate(0);
  loadfiles();
  
  _filename0->cd();
  slopexy("hET",1.076,-2.395);
  TH2F * hOutput1=(TH2F *) gROOT->FindObject("hET_slope");
  cFit->Close();
  
  hOutput1->SetAxisRange(0,82,"X");
  hOutput1->SetAxisRange(minE,maxE,"Y");
  hOutput1->SetYTitle(ytitle);
  hOutput1->GetYaxis()->CenterTitle(1);
  hOutput1->SetXTitle("Time (ns)");
  hOutput1->GetXaxis()->CenterTitle(1);
  hOutput1->SetStats(kFALSE);
  hOutput1->SetTitle();
  hOutput1->SetMinimum(32);

  mkCanvas2("cET_big","cET_big");
  cET_big->SetTopMargin(.02);
  cET_big->SetRightMargin(.03);  
  cET_big->SetCanvasSize(width,(UInt_t)(width*ratio2));
  cET_big->SetLogz();
  hOutput1->Draw("col");

  _filename3=new TFile("cut3.root");//missing!
  CUTG->SetLineStyle(2);
  CUTG->SetLineWidth(2);
  CUTG->Draw("same");

  cET_big->SaveAs("cET_big.eps");
}

//-----------------------------------------------------------------------------------------------
//macros for checking and setting calibration for EMMA PGAC test data----------------------------
void showX(float ymin=-100,float ymax=100)
{//test extent of x-positions (calibrated)
  plotvlines(0.001-6,6.001-6);
  plotvlines(166-6,160-6);
}

void showY()
{//test extent of y-positions (calibrated) and show anode segments
  //setvlines(-300,100);
  plotvlines(0.001-6,6.001-6);
  plotvlines(22-6,0,0,0,2);
  plotvlines(44-6,0,0,0,2);
  plotvlines(66-6);
  plotvlines(0,60-6);
}

void showXY(float xmin=55,float xmax=180,float ymin=-100,float ymax=100)
{
  setvlines(ymin,ymax);
  showX();
  plothlines(0.001-6,6.001-6);
  plothlines(22-6,0,0,xmin,xmax,2);
  plothlines(44-6,0,0,xmin,xmax,2);
  plothlines(66-6);
  plothlines(0,60-6);
  /*  
      TEllipse *ellipse = new TEllipse(118,54/2,99.5/2.);
      ellipse->SetLineWidth(2);
      ellipse->SetLineStyle(4); 
      ellipse->SetFillStyle(0);
      ellipse->Draw();
  */
}

void showX1(float ymin=-100,float ymax=100)
{
  setvlines(ymin,ymax);
  showX();
  loadcal("cal/X1.lst");
  plotvlines(0,0,sorted[0]);
  plotvlines(0,0,sorted[1]);
  plotvlines(0,0,sorted[2]);
  plotvlines(0,0,sorted[3]);
  plotvlines(0,0,sorted[4]);
}

void showY1x(float ymin=-100,float ymax=100)
{
  setvlines(ymin,ymax);  
  showY();
  loadcal("cal/Y1.lst");  
  plotvlines(0,0,sorted[1]);
  plotvlines(0,0,sorted[2]);
  plotvlines(0,0,sorted[3]);
  plotvlines(0,0,sorted[4]);
}

void showY1y(float xmin=55,float xmax=180,float ymin=-100,float ymax=100)
{
  showXY();
  loadcal("cal/Y1.lst");
  setvlines(ymin,ymax);  
  plothlines(0,0,sorted[1],xmin,xmax);
  plothlines(0,0,sorted[2],xmin,xmax);
  plothlines(0,0,sorted[3],xmin,xmax);
  plothlines(0,0,sorted[4],xmin,xmax);
  showX1();
}

void showX2(float ymin=-100,float ymax=100)
{
  setvlines(ymin,ymax);  
  showX();
  loadcal("cal/X2.lst");
  plotvlines(0,0,sorted[0]);
  plotvlines(0,0,sorted[1]);
  plotvlines(0,0,sorted[2]);
  plotvlines(0,0,sorted[3]);
  plotvlines(0,0,sorted[4]);
}

void showY2x(float ymin=-100,float ymax=100)
{
  setvlines(ymin,ymax);  
  showY();
  loadcal("cal/Y2.lst");  
  plotvlines(0,0,sorted[1]);
  plotvlines(0,0,sorted[2]);
  plotvlines(0,0,sorted[3]);
  plotvlines(0,0,sorted[4]);
}

void showY2y(float xmin=55,float xmax=180,float ymin=-100,float ymax=100)
{
  showXY();
  loadcal("cal/Y2.lst");  
  plothlines(0,0,sorted[1],xmin,xmax);
  plothlines(0,0,sorted[2],xmin,xmax);
  plothlines(0,0,sorted[3],xmin,xmax);
  plothlines(0,0,sorted[4],xmin,xmax);
  showX2(ymin,ymax);
}

static const int cpeaks=6;  
Float_t cpositions[cpeaks]={};
Float_t sorted[cpeaks]={};

void loadcal(Char_t *filename="cal/X1.lst")
{
  printf("Loading file %s...\n",filename);
  Float_t max=0;
  ifstream infile(filename);
  int i =0;  
  while (infile >> cpositions[i]) {
    printf(" %13.6f ",cpositions[i]);
    i++;
  }
  printf("\n");
    
  printf("Sorting points... ");
  sorted[0]=cpositions[0];//initializes sorted array with a valid position
  for(Int_t i=0;i<cpeaks;i++){
    if((cpositions[i])<(sorted[0])){
      sorted[0]=cpositions[i];//locates minimum
    }
    if((cpositions[i])>max){
      max=cpositions[i];//locates maximum
    }
  }
  for(Int_t i=1;i<cpeaks;i++){
    sorted[i]=max;    
    for(Int_t j=0;j<cpeaks;j++){
      if(((cpositions[j])<=(sorted[i]))&&((cpositions[j])>(sorted[i-1]))){
	sorted[i]=cpositions[j];//locates next-smallest position
      }
    }
  }
  printf("sort complete.\n");
  for(int i=0;i<cpeaks;i++){
    cpositions[i]=sorted[i];
    printf(" %13.6f ",cpositions[i]);
  }
  printf("\n");
}

TF1 *ls1,*ls2,*ls3,*lstotal;

void funcload(Char_t *filename="cal/X1.lst",Float_t slope=1, Float_t offset=0)
{//loads functions for funcfit() base on ranges in cal file
  TF1 *fnew = new TF1("fnew","(x>0.46)*(x<0.463)*pol1+(x>0.467)*(x<0.492)*pol1+(x>0.497)*(x<0.505)*pol1");
  
  loadcal(filename);

  ls1 = new TF1("ls1","pol3",(cpositions[0]-offset)/slope,(cpositions[1]-offset)/slope);
  ls2 = new TF1("ls2","pol3",(cpositions[2]-offset)/slope,(cpositions[3]-offset)/slope);
  ls3 = new TF1("ls3","pol3",(cpositions[4]-offset)/slope,(cpositions[5]-offset)/slope);
  lstotal = new TF1("lstotal","pol3(0)+pol3(4)+pol3(7)",(cpositions[0]-offset)/slope,(cpositions[5]-offset)/slope);

  lstotal->SetLineColor(2);

  TString name = "ls";
  for(int i=0;i<3;i++){
    name = "ls";
    name+=i+1;
    TF1 * finput=gROOT->FindObject(name.Data());

    printf("%s range is %f %f\n",name.Data(),finput->GetXaxis()->GetXmin(),finput->GetXaxis()->GetXmax());
  }
}

void funcfit(Char_t *histin="hposc0", Char_t *filename="cal/X1.lst",Float_t slope=1, Float_t offset=0)
{//based on multifit(); fits histogram over ranges defined in loadfunc
  hname=histin;
  funcload(filename,slope,offset);
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  cFit->Clear();
  hProj=(TH1F *) gROOT->FindObject(histin); 
  hProj->Draw("COL");
  
  Double_t par[12];
  Double_t chi[4];
  Double_t ndf[4];
  
  // Fit each function and add it to the list of functions 
  hProj->Fit(ls1,"R");//the option of the first function must not include "+" so that the fuctions are reset
  chi[1]=ls1->GetChisquare();
  ndf[1]=ls1->GetNDF();
  hProj->Fit(ls2,"R+");
  chi[1]=ls2->GetChisquare();
  ndf[1]=ls2->GetNDF();
  hProj->Fit(ls3,"R+");
  chi[1]=ls3->GetChisquare();
  ndf[1]=ls3->GetNDF();

  float chi_sum=0;
  float ndf_sum=0;
  TString name = "ls";
  for(int i=0;i<3;i++){
    name = "ls";
    name+=i+1;
    TF1 * finput=gROOT->FindObject(name.Data());
    chi[i]=finput->GetChisquare();
    ndf[i]=finput->GetNDF();
    chi_sum+=chi[i];
    ndf_sum+=(chi[i]/ndf[i]);
    printf("%s: Chi squared = %f, Ndf = %f, (Chi squared)/Ndf = %f\n",name.Data(),chi[i],ndf[i],chi[i]/ndf[i]);  
    printf("Chi squared sum is %f, (Chi squared)/Ndf sum is %f\n",chi_sum,ndf_sum);
  }
  // Get the parameters from the fit 
  ls1->GetParameters(&par[0]); 
  ls2->GetParameters(&par[4]); 
  ls3->GetParameters(&par[7]); 
 
  // Use the parameters on the sum 
  lstotal->SetParameters(par); 
  hProj->Fit(lstotal,"R+");
}

TF1 *ruth;
void ruthdef()
{//Rutherford scattering cross section function with offsets and scale
  //  ruth =new TF1("ruth","[0]+[1]*TMath::Power(TMath::Sin([2]*(x*TMath::DegToRad()/2)-[3]),-4)");
  ruth =new TF1("ruth","[0]*TMath::Power(TMath::Sin(((x)*TMath::DegToRad())/2),-4)");
  //ruth->SetParLimits(1,0,1e+06);
  //ruth->SetParLimits(2,-1e-01,1e-01);
  ruth->SetParName(0,"y scale");  
  //ruth->SetParName(1,"y offset");//physical meaning?
  //ruth->SetParName(2,"x scale");//does't work  
  //ruth->SetParName(1,"x offset");  
  //  ruth->SetParameters(1,1);
}

void ruthtest()
{
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();    
  if(gROOT->FindObject("hruth"))hruth->Delete();
  h1("hruth","Rutherford Scattering Data",1000,1,180);
  fillhist("ruth.lst","hruth",0);
  ruthdef();
  hruth->Fit("ruth","m","",10,160);
  cFit->SetLogy();
  hruth->SetMarkerStyle(2);
  hruth->SetMarkerColor(1);
  hruth->SetMarkerSize(2);
  hruth->Draw("P");
}

TF1 *ruther;
Float_t ruther_min=1e-1;
void rutherdef(Float_t min=ruther_min, Float_t max=180)
{//Rutherford scattering cross section function (simple)
  ruther =new TF1("ruther","[0]*TMath::Power(TMath::Sin(TMath::DegToRad()*x/2),-4)",1e-1,180);
  ruther->SetParLimits(0,0,1e+05);
  ruther->SetParName(0,"y scale");  
  ruther->SetParameter(0,1);
  ruther->SetRange(min,max);
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  cFit->SetLogy();  
  ruther->Draw();
}

void rutherfill(Int_t nfill=1e4, Char_t *histin="h1")
{
  if(!(gROOT->FindObject("ruther")))
    rutherdef();
  hname=histin;
  if(!gROOT->FindObject(histin))
    mkhist1(hname,360,180);
  hProj=(TH1F *) gROOT->FindObject(histin); 
  for(int i=0;i<nfill;i++){
    hProj->Fill(ruther->GetRandom());
  }
  hProj->Fit("ruther","q","",1,180);
}

void rlindef(Float_t min=ruther_min, Float_t max=180)
{//Rutherford scattering cross section function (simple)
  rlin =new TF1("rlin","[0]+[1]*x",1e-1,180);
  rlin->SetParName(0,"offset");  
  rlin->SetParName(1,"slope");  
  rlin->SetParameter(0,1);
  rlin->SetRange(min,max);
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  cFit->SetLogy();  
  rlin->Draw();
}

void ruthset()
{
  ruthdef();
  float sets[4]={-1.87e+02,6.41e-03,6.60e-05,-5.92e-02};
  printf("Setting ruth with default parameters...\n");
  for(int i=0;i<4;i++){
    ruth->SetParameter(i,sets[i]);
    printf(" parameter %d = %g\n",i,ruth->GetParameter(i));
  }
}

void ruthfitsimp(Char_t *histin="hposc0",Float_t min=101,Float_t max=148)
{
  ruthset();
  hname=histin;
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  cFit->Clear();
  hProj=(TH1F *) gROOT->FindObject(histin); 
  hProj->Draw("COL");
  hProj->Fit("ruth","","",min,max);
}

void ruthload(Char_t *filename="cal/X1.lst",Float_t slope=1, Float_t offset=0)
{
  loadcal(filename);
  ruthset();
  
  ls1 = new TF1("ls1","ruth",(cpositions[0]-offset)/slope,(cpositions[1]-offset)/slope);
  ls2 = new TF1("ls2","ruth",(cpositions[2]-offset)/slope,(cpositions[3]-offset)/slope);
  ls3 = new TF1("ls3","ruth",(cpositions[4]-offset)/slope,(cpositions[5]-offset)/slope);
  lstotal = new TF1("lstotal","ruth",(cpositions[0]-offset)/slope,(cpositions[5]-offset)/slope);

  lstotal->SetLineColor(2);

  TString name = "ls";
  for(int i=0;i<3;i++){
    name = "ls";
    name+=i+1;
    TF1 * finput=gROOT->FindObject(name.Data());

    printf("%s range is %f %f\n",name.Data(),finput->GetXaxis()->GetXmin(),finput->GetXaxis()->GetXmax());
  }
}

void ruthfit(Char_t *histin="hposc0", Char_t *filename="cal/X1.lst",Float_t slope=1, Float_t offset=0)
{//calls ruthload() and follows multifit
  hname=histin;
  ruthload(filename,slope,offset);
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  cFit->Clear();
  hProj=(TH1F *) gROOT->FindObject(histin); 
  hProj->Draw("COL");
  
  Double_t par[12];
  Double_t chi[4];
  Double_t ndf[4];
  
  // Fit each function and add it to the list of functions 
  hProj->Fit(ls1,"BR");
  chi[1]=ls1->GetChisquare();
  ndf[1]=ls1->GetNDF();
  ls1->GetParameters(&par[0]); 
  ls2->SetParameters(par); 
  //  ls3->SetParameters(par); 
  hProj->Fit(ls2,"BR+");
  chi[1]=ls2->GetChisquare();
  ndf[1]=ls2->GetNDF();
 
  hProj->Fit(ls3,"BR+");
  chi[1]=ls3->GetChisquare();
  ndf[1]=ls3->GetNDF();

  float chi_sum=0;
  float ndf_sum=0;
  TString name = "ls";
  for(int i=0;i<3;i++){
    name = "ls";
    name+=i+1;
    TF1 * finput=gROOT->FindObject(name.Data());
    chi[i]=finput->GetChisquare();
    ndf[i]=finput->GetNDF();
    chi_sum+=chi[i];
    ndf_sum+=(chi[i]/ndf[i]);
    printf("%s: Chi squared = %f, Ndf = %f, (Chi squared)/Ndf = %f\n",name.Data(),chi[i],ndf[i],chi[i]/ndf[i]);  
    printf("Chi squared sum is %f, (Chi squared)/Ndf sum is %f\n",chi_sum,ndf_sum);
  }
  // Get the parameters from the fit 
 
  ls2->GetParameters(&par[4]); 
  ls3->GetParameters(&par[7]); 
 
  // Use the parameters on the sum 
  lstotal->SetParameters(par);
  hProj->Fit(lstotal,"BR+");
}

void ruthload2(int pos=0, Char_t *filename="cal/X1.lst",Float_t slope=1, Float_t offset=0)
{//testing variable range
  loadcal(filename);
  ruthset();
  
  ls1 = new TF1("ls1","ruth * (x>(([0]-[3])/[2]))*(x<(([1]-[3])/[2]))");
  ls1->SetParameters(101,148,1,10,-1.87e+02,6.41e-03,6.60e-05,-5.92e-02);//,101,148); 
  int off=4;//number of parameters before ruth is called
  ls1->SetParLimits(2,0.95,1.05);
  ls1->SetParLimits(off+1,1e-05,1e+05);
  ls1->SetParLimits(off+2,0,1e+05);

  for(int i=0;i<2;i++){
    ls1->FixParameter(i,cpositions[i+pos*2]);
    printf("i=%d cpositions[%d]=%f parameter[%d]=%f\n",i,i,cpositions[i],i,ls1->GetParameter(i));
  }
}

void lsload(Char_t *filename="cal/X1.lst",Float_t slope=1, Float_t offset=0)
{//defines fuction ls1 and fixes parameters
  loadcal(filename);
     
  ls1 = new TF1("ls1","([8]+[9]*x+[10]*x*x)*(x>(([0]-[7])/[6]))*(x<(([1]-[7])/[6]))+([8]+[9]*x+[10]*x*x)*(x>(([2]-[7])/[6]))*(x<(([3]-[7])/[6]))+([8]+[9]*x+[10]*x*x)*(x>(([4]-[7])/[6]))*(x<(([5]-[7])/[6]))");
  ls1->SetParName(6,"p6 range slope");
  ls1->SetParName(7,"p7 range offset");
  ls1->SetParName(8,"p8 constant");
  ls1->SetParName(9,"p9 linear");
  ls1->SetParName(10,"p10 quadratic");
  printf("Fixing parameters to match positions from %s...\n",filename);
  for(int i=0;i<6;i++){
    ls1->FixParameter(i,cpositions[i]);
    printf(" cpositions[%d]=%5.1f => parameter[%d]=%5.1f\n",i,cpositions[i],i,ls1->GetParameter(i));
  }
  ls1->SetParameter(6,slope);
  ls1->SetParameter(7,offset);
}

void ls2load(Char_t *filename="cal/X1.lst",Float_t slope=1, Float_t offset=0)
{//defines fuction ls2 and fixes parameters
  loadcal(filename);
     
  ls2 = new TF1("ls2","([8]+[9]*x+[10]*x*x)*(x>(([0]-[7])/[6]))*(x<(([1]-[7])/[6]))+([11]+[12]*x+[13]*x*x)*(x>(([2]-[7])/[6]))*(x<(([3]-[7])/[6]))+([14]+[15]*x+[16]*x*x)*(x>(([4]-[7])/[6]))*(x<(([5]-[7])/[6]))");
  ls2->SetParName(6,"p6 range slope");
  ls2->SetParName(7,"p7 range offset");
  ls2->SetParName(8,"p8 constant");
  ls2->SetParName(9,"p9 linear");
  ls2->SetParName(10,"p10 quadratic");

  for(int i=0;i<6;i++){
    ls2->FixParameter(i,cpositions[i]);
    printf("i=%d cpositions[%d]=%5.1f parameter[%d]=%5.1f\n",i,i,cpositions[i],i,ls1->GetParameter(i));
  }
  ls2->SetParameter(6,slope);
  ls2->SetParameter(7,offset);
}

void lssetc()
{//set (temporary) parameters to aid fit for calibrated spectrum
  lsload();
  ls1->SetParameter(6,1);
  ls1->SetParameter(7,0);
  ls1->SetParameter(8,0);
  ls1->SetParameter(9,100);
}

void lsset()
{//set (temporary) parameters to aid fit
  lsload();
  ls1->SetParameter(6,1914);
  ls1->SetParameter(7,-882);
  ls1->SetParameter(8,0);
  ls1->SetParameter(9,100);
  ls1->SetParameter(10,0);
  //ls1->SetParLimits(6,);
  //ls1->SetParLimits(7,);
}

void lsreset()
{//set (temporary) parameters to aid fit
 // lsload();
  ls1->SetParameter(6,617);
  ls1->SetParameter(7,-170);
}

void lsfix()
{//fix fit parameters (before adjusting the range parameters)
  double temp[3]={};
  for(int i=0;i<3;i++){
    temp[i]=ls1->GetParameter(i+8);
    ls1->FixParameter(i+8,temp[i]);
  }
  for(int i=0;i<11;i++){
    printf("  %2d  %-20s %10.2f\n",i,ls1->GetParName(i),ls1->GetParameter(i));
  }
}

void lsfixlin()
{//fix fit parameters (before adjusting the range parameters)
  double temp[3]={};
  for(int i=0;i<3;i++){
    temp[i]=pol1->GetParameter(i);
    ls1->FixParameter(i+8,temp[i]);
  }
  for(int i=0;i<11;i++){
    printf("  %2d  %-20s %10.2f\n",i,ls1->GetParName(i),ls1->GetParameter(i));
  }
}

void lsfixquad()
{//fix fit parameters (before adjusting the range parameters)
  double temp[3]={};
  for(int i=0;i<3;i++){
    temp[i]=pol2->GetParameter(i);
    ls1->FixParameter(i+8,temp[i]);
  }
  for(int i=0;i<11;i++){
    printf("  %2d  %-20s %10.2f\n",i,ls1->GetParName(i),ls1->GetParameter(i));
  }
}

void lssetquad()
{//fix fit parameters (before adjusting the range parameters)
  double temp[3]={};
  for(int i=0;i<3;i++){
    temp[i]=pol2->GetParameter(i);
    ls1->SetParameter(i+8,temp[i]);
  }
  for(int i=0;i<11;i++){
    printf("  %2d  %-20s %10.2f\n",i,ls1->GetParName(i),ls1->GetParameter(i));
  }
}

void ruthload3(Char_t *filename="cal/X1.lst",Float_t slope=1, Float_t offset=0)
{
  loadcal(filename);
  ruthset();
    
  ls1 = new TF1("ls1","([8]+[9]*TMath::Power(sin([10]*x-[11]),-4))*(x>(([0]-[7])/[6]))*(x<(([1]-[7])/[6]))+([8]+[9]*TMath::Power(sin([10]*x-[11]),-4))*(x>(([2]-[7])/[6]))*(x<(([3]-[7])/[6]))+([8]+[9]*TMath::Power(sin([10]*x-[11]),-4))*(x>(([4]-[7])/[6]))*(x<(([5]-[7])/[6]))");
  ls1->SetParameters(67,88,94,145,150,154,slope,offset,-1.87e+02,6.41e-03,6.60e-05);//,-5.92e-02);
  ls1->SetParLimits(9,0,1e+05);
  ls1->SetParLimits(10,0,1e+05);
  ls1->SetParName(6,"p6 range slope");
  ls1->SetParName(7,"p7 range offset");
  ls1->SetParName(8,"fit y-offset");
  ls1->SetParName(9,"fit y-scale");
  ls1->SetParName(10,"fit x-scale");
  ls1->SetParName(11,"fit x-offset");

  for(int i=0;i<6;i++){
    ls1->FixParameter(i,cpositions[i]);
    printf("i=%d cpositions[%d]=%5.1f parameter[%d]=%5.1f\n",i,i,cpositions[i],i,ls1->GetParameter(i));
  }
}

void ruthfit3(Char_t *histin="hposc0", Char_t *filename="cal/X1.lst",Float_t slope=1, Float_t offset=0)
{
  hname=histin;
  ruthload3(filename,slope,offset);
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  cFit->Clear();
  hProj=(TH1F *) gROOT->FindObject(histin); 
  hProj->Draw("COL");

  hProj->Fit("ls1","MB","");
}

void qfitc(Char_t *histname,Float_t center=0, Float_t width=10)
{//fit a quadratic function about a central point
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
  TH1F *hist1=(TH1F*) gROOT->FindObject(histname);
  hist1->Draw();
  Float_t xmin=center-width;
  Float_t xmax=center+width;
  hist1->Fit("pol2","","",xmin,xmax);
  Float_t a=0,b=0,c=0;
  double temp[3]={};
  for(int i=0;i<3;i++){
    temp[i]=pol2->GetParameter(i);
  }
  c=temp[0];
  b=temp[1];
  a=temp[2];
  printf("Fit range is %f to %f\n",xmin,xmax);
  printf("Fit range center is %f\n",center);
  printf("     Fit minimum is %f\n",-b/(2*a));
}

void qfitc2(Char_t *histname,Float_t center1=0, Float_t center2=0,Float_t width=0.0021)
{//fit two quadratic functions, relative to two central points
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
  cFit->Clear();
  TH1F *hist1=(TH1F*) gROOT->FindObject(histname);
  Float_t xmin=center1-width;
  Float_t xmax=center1+width;
  hist1->Fit("pol2","q","",xmin,xmax);
  hist1->GetFunction("pol2")->SetLineColor(2);
  hist1->GetFunction("pol2")->SetLineStyle(1);
  pol2->SetLineColor(2);
  pol2->SetLineStyle(1);
  Float_t a=0,b=0,c=0;
  double temp[12]={};
  for(int i=0;i<3;i++){
    temp[i]=pol2->GetParameter(i);
  }
  c=temp[0];
  b=temp[1];
  a=temp[2];
  double min[4]={};
  min[0]=-b/(2*a);
  printf("Fit range is %f to %f\n",xmin,xmax);
  printf("Fit range center is %f\n",center1);
  printf("     Fit minimum is %f\n",min[0]);

  xmin=center2-width;
  xmax=center2+width;
  hist1->Fit("pol2","+q","",xmin,xmax);
 
  for(int i=0;i<3;i++){
    temp[i+3]=pol2->GetParameter(i);
  }
  c=temp[3];
  b=temp[4];
  a=temp[5];
  min[1]=-b/(2*a);
  printf("Fit range is %f to %f\n",xmin,xmax);
  printf("Fit range center is %f\n",center2);
  printf("     Fit minimum is %f\n",min[1]);
}

double mid[2]={};

void loadcalm(Char_t *filename="cal/X1.lst")
{//load cal and find midpoints of masks
  loadcal(filename);

  printf("The midpoint of the mask projections are:\n"); 
  for(int i=0;i<2;i++){
    mid[i]=(cpositions[(i+1)*2]-cpositions[(i+1)*2-1])/2.;
    printf(" The midpoint of %13.6f and %13.6f is (%13.6f wide)",cpositions[(i+1)*2-1],cpositions[(i+1)*2],mid[i]);
    mid[i]+=cpositions[(i+1)*2-1];
    printf(" at %13.6f\n",mid[i]);
  }
  printf(" The difference in the midpoints is       %f\n",mid[1]-mid[0]);
}

void qfitc2m(Char_t *histname,Float_t center1=0, Float_t center2=0,Float_t width=0.0021,Char_t *filename="cal/X1.lst")
{//fit two quadratic functions, relative to two central points & set calibration scale based on mask position from file
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
  cFit->Clear();
  TH1F *hist1=(TH1F*) gROOT->FindObject(histname);
  Float_t xmin=center1-width;
  Float_t xmax=center1+width;
  hist1->Fit("pol2","q","",xmin,xmax);  
  hist1->GetFunction("pol2")->SetLineColor(2);
  hist1->GetFunction("pol2")->SetLineStyle(1);
  pol2->SetLineColor(2);
  pol2->SetLineStyle(1);
    
  Float_t a=0,b=0,c=0;
  double temp[12]={};
  
  for(int j=0;j<3;j++) {//iterate center-finding
    //calcualte fit range based on center and width
    xmin=center1-width;
    xmax=center1+width;
    hist1->Fit("pol2","q","",xmin,xmax);  

    for(int i=0;i<3;i++){//read in fit parameters and calculate min
      temp[i]=pol2->GetParameter(i);
    }
    c=temp[0];
    b=temp[1];
    a=temp[2];
    double min[4]={};
    min[0]=-b/(2*a);
    printf("Fit range is %f to %f\n",xmin,xmax);
    printf("Fit range center is %f\n",center1);
    printf("     Fit minimum is %f (%f)\n",min[0],center1-min[0]);
    
    //calcualte fit range based on center and width
    xmin=center2-width;
    xmax=center2+width;
    hist1->Fit("pol2","+q","",xmin,xmax);
    
    for(int i=0;i<3;i++){//read in fit parameters and calculate min
      temp[i+3]=pol2->GetParameter(i);
    }
    c=temp[3];
    b=temp[4];
    a=temp[5];
    min[1]=-b/(2*a);
    printf("Fit range is %f to %f\n",xmin,xmax);
    printf("Fit range center is %f\n",center2);
    printf("     Fit minimum is %f      (%f)\n",min[1],center2-min[1]);
    
    center1=min[0];
    center2=min[1];
  }

  printf("Histogram is %s with title %s.\n",histname,hist1->GetTitle());
  hname=histname;
  for(Int_t i=0;i<hname.Length();i++){//loop added by Jack
    TString tempst="";
    tempst=hname(i,hname.Length()-i);
    if(tempst.IsFloat())
      {
	det=tempst.Atoi();
	break;
      }
  }
  printf("Detector number is %d\n",det);
  switch(det){
  case 0:
    filename="cal/X1.lst";
    break;
  case 1:
    filename="cal/Y1.lst";
    break;
  case 2:
    filename="cal/X2.lst";
    break;
  case 3:
    filename="cal/Y1.lst";
    break;
  default:
    printf("Detector %d not recognized!\n",det);
    
    break; 
  }
  //printf("Loading calibration file %s...\n",filename);
  loadcalm(filename);  
  double slope=(mid[1]-mid[0])/(min[1]-min[0]);
  printf(" The difference in the measured minima is %f\n",min[1]-min[0]);
  printf(" slope is %f\n",slope);
  double offset=slope*-min[0]+mid[0];
  printf(" offset is %f\n",offset);
  
  printf(" [y=m*x+b  ] Fit parameters are:         Slope = %8.3f, Offset = %8.3f\n",slope,offset);
  printf(" [x=(y-b)/m] Inverse fit parameters are: Slope = %8.3f, Offset = %8.3f\n",1/slope,-offset/slope);
  
  FILE * outfile;
  outfile=fopen("temp.lst","w");
  fprintf(outfile,"%g, %g\n",slope,offset);
  fclose(outfile);

  double edge[6]={};
  for(int i=0;i<6;i++){
    edge[i]=(cpositions[i]-offset)/slope;
    //printf("Edge %d estimated location is %f\n",i,edge[i]);
  }
  // printf("estimated edges are located at %f, %f\n",(edge[1]+edge[2])/2,(edge[3]+edge[4])/2);
}

void qfitc2mr(Char_t *histname,Float_t center1=0, Float_t center2=0,Float_t width=0.0021,Char_t *filename="cal/X1.lst")
{//fit two quadratic functions, relative to two central points & set calibration scale based on mask position from file
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
  cFit->Clear();
  TH1F *hist1=(TH1F*) gROOT->FindObject(histname);
  Float_t xmin=center1-width;
  Float_t xmax=center1+width;
  hist1->Fit("pol2","q","",xmin,xmax);  
  hist1->GetFunction("pol2")->SetLineColor(2);
  hist1->GetFunction("pol2")->SetLineStyle(1);
  pol2->SetLineColor(2);
  pol2->SetLineStyle(1);
  
  
  Float_t a=0,b=0,c=0;
  double temp[12]={};
  for(int i=0;i<3;i++){
    temp[i]=pol2->GetParameter(i);
  }
  c=temp[0];
  b=temp[1];
  a=temp[2];
  double min[4]={};
  min[0]=-b/(2*a);
  printf("Fit range is %f to %f\n",xmin,xmax);
  printf("Fit range center is %f\n",center1);
  printf("     Fit minimum is %f\n",min[0]);

  xmin=center2-width;
  xmax=center2+width;
  hist1->Fit("pol2","+q","",xmin,xmax);
 
  for(int i=0;i<3;i++){
    temp[i+3]=pol2->GetParameter(i);
  }
  c=temp[3];
  b=temp[4];
  a=temp[5];
  min[1]=-b/(2*a);
  printf("Fit range is %f to %f\n",xmin,xmax);
  printf("Fit range center is %f\n",center2);
  printf("     Fit minimum is %f\n",min[1]);

  printf("Hisogram is %s with title %s.\n",histname,hist1->GetTitle());
  loadcalm(filename);
  hname=histname;
  for(Int_t i=0;i<hname.Length();i++){//loop added by Jack
    TString tempst="";
    tempst=hname(i,hname.Length()-i);
    if(tempst.IsFloat())
      {
	det=tempst.Atoi();
	break;
      }
  }
  printf("Detector number is %d\n",det);
  double slope=(mid[1]-mid[0])/(min[1]-min[0]);
  printf(" The difference in the measured minima is %f\n",min[1]-min[0]);
  printf(" The midpoint between the measured minima is %f\n",(min[1]-min[0])/2+min[0]);
  printf(" slope is %f\n",slope);
  double offset=slope*-min[0]+mid[0];
  printf(" offset is %f\n",offset);
  
  printf(" [y=m*x+b  ] Fit parameters are:         Slope = %8.3f, Offset = %8.3f\n",slope,offset);
  printf(" [x=(y-b)/m] Inverse fit parameters are: Slope = %8.3f, Offset = %8.3f\n",1/slope,-offset/slope);
  
  FILE * outfile;
  outfile=fopen("temp.lst","w");
  fprintf(outfile,"%g, %g\n",slope,offset);
  fclose(outfile);

  double edge[6]={};
  for(int i=0;i<6;i++){
    edge[i]=(cpositions[i]-offset)/slope;
  }
  
  gaus->SetLineColor(1);
  if(det%2){//y-positions
    for(i=0;i<6;i++){
      printf("Edge %d estimated location is %7.3f",i,edge[i]);
      if(i%2)//right-hand sides
	hist1->Fit("gaus","q+","",edge[i],edge[i]+0.85*width);
      else//left-hand sides
	hist1->Fit("gaus","q+","",edge[i]-width*0.85,edge[i]);
      temp[i+6]=gaus->GetParameter(2);
      printf(", width is %5.3f\n",temp[i+6]);
    }
  }
  else{//x-positions
    for(i=0;i<6;i++){
      printf("Edge %d estimated location is %7.3f",i,edge[i]);
      if(i%2)//right-hand sides
	hist1->Fit("gaus","q+","",edge[i]-width/3,edge[i]+width/2);
      else//left-hand sides
	hist1->Fit("gaus","q+","",edge[i]-width/2,edge[i]+width/3);
      temp[i+6]=gaus->GetParameter(2);
      printf(", width is %5.3f\n",temp[i+6]);
    }
  }
  double temp_sum=0;
  for(i=1;i<5;i++){
    temp_sum+=temp[i+6];
    //printf("width is %5.3f, temp sum is %f\n",temp[i+6],temp_sum);
  }
  printf("Average width of central edges is %f mm (%f mm FWHM)\n",temp_sum/4,temp_sum/4*2.35482);
}

void lsgetslope()
{
  double slope=0,offset=0;
  slope=ls1->GetParameter(6);
  offset=ls1->GetParameter(7);
  printf("[y=m*x+b  ] Fit parameters are:         Slope = %8.3f, Offset = %8.3f\n",slope,offset);
  printf("[x=(y-b)/m] Inverse fit parameters are: Slope = %8.3f, Offset = %8.3f\n",1/slope,-offset/slope);
  
  FILE * outfile;
  outfile=fopen("temp.lst","w");
  fprintf(outfile,"%g, %g\n",slope,offset);
  fclose(outfile);
  outfile=fopen("temp_inv.lst","w");
  fprintf(outfile,"%g, %g\n",1/slope,-offset/slope);
  fclose(outfile);
}

void qfitc2mls(Char_t *histname,Float_t center1=0, Float_t center2=0,Float_t width=0.0021,Char_t *filename="cal/X1.lst")
{//fit two quadratic functions, relative to two central points & set calibration scale based on mask position from file, then calculated resolution based on linear segment fit.
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
  cFit->Clear();
  cFit->Divide(1,2);
  cFit->cd(1);
  TH1F *hist1=(TH1F*) gROOT->FindObject(histname);
  Float_t xmin=center1-width;
  Float_t xmax=center1+width;
  hist1->Fit("pol2","q","",xmin,xmax);  
  hist1->GetFunction("pol2")->SetLineColor(2);
  hist1->GetFunction("pol2")->SetLineStyle(1);
  pol2->SetLineColor(2);
  pol2->SetLineStyle(1);
  
  
  Float_t a=0,b=0,c=0;
  double temp[12]={};
  for(int i=0;i<3;i++){
    temp[i]=pol2->GetParameter(i);
  }
  c=temp[0];
  b=temp[1];
  a=temp[2];
  double min[4]={};
  min[0]=-b/(2*a);
  printf("Fit range is %f to %f\n",xmin,xmax);
  printf("Fit range center is %f\n",center1);
  printf("     Fit minimum is %f\n",min[0]);

  xmin=center2-width;
  xmax=center2+width;
  hist1->Fit("pol2","+q","",xmin,xmax);
 
  for(int i=0;i<3;i++){
    temp[i+3]=pol2->GetParameter(i);
  }
  c=temp[3];
  b=temp[4];
  a=temp[5];
  min[1]=-b/(2*a);
  printf("Fit range is %f to %f\n",xmin,xmax);
  printf("Fit range center is %f\n",center2);
  printf("     Fit minimum is %f\n",min[1]);

  printf("Hisogram is %s with title %s.\n",histname,hist1->GetTitle());
  loadcalm(filename);
  double slope=(mid[1]-mid[0])/(min[1]-min[0]);
  printf(" The difference in the measured minima is %f\n",min[1]-min[0]);
  printf(" slope is %f\n",slope);
  double offset=slope*-min[0]+mid[0];
  printf(" offset is %f\n",offset);
  double edge[2]={};
  edge[0]=(cpositions[0]-offset)/slope;
  edge[1]=(cpositions[5]-offset)/slope;

  pol2->SetLineStyle(2);
  hist1->Fit("gaus","+q","",edge[0]-width,edge[0]+width);
  hist1->Fit("pol2","+q","",edge[0]-width,edge[0]+width);
  for(int i=0;i<3;i++){
    temp[i+6]=pol2->GetParameter(i);
  }
  c=temp[6];
  b=temp[7];
  a=temp[8];
  min[2]=-b/(2*a);
  printf("estimated edge is %f\n",edge[0]);
  printf("       minimum is %f\n",min[2]);

  hist1->Fit("gaus","+q","",edge[1]-width,edge[1]+width);
  hist1->Fit("pol2","+q","",edge[1]-width,edge[1]+width);
  for(int i=0;i<3;i++){
    temp[i+9]=pol2->GetParameter(i);
  }
  c=temp[9];
  b=temp[10];
  a=temp[11];
  min[3]=-b/(2*a);
  printf("estimated edge is %f\n",edge[1]);
  printf("       minimum is %f\n",min[3]);

  cFit->cd(2);
  hname=histname; 
  hname+="_copy"; 
  if ((TH1F *) gROOT->FindObject(hname)) {
    gROOT->FindObject(hname)->Delete();  
    printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
  }
  hist1->Clone(hname.Data());

  TH1F *hist2=(TH1F*) gROOT->FindObject(hname.Data());
  printf("Hisogram is %s with title %s.\n",histname,hist1->GetTitle());
  lsload(filename,slope,offset);
  //ls2load(filename,slope,offset);
  hist2->Fit("ls1","Q");
  hist2->Fit("ls1","Q");
  hist2->Fit("ls1","");
  lsgetslope();
}

void getressum(Char_t *hsum, Char_t *hcorrelation, Float_t coef_x=1., Float_t coef_y=1.)
{
  Char_t *histin1=hsum;
  Char_t *histin2=hcorrelation;
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  TH1F *hist1=(TH1F *) gROOT->FindObject(histin1);
  TH2F *hist2=(TH2F *) gROOT->FindObject(histin2);
  if(!gROOT->FindObject(histin1)) printf("Histogram \"%s\" not found.\n",histin1);
  if(!gROOT->FindObject(histin2)) printf("Histogram \"%s\" not found.\n",histin2);
  cFit->Clear();
  cFit->Divide(1,2);
  
  double var_sum = 0;
  cFit->cd(1);
  hist1->GetXaxis()->UnZoom();
  hist1->Draw("");
  printf("Statistics of %s\n",histin1);
  hist1->Fit("gaus");
  double sigma_sum = 0;
  sigma_sum = gaus->GetParameter(2);
  double mean_sum = 0;
  mean_sum = gaus->GetParameter(1);
  var_sum=TMath::Power(sigma_sum,2);
  printf(" Gassian: The width of %s is %f (FWHM = %f); variance = %f\n",histin1,sigma_sum,sigma_sum*2.35482,var_sum);
  printf(" Scaled:  The width of %s is %f (FWHM = %f); variance = %f\n",histin1,sigma_sum*0.70711,sigma_sum*2.35482*0.70711,var_sum*0.70711);
  printf(" Full:  The std dev of %s is %f (FWHM = %f); variance = %f\n",histin1,hist1->GetStdDev(),(hist1->GetStdDev())*2.35482,TMath::Power(hist1->GetStdDev(),2));
  //  printf(" The variance of the sum is %f (%f full)\n\n",var_sum,TMath::Power(hist1->GetStdDev(),2));
  hist1->GetXaxis()->SetRangeUser(mean_sum-4*sigma_sum,mean_sum+4*sigma_sum);
  hist1->Draw();

  double covar = 0;
  cFit->cd(2);
  hist2->Draw("COL"); 
  covar=hist2->GetCovariance();
  double corr=hist2->GetCorrelationFactor();
  double sig_x=hist2->GetStdDev(1);
  double sig_y=hist2->GetStdDev(2);
  double var_x=sig_x*sig_x;
  double var_y=sig_y*sig_y;
  printf("\nStatistics of %s\n",histin2);
  printf(" The standard deviations of X and Y are %f and %f (%f average)\n",sig_x,sig_y,(sig_x+sig_y)/2.); 
  printf(" The covariance of the variables is %f\n",covar);
  printf(" The correlation of the variables is  %f\n",corr);
  printf(" The calculated correlation factor is %f\n",covar/(sig_x*sig_y));
  printf(" The covariance term is %f\n",2*coef_x*coef_y*covar);
  printf(" The variances of X and Y are %f and %f (%f average)\n",var_x,var_y,(var_x+var_y)/2.);  
  printf(" The sum of the variances of X and Y is %f (direct) \n",var_x+var_y);  
  printf(" The calculated variance of the sum is  %f \n",var_x+var_y+2*coef_x*coef_y*covar);
  printf(" The calculatded standard deviation of the sum is  %f (%f FWHM)\n",TMath::Power(var_x+var_y+2*coef_x*coef_y*covar,0.5),2.35482*TMath::Power(var_x+var_y+2*coef_x*coef_y*covar,0.5));

  printf("\n The square of the sum of the variances of X and Y is %f (direct) \n",TMath::Power((var_x+var_y),0.5));  
  printf(" The square of one half the sum of the variances of X and Y is %f (direct) \n",TMath::Power((var_x+var_y)/2.,0.5));  

  printf("The square of one half the calculated sum of the variances is %f (from Gaussian)\n",TMath::Power((var_sum-2*coef_x*coef_y*covar)/2.,0.5));
  printf("The square of one half the calculated sum of the variances is %f (full)\n",TMath::Power((TMath::Power(hist1->GetStdDev(),2)-2*coef_x*coef_y*covar)/2.,0.5));
}

void grantplot(Float_t shift=-0.85, Bool_t set_log=1, Int_t plots=1, Int_t test=1)
{
  //set canvas appearence
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  cFit->SetLogz(set_log);
 
  //set hitogram appearence
  Float_t nbins=hhitc0->GetNbinsX();
  Float_t set_bins=375;
  Float_t rebin=nbins/set_bins;
  hhitc0->Rebin2D(rebin,rebin);
  TH2F * hOutput1=(TH2F *) gROOT->FindObject("hhitc0");
  hOutput1->SetYTitle("Vertical Position (mm)");
  //hOutput1->GetYaxis()->CenterTitle(1);
  hOutput1->SetXTitle("Horizontal Position (mm)");
  // hOutput1->GetXaxis()->CenterTitle(1);
  hOutput1->SetTitle("Detector 1 Y vs. X, Calibrated");
 
  //draw overlays
  shadowzc("hhitc0",z_A1); //set line spans
  dr("hhitc0"); //redraw histogram to clear lines
  //if(test==1)
  maskz(z_A1,1);
  y_cal=27+shift;  //change position lines
  if(set_log)
    shieldz(z_A1,1,4);
  else
    shieldz(z_A1,1,2);
  prop(1);
  if(test==1) {
    cFit->SaveAs("figures/run_430_hhitc0_overlay_c_log2.pdf");
    cFit->SaveAs("figures/run_430_hhitc0_overlay_c_log2.eps");
  }
  else {
    cFit->SaveAs("figures/run_606_hhitc0_overlay_c_log2.pdf");
    cFit->SaveAs("figures/run_606_hhitc0_overlay_c_log2.jpg");
  }
    
  
  if(plots>1) {
    //plot with scale
    cFit->SetRightMargin(0.15);
    hOutput1->SetZTitle("Counts per 0.032 mm^{2}");
    hOutput1->SetZTitle("Counts per 0.128 mm^{2}");
    hOutput1->GetZaxis()->SetTitleOffset(1.25);
    //hOutput1->GetZaxis()->CenterTitle(1);
    hOutput1->Draw("col");
    hOutput1->Draw("colz");
    y_cal=27;  //change position lines
    maskz(z_A1,1);
    y_cal=27+shift;  //change position lines
    if(set_log)
      shieldz(z_A1,1,4);
    else 
      shieldz(z_A1,1,2);
    cFit->SaveAs("figures/run_430_hhitc0_overlay_c_log2z.pdf");
    cFit->SaveAs("figures/run_430_hhitc0_overlay_c_log2z.eps");
  }
}

void resplot()
{
  setbeam();
  setangles(0,4.48);
  compX(0.43,0.608,1e5,50,1);
  prop(16,9);
  cFit->SaveAs("figures/run_430_compX2.pdf");
}

void printruns(Int_t run_start=480, Int_t run_stop=480, Bool_t plots=kFALSE)
{
  Float_t counts=0;
  FILE * outfile;
  TPaveText *pt;
  TText *text;  
  TString rname="";
  outfile=fopen("runs.lst","w");
  fprintf(outfile,"generated with the command printruns(%d,%d,%s)\n",run_start,run_stop,plots ? "kTRUE" : "kFALSE");
  fprintf(outfile,"Run, Anode Coincidences, Time width, Time peak, X1 width, X1, Y1 width, Y1, X2 width, X2, Y2 width, Y2\n");
  fprintf(outfile,"Run, Anode Coincidences, Time B width; Time B peak; Time M width; Time M peak; Time T width; Time T peak, X1 width, X1, Y1 width, Y1, X2 width, X2, Y2 width, Y2\n");
  fclose(outfile);
  
  FILE * outfile2;
  outfile2=fopen("nevents.lst","w");
  fprintf(outfile2,"generated with the command printruns(%d,%d,%s)\n",run_start,run_stop,plots ? "kTRUE" : "kFALSE");
  if(run_start==397)
    fprintf(outfile2,"const Int_t nevents[%d] = {",run_stop-run_start+1);
  fclose(outfile2);
  
  for(Int_t i=run_start;i<(run_stop+1);i++) {
    hname="output/emma_ana_00";
    hname+=i;
    hname+=".root";
    printf("Input file is %s\n",hname.Data());
    outfile=fopen("runs.lst","a");
    outfile2=fopen("nevents.lst","a");
    TFile *_file0 = TFile::Open(hname.Data());
    if(_file0) {
      counts=  hcounts15->GetBinContent(1);//added for 2015 tests
      counts+=  hcounts15->GetBinContent(2);
      counts+=  hcounts15->GetBinContent(3);
      printf(" Total anode coincidences is %f\n",counts);
      if(i<481) {
	gfitcp("hdiffz1",-1,10,2,0.05,"q");
	fprintf(outfile,"%d, %.1f, %f, %f, ",i,counts,gaus->GetParameter(2),gratio);
	fprintf(outfile2,"%.0f",counts);
      }
      else {
	fprintf(outfile,"%d, %.1f, ",i,counts);
	fprintf(outfile2,"%.0f",counts);	
	gfitcp("htimez0",-1,10,2,0.05,"q");
	fprintf(outfile,"%f; %f; ",gaus->GetParameter(2),gratio);
	gfitcp("htimez4",-1,10,2,0.05,"q");
	fprintf(outfile,"%f; %f; ",gaus->GetParameter(2),gratio);
	gfitcp("htimez8",-1,10,2,0.05,"q");
	fprintf(outfile,"%f; %f, ",gaus->GetParameter(2),gratio);
      }
      gfitcp("hsumz3",-1,30,2,0.05,"q");
      fprintf(outfile,"%f, %f, ",gaus->GetParameter(2),gratio);
      gfitcp("hsumz4",-1,30,2,0.05,"q");
      fprintf(outfile,"%f, %f, ",gaus->GetParameter(2),gratio);
      gfitcp("hsumz5",-1,30,2,0.05,"q");
      fprintf(outfile,"%f, %f, ",gaus->GetParameter(2),gratio);
      gfitcp("hsumz6",-1,30,2,0.05,"q");
      fprintf(outfile,"%f, %f\n",gaus->GetParameter(2),gratio);
      if(plots) {
	rname="";
	plotall("hxx");
	prop(4,3);
	pt = new TPaveText(0.937,0.008,0.988,0.072,"NDC");
	pt->SetFillColor(0); 
	pt->SetShadowColor(0);
	rname+=i;
	printf("%s",rname.Data());
	text = pt->AddText(rname.Data());
	pt->Draw();
	hname="figures/run_";
	hname+=i;
	hname+="_hxx.pdf";	
	cFit->SaveAs(hname.Data());
	plotall("hdiffz");
	pt->Draw();
	hname="figures/run_";
	hname+=i;
	hname+="_hdiffz.pdf";	
	cFit->SaveAs(hname.Data());
	plotall("htsum");
	pt->Draw();
	hname="figures/run_";
	hname+=i;
	hname+="_htsum.pdf";	
	cFit->SaveAs(hname.Data());
	plotall("hang");
	pt->Draw();
	hname="figures/run_";
	hname+=i;
	hname+="_hang.pdf";	
	cFit->SaveAs(hname.Data());
      }
    }
    else {
      printf("File %s not found\n",hname.Data());
      fprintf(outfile,"%d, %.1f, %f, 0, 0, 0, 0\n",i,0,0);
      fprintf(outfile2,"0");
    }
    if(i==run_stop)
      fprintf(outfile2,"};");
    else
      fprintf(outfile2,", ");
    fclose(outfile);
    fclose(outfile2);
  }
}

//-----------------------------------------------------------------------------------------------
//macros for the EMMA IC-------------------------------------------------------------------------

void Gaus(double mean1, double sigma1, double mean2, double sigma2)
{//developed by Lily
 //creates a 2D histogram of a Gaussian distribution with given means and widths
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gRandom = new TRandom3();

  TH2F * hist = new TH2F("data", "hist", 100, 0.0, 100.0,100,0.0,100.0);

  for (int i = 0; i < 10000; ++i)
    hist->Fill(gRandom->Gaus(mean1, sigma1),gRandom->Gaus(mean2, sigma2));

  TCanvas * c1 = new TCanvas ("c1", "fitted data", 5, 5, 800, 600);

  // hist->Fit("gaus");
  hist->Draw();
  hist->SaveAs("fit.eps");
}

void mkgaus2d(Char_t *histin,Int_t points=10000,double mean1, double sigma1, double mean2, double sigma2, Int_t color=2, Float_t blur=0.0,Int_t bins=1000)
{//creates a 2D histogram of a Gaussian distribution with given means and widths

  gRandom = new TRandom3();

  hname=histin; 
  if ((TH1F *) gROOT->FindObject(hname)) {
    printf("Histogram \"%s\" exists. ",hname.Data());
  }
  else {
  
  TH2F * hist = new TH2F(hname.Data(),hname.Data(),bins,mean1-(5*sigma1),mean1+(5*sigma1),bins,mean2-(5*sigma2),mean2+(5*sigma2));
  //TH2F * hist = new TH2F(hname.Data(), hname.Data(), bins, 0.,100.,bins,0.,100.);
  }
   hInput=(TH2F*)gROOT->FindObject(hname.Data());
   hInput->SetMarkerColor(color);

   printf("Distribution is %f wide (x) and %f wide (y)",sigma1/mean1,sigma2/mean2);
   if(blur==0) 
     printf("No blurring applied.\n");
   else
     printf("%f%% blurring applied.\n",blur*100);
   blur/=2.354820045;
   for (int i = 0; i < points; ++i) {
     if(blur==0) 
       hInput->Fill(gRandom->Gaus(mean1,sigma1),gRandom->Gaus(mean2,sigma2));
     else { 
       if(sigma1==0) {
	 if(sigma2==0) //both 0
	   hInput->Fill(gRandom->Gaus(mean1,mean1*blur),
			gRandom->Gaus(mean2,mean2*blur));
	 else //only sigma1=0
	   hInput->Fill(gRandom->Gaus(mean1,mean1*blur),
			gRandom->Gaus(mean2, sigma2)+gRandom->Gaus(0,mean2*blur));
       }
       else {//sigma1 non-zero
	 if(sigma2==0)//only sigma2=0
	   hInput->Fill(gRandom->Gaus(mean1,sigma1)+gRandom->Gaus(0,mean1*blur),
			gRandom->Gaus(mean2, mean2*blur));
	 else //neither
	   hInput->Fill(gRandom->Gaus(mean1,sigma1)+gRandom->Gaus(0,mean1*blur),
			gRandom->Gaus(mean2, sigma2)+gRandom->Gaus(0,mean2*blur));
       }
     }
   }

   if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();//added
   //  // hist->Fit("gaus");
   //if(bins<1000)
    hInput->Draw();
  //  // hist->SaveAs("fit.eps");
}

void smearedGaus(double mean1, double sigma1, double mean2, double sigma2, double smear) 
{//developed by Lily 
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gRandom = new TRandom3();

  TH2F * hist = new TH2F("data", "hist", 100, 0.0, 100.0,100,0.0,100.0);

  for (int i = 0; i < 10000; ++i)
    hist->Fill(gRandom->Gaus(mean1, sigma1)+gRandom->Gaus(0,2.3548*smear),gRandom->Gaus(mean2, sigma2)+gRandom->Gaus(0,2.3548*smear));

  TCanvas * c1 = new TCanvas ("c1", "fitted data", 5, 5, 800, 600);

  // hist->Fit("gaus");
  hist->Draw();
  hist->SaveAs("fit.eps");
}

void plotgaus(Int_t bins=100, Int_t counts=10000, Bool_t set_strag=kTRUE, Float_t det_res=0.036, Float_t x_low=50, Float_t x_high=72, Float_t y_low=35, Float_t y_high=48)
{
  mkhist2d("hist1", bins,x_low,x_high,bins,y_low,y_high,kFALSE);
  mkhist2d("hist2", bins,x_low,x_high,bins,y_low,y_high,kFALSE);
  mkhist2d("hist3", bins,x_low,x_high,bins,y_low,y_high,kFALSE);
  printf("Axis proportion scale is %f\n",(x_high-x_low)/(y_high-y_low));
  
  Int_t mark_style=20;
  hist1->SetMarkerStyle(mark_style);
  hist2->SetMarkerStyle(mark_style);
  hist3->SetMarkerStyle(mark_style);
  Float_t mark_size=0.4;
  hist1->SetMarkerSize(mark_size);  
  hist2->SetMarkerSize(mark_size);
  hist3->SetMarkerSize(mark_size);

   Int_t steps=counts;

   Int_t strag=0;
   if(set_strag)
     strag=1;
   //   Int_t blur=0;
   //if(det_res)
   //  blur=0.036;

   //no straggling
   //mkgaus2d("hist1",steps/10, 64.120, 0, 39.585, 0,1,det_res,bins); //Rb
   //mkgaus2d("hist2",steps,58.211, 0, 41.763, 0,2,det_res,bins); //Sr
   //mkgaus2d("hist3",steps/10,57.472, 0, 41.526, 0,4,det_res,bins); //Y

   //calculated straggling (supposed to be 3.6%)
   //mkgaus2d("hist1",steps/10, 64.120, 1.359, 39.585, 0.839,1,det_res,bins); //Rb
   //mkgaus2d("hist2",steps,58.211, 1.234, 41.763, 0.885,2,det_res,bins); //Sr
   //mkgaus2d("hist3",steps/10,57.472, 1.218, 41.526, 0.880 ,4,det_res,bins); //Y
   
   //real straggling, parameterized
   //mkgaus2d("hist1",steps/10, 64.120, strag*1.186, 39.585, strag*0.903,1,det_res,bins); //Rb
   //mkgaus2d("hist2",steps,58.211, strag*1.020, 41.763, strag*1.154,2,det_res,bins); //Sr
   //mkgaus2d("hist3",steps/10,57.472, strag*0.616, 41.526, strag*0.670,4,det_res,bins); //Y

   //updated energies, symmetric straggling
   //mkgaus2d("hist1",steps/10, 64.117, strag*1.065, 39.562, strag*1.065,1,det_res,bins); //Rb
   //mkgaus2d("hist2",steps,58.032, strag*1.300, 41.846, strag*1.300,2,det_res,bins); //Sr
   //mkgaus2d("hist3",steps/10,57.308, strag*0.588, 41.473, strag*0.588,4,det_res,bins); //Y
   
   //Butane only, from SRIM
   //mkgaus2d("hist1",steps/10,69.919, strag*0.896, 33.794, strag*0.896 ,1,det_res,bins); //Rb
   //mkgaus2d("hist2",steps,64.186, strag*0.997, 35.670, strag*0.997 ,2,det_res,bins); //Sr
   //mkgaus2d("hist3",steps/10,63.246, strag*0.597, 35.502, strag*0.597 ,4,det_res,bins); //Y

   //updated energies, symmetric straggling, scaled energy loses
   mkgaus2d("hist1",steps/10, 56.624, strag*1.269, 40.830, strag*1.269
	    ,1,det_res,bins); //Rb
   mkgaus2d("hist2",steps,58.032, strag*1.300, 41.846, strag*1.300
	    ,2,det_res,bins); //Sr
   mkgaus2d("hist3",steps/10,59.413, strag*1.331, 42.841, strag*1.331
	    ,4,det_res,bins); //Y
   
  hist1->SetTitle("EMMA IC");
  hist1->GetXaxis()->SetTitle("E_{res} (MeV)");
  hist1->GetYaxis()->SetTitle("#\DeltaE (MeV)");
  hist1->Draw();
  hist2->Draw("same");
  hist3->Draw("same");
  
 leg = new TLegend(0.1,0.75,0.2,0.9);
 leg->AddEntry(hist1,"^{84}Rb","p");   
 leg->AddEntry(hist2,"^{84}Sr","p");
 leg->AddEntry(hist3,"^{84}Y","p");
 leg->Draw();
}

//-----------------------------------------------------------------------------------------------
//macros for EMMA tree analysis -----------------------------------------------------------------
void treecomp(Int_t detno=0)
{
  TString names[7]={ "ab","am","at","xf","xn","yn","yf"};//please note that here the order of signals matches the readout and now the tree order
  //  gEnv->SetValue("Hist.Binning.1D.x",100);
  hname="hdata";
  hname+=detno;
  printf(" Input histogram name is %s\n",hname.Data());
  dr(hname.Data());
  hProj=(TH1F *) gROOT->FindObject(hname.Data());
  //Int_t set_bins=hProj->GetXaxis()->GetNbins();
  //gEnv->SetValue("Hist.Binning.1D.x",set_bins);
  hname="hdata0_copy";
  if(!(TH1F*)gROOT->FindObject(hname.Data())){
    hProj->Clone(hname.Data());
  }
  hFit=(TH1F *) gROOT->FindObject(hname.Data());
  hFit->Reset();
  cFit->SetLogy();
  hFit->SetLineColor(2);
  hFit->SetLineStyle(2);
  Int_t signo=0;
  Bool_t pgac1=kTRUE;
  if(detno<6) {//anode data
    if(detno<3) //detector 1
      signo=detno;
    else {//detector 2
      signo=detno-3;
      pgac1=kFALSE;
  }
  }
  else {//cathode data 
    if(detno<10)//detector 1
      signo=detno-3;
    else {
      signo=detno-7;
      pgac1=kFALSE;
    }
  }
  
  printf("signal number is %d on detector %d\n",signo,pgac1+1);
  //TString tname;
  if(pgac1)
    hname="PGAC1.";
  else
    hname="PGAC2.";  

  hname+=names[signo];
  // hname+=" \x003E\x003E hdata0_copy";
  hname+=">";
  hname+=">";
  hname+="hdata0_copy";
  printf("Target branch is %s\n",hname.Data());
  TString cname=names[signo];
  cname+=">0";
  TCut c1 = cname.Data();
  printf("Cut is %s\n",cname.Data());

  t1->Draw(hname.Data(),c1,"same");
}

void anode (Int_t detno=0)
{
  TString names[7]={ "ab","am","at","xf","xn","yn","yf"};//please note that here the order of signals matches the readout and now the tree order
  Int_t signo=0;
  Bool_t pgac1=kTRUE;
  if(detno<6) {//anode data
    if(detno<3) //detector 1
      signo=detno;
    else {//detector 2
      signo=detno-3;
      pgac1=kFALSE;
    }
  }
  else {//cathode data 
    if(detno<10)//detector 1
      signo=detno-3;
    else {
      signo=detno-7;
      pgac1=kFALSE;
    }
  }
  
  printf("signal number is %d on detector %d\n",signo,pgac1+1);
  //TString tname;
  if(pgac1)
    hname="PGAC1.";
  else
    hname="PGAC2.";  

  hname+=names[signo];
  hname+=">";
  hname+=">";
  hname+="hdata0_copy";
  printf("Target branch is %s\n",hname.Data());
  TString cname=names[signo];
  cname+=">0";
  TCut c1 = cname.Data();
  printf("Cut is %s\n",cname.Data());

  t1->Draw(hname.Data(),c1,"same");
}

//-----------------------------------------------------------------------------------------------
//macros for SE-SPS GFM calculations (September 2016)--------------------------------------------
void fpplot()
{//from June 2016
  mkhist2d("h",3200,-30,100,100,-1,1);
  fill2("18O_He_1e5.out.txt","h",1);
  h->Clone("new");
  fill2("224Pa_He_xy.out.csv","new",1);
  new->SetTitle("Focal Plane Position");
  new->SetXTitle("x-position (cm)");
  new->SetYTitle("y-position (cm)");
  h->Draw("same");
  leg = new TLegend(0.1,0.75,.2,.9);
  leg->AddEntry(new,"^{224}Pa","p");
  leg->AddEntry(h,"^{18}O","p");
  leg->Draw();  
}

void Gaus1(double mean1, double sigma1,int counts=10000)
{//extension to Gaus(), creates a 1D histogram of a Gaussian distribution with given mean and width
  gRandom = new TRandom3();
  gRandom->SetSeed(0);
  
  hname="data"; 
  if ((TH1F *) gROOT->FindObject(hname)) {
    gROOT->FindObject(hname)->Delete();  
    printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
  }

  TH1F * hist = new TH1F(hname, "hist", 100, 0.0, 100.0);

  for (int i = 0; i < counts; ++i)
    hist->Fill(gRandom->Gaus(mean1, sigma1));

  if(!((TCanvas *) gROOT->FindObject("c1")))
    TCanvas * c1 = new TCanvas ("c1", "fitted data", 5, 5, 800, 600);
  
  hist->Draw();
  gfitc(hname,mean1,4*sigma1);
}

void Gaus1b(double mean1, double sigma1, int bins, int width=5,int counts=1000)
{//creates a Gaussian distribution with a certain mean, number of bins
  gRandom = new TRandom3();
  gRandom->SetSeed(0);
  bins=(2*bins)+1;
  printf("Number of bins is 2N+1=%d\n",bins);
  Float_t xmin=mean1-width*sigma1;
  Float_t xmax=mean1+width*sigma1;
  Float_t bwid=(xmax-xmin)/bins;
  printf("Calculated bin width is %f\n",bwid);

  hname="data"; 
  if ((TH1F *) gROOT->FindObject(hname)) {
    gROOT->FindObject(hname)->Delete();  
    printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
  }
  TH1F * hist = new TH1F(hname, "Rand Gaus w/ spec N bins", bins, xmin, xmax);
  printf("    Actual bin width is %f\n",hist->GetBinWidth(1));
  
  for (int i = 0; i < counts; ++i)
    hist->Fill(gRandom->Gaus(mean1, sigma1));
  if(!((TCanvas *) gROOT->FindObject("c1"))) 
    TCanvas * c1 = new TCanvas ("c1", "fitted data", 5, 5, 800, 600);
  
  hist->Draw();
  gfitc(hname,mean1,width*sigma1);

  FILE * outfile;
  outfile=fopen("temp.lst","w");
  
  Float_t center=0;
  Float_t cont=0;
  Float_t inte=0;

    for (int i = 1; i < bins+1; ++i) {
      center=hist->GetBinCenter(i);
      cont=hist->GetBinContent(i);
      inte+=cont;
      if(cont>0&&center>0) {
	printf("Bin %9.4f has %6.0f counts\n",center,cont);
	fprintf(outfile,"%10g %d\n",center,cont);
      }
    }
    printf("  Total counts is %6.0f counts\n",inte);
    fclose(outfile);
}

void Gaus1bw(double mean1, double sigma1, double bwid, int width=5,int counts=1000)
{//creates a Gaussian distribution with a certain mean, sigma, and bin width
  gRandom = new TRandom3();
  gRandom->SetSeed(0);
  Float_t xmin=mean1-width*sigma1;
  Float_t xmax=mean1+width*sigma1;
  Float_t bins=(Int_t)((xmax-xmin)/bwid);
  xmax=xmin+bins*bwid;
  printf("Bin width is %f\n",bwid);
  printf("%d bins to cover +/-%.1f sigma\n",bins,width);

  hname="data"; 
  if ((TH1F *) gROOT->FindObject(hname)) {
    gROOT->FindObject(hname)->Delete();  
    printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
  }
  
  TH1F * hist = new TH1F(hname, "Rand Gaus w/ fixed bin width", bins, xmin, xmax);

  printf("bin width is %f\n",hist->GetBinWidth(1));
  
  for (int i = 0; i < counts; ++i)
    hist->Fill(gRandom->Gaus(mean1, sigma1));
  if(!((TCanvas *) gROOT->FindObject("c1"))) 
    TCanvas * c1 = new TCanvas ("c1", "fitted data", 5, 5, 800, 600);
  
  hist->Draw();
  gfitc(hname,mean1,width*sigma1);

  FILE * outfile;
  outfile=fopen("temp.lst","w");
  Float_t center=0;
  Float_t cont=0;
  Float_t inte=0;
  Float_t intec=0;
  for (int i = 1; i < bins+1; ++i) {
    center=hist->GetBinCenter(i);
    cont=hist->GetBinContent(i);
    inte+=cont;
    if(cont>(counts/100)&&center>0) {
      printf("Bin %9.4f has %6.0f counts\n",center,cont);
      fprintf(outfile,"%10g %d\n",center,cont);
      intec+=cont;
    }
  }
  printf("  Total counts is %6.0f counts\n",inte);
  printf("  Gated counts is %6.0f counts\n",intec);
   fclose(outfile);
}

void uni1bw(double mean1, double sigma1, double bwid, int width=5,int counts=1000)
{//creates a uniform distribution with a certain mean, sigma, and bin width
  gRandom = new TRandom3();
  gRandom->SetSeed(0);
  Float_t xmin=mean1-width*sigma1;
  Float_t xmax=mean1+width*sigma1;
  Float_t bins=(Int_t)((xmax-xmin)/bwid);
  xmax=xmin+bins*bwid;
  printf("Bin width is %f\n",bwid);
  printf("%d bins to cover +/-%.1f * %.1f\n",bins,width,sigma1);

  hname="data"; 
  if ((TH1F *) gROOT->FindObject(hname)) {
    gROOT->FindObject(hname)->Delete();  
    printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
  }
  
  TH1F * hist = new TH1F(hname, "Uniform random w/ fixed bin width", bins, xmin, xmax);

  printf("bin width is %f\n",hist->GetBinWidth(1));

  for (int i = 0; i < counts; ++i)
    hist->Fill(gRandom->Uniform(xmin,xmax));

  if(!((TCanvas *) gROOT->FindObject("c1")))
    TCanvas * c1 = new TCanvas ("c1", "fitted data", 5, 5, 800, 600);
  
  hist->Draw();
  hist->GetYaxis()->SetRangeUser(0,(1.2*(Float_t)counts/bins));
   
  FILE * outfile;
  outfile=fopen("temp.lst","w");
  Float_t center=0;
  Float_t cont=0;
  Float_t inte=0;
  Float_t intec=0;
  for (int i = 1; i < bins+1; ++i) {
    center=hist->GetBinCenter(i);
    cont=hist->GetBinContent(i);
    inte+=cont;
    if(cont>(0)&&center>0) {
      //      printf("Bin %9.4f has %6.0f counts\n",center,cont);
      fprintf(outfile,"%10g %d\n",center,cont);
      intec+=cont;
    }
  }
  printf("  Total counts is %6.0f counts\n",inte);
  printf("  Gated counts is %6.0f counts\n",intec);
   fclose(outfile);
}

Int_t  bins=512*2;
  Float_t xmin=-85;
  Float_t xmax=130;
  Float_t ymin=-2;
  Float_t ymax=2;
Float_t tmin=400;//0;
Float_t tmax=750;//1400;
  Float_t emin=0;
  Float_t emax=8;

  TH2F *h1;
  TH2F *h2;
  TH2F *h3;
  TH2F *h4;
TTree *t = new TTree();

Float_t a,b,c;

mca2root(TString fname="output.mca")
{
  TFile *f = new TFile("output.root","recreate");
  
  //t->ReadFile(fname,"x/F:y/F:e/F:q/F:t/F");
  t->ReadFile(fname,"x/F:y/F:e/F:q/F:t/F:theta/F:phi/F");
  t->Write("");
  if(!((TCanvas *) gROOT->FindObject("cFit"))) {
    TCanvas * cFit=new TCanvas("cFit","cFit",0,0,675,615);
    if(!(cFit->GetShowEventStatus()))cFit->ToggleEventStatus();
    if(!(cFit->GetShowToolBar()))cFit->ToggleToolBar();
  }
  
  for(int i=1;i<5;i++){
    TString hname="h";
    hname+=i;
    if ((TH2F *) gROOT->FindObject(hname.Data())) {
      gROOT->FindObject(hname.Data())->Delete();  
      printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
    }
  }
  
  h1 = new TH2F("h1","",bins,xmin,xmax,bins,ymin,ymax);
  h2 = new TH2F("h2","",bins,xmin,xmax,bins,tmin,tmax);
  h3 = new TH2F("h3","",bins,xmin,xmax,bins,emin,emax);
  h4 = new TH2F("h4","",bins,emin,emax,bins,tmin,tmax);

  h1->SetXTitle("x-position (cm)");
  h2->SetXTitle("x-position (cm)");
  h3->SetXTitle("x-position (cm)");
  h4->SetXTitle("Energy (MeV)");
  h1->SetYTitle("y-position (cm)");
  h2->SetYTitle("time-of-flight (ns)");
  h3->SetYTitle("Energy (MeV)");
  h4->SetYTitle("time-of-flight (ns)");
  
  cFit->Divide(2,2);
  cFit->cd(1);
  t->Draw("y:x>>h1","","col");
  cFit->cd(2);
  t->Draw("t:x>>h2","","col");
  cFit->cd(3);
  t->Draw("e:x>>h3","","col");
  cFit->cd(4);
  t->Draw("t:e>>h4","","col");
  
  //calculate time correction
  //t->Draw("t:theta>>htemp","","col");
  //t->Draw("t:x>>htemp","","col");
  //htemp->ProfileX();
  //htemp_pfx->Fit("pol2","m");
   
  //h4->Draw("col");
  //h4->Fit("pol2","m");
  //copy2("h4",31);
  //h4_copy->Fit("pol2","m");
  pfx("h4");
  xprof->Fit("pol2","m");
      
  //Float_t a=pol2->GetParameter(2),b=pol2->GetParameter(1),c=pol2->GetParameter(0);
  a=pol2->GetParameter(2);b=pol2->GetParameter(1);c=pol2->GetParameter(0);
  //  a=pol3->GetParameter(2);b=pol3->GetParameter(1);c=pol3->GetParameter(0);
  //  Float_t cc=0;
  //  cc=pol3->GetParameter(3);
  Float_t mean=h4->GetMean(2);

  plotall("h");

  h1->Clone("hc1");
  h1->Clone("hcc1");
  h1->Clone("hccc1");

  h3->Clone("hc3");
  h3->Clone("hcc3");
  h3->Clone("hccc3");
  
  cFit->cd(2);
  h2->Clone("hc2");
  h2->Clone("hcc2");
  hcc2->SetMarkerColor(2);
  h2->Clone("hccc2");
  hccc2->SetMarkerColor(3);
  //t->Draw("t-pol2->GetParameter(2)*theta*theta-pol2->GetParameter(1)*theta:x>>hc2","","col");
  //t->Draw(TString::Format("t-%f*theta*theta-%f*theta:x>>hc2",a,b),"","col");
  //t->Draw("t-0.00442252*theta*theta+9.13593*theta:e>>hc2","","col");
  //t->Draw(TString::Format("t-((%f)*theta*theta+(%f)*theta+(%f))+700:x>>hc2",a,b,c),"","col");
  //t->Draw(TString::Format("t-%f-%f*TMath::Power(theta,2)-%f*theta:x>>hc2",c-mean,a,b),"","col");
  t->Draw(TString::Format("t-%f-%f*TMath::Power(e,2)-%f*e:x>>hc2",c-mean,a,b),"","same");
  //t->Draw(TString::Format("t-%f-%f*TMath::Power(x,2)-%f*x:x>>hc2",c-mean,a,b),"","same");

  TCanvas::MakeDefCanvas();
  c1->cd();
  pfx("hc2");
  xprof->Fit("pol1","m");
  
  cFit->cd(2);
  h2->Draw("col");
  hc2->Draw("same");
  Float_t d=pol1->GetParameter(1), e=pol1->GetParameter(0);
  t->Draw(TString::Format("(t-%f-%f*TMath::Power(e,2)-%f*e)-%f*x-%f:x>>hcc2",c-mean,a,b,d,e-mean),"","same");

  c1->cd();
  t->Draw(TString::Format("(t-%f-%f*TMath::Power(e,2)-%f*e)-%f*x-%f:theta",c-mean,a,b,d,e-mean),"","col");
  htemp->Draw("col");
  pfx("htemp");
  xprof->Fit("pol1","m");
  
  Float_t ff=pol1->GetParameter(1), g=pol1->GetParameter(0);
  cFit->cd(2);
  t->Draw(TString::Format("((t-%f-%f*TMath::Power(e,2)-%f*e)-%f*x-%f)-%f*theta-%f:x>>hccc2",c-mean,a,b,d,e-mean,ff,g-mean),"","same");
  
  cFit->cd(4);
  h4->Clone("hc4");
  h4->Clone("hcc4");
  hcc4->SetMarkerColor(2);
  h4->Clone("hccc4");
  hccc4->SetMarkerColor(3);
  //t->Draw("t-pol2->GetParameter(2)*theta*theta-pol2->GetParameter(1)*theta:e>>hc4","","col");//works
  //t->Draw("t-0.004423*theta*theta+9.13593*theta:e>>h4","","col");//doesn't work!?
  //t->Draw(TString::Format("t-%f*TMath::Power(e,2)-%f*e:x>>hc4",a,b),"","col");
  //t->Draw(TString::Format("t-%f*TMath::Power(theta,2)-%f*theta:e>>hc4",a,b),"","col");
  t->Draw(TString::Format("t-%f-%f*TMath::Power(e,2)-%f*e:e>>hc4",c-mean,a,b),"","same");
  //t->Draw(TString::Format("t-%f-%f*TMath::Power(x,2)-%f*x:e>>hc4",c-mean,a,b),"","same");
  t->Draw(TString::Format("(t-%f-%f*TMath::Power(e,2)-%f*e)-%f*x-%f:e>>hcc4",c-mean,a,b,d,e-mean),"","same");
  t->Draw(TString::Format("((t-%f-%f*TMath::Power(e,2)-%f*e)-%f*x-%f)-%f*theta-%f:e>>hccc4",c-mean,a,b,d,e-mean,ff,g-mean),"","same");

  c1->cd();
  //pjx("h2");gfitcp("xproj");
  pjy("h2");gfitcp("yproj");
  pjy("hc2");gfitcp("yproj");
  pjy("hcc2");gfitcp("yproj");
  pjy("hccc2");gfitcp("yproj");  
  
  //h4->ProfileX();
  //h4_pfx->Fit("pol2","m");
  //a=pol2->GetParameter(2);b=pol2->GetParameter(1);c=pol2->GetParameter(0);
  //t->Draw(TString::Format("t-%f-%f*TMath::Power(e,2)-%f*e:e>>hc4",c-mean,a,b),"","same");
  //  pjy("hc4");gfitcp("yproj"); 
  
  
  //h2->Clone("hcc2");
  //hc2->ProfileX();
  //hc2_pfx->Fit("pol1","m");
  
  //d=pol1->GetParameter(2);
  //e=pol1->GetParameter(1);
  //f=pol1->GetParameter(0);

  /*
  
  float tt;
  float time,theta;
  TBranch *btt = t->Branch("tt",&tt,"tt/F");
  t->SetBranchAddress("t",&time);
  t->SetBranchAddress("theta",&theta);
  Long64_t nentries = t->GetEntries();
  for (Long64_t j=0;j<nentries;j++) {
    t->GetEntry(j);
        tt = time-pol2->GetParameter(2)*theta*theta-pol2->GetParameter(1)*theta;
    //tt = time-a*theta*theta-b*theta;
    btt->Fill();
  }
   //   t->Print();
  t->Write();
  */
}

Int_t mark_style=1;//20;
Float_t mark_size=0.4;
Int_t mark_col=2;

TTree *t2 = new TTree();

void mca2root2(TString fname="output2.mca")
{//assumes mca2root.C has been run
 
for(int i=1;i<5;i++){
    TString hname="h";
    hname+=i;
    TH2F * hInput=(TH2F *) gROOT->FindObject(hname.Data());
    hInput->SetMarkerStyle(mark_style);
    hInput->SetMarkerSize(mark_size);
    hInput->SetMarkerColor(mark_col);
    cFit->cd(i);
    hInput->Draw();
  }

//t2->ReadFile(fname,"x/F:y/F:e/F:q/F:t/F");
 t2->ReadFile(fname,"x/F:y/F:e/F:q/F:t/F:theta/F:phi/F");
 t2->Write("");

 for(int i=5;i<9;i++){
   TString hname="h";
   hname+=i;
   if ((TH2F *) gROOT->FindObject(hname.Data())) {
     gROOT->FindObject(hname.Data())->Delete();  
     printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
   }
 }

  TH2F *h5;
  TH2F *h6;
  TH2F *h7;
  TH2F *h8;
  
  h5 = new TH2F("h5","",bins,xmin,xmax,bins,ymin,ymax);
  h6 = new TH2F("h6","",bins,xmin,xmax,bins,tmin,tmax);
  h7 = new TH2F("h7","",bins,xmin,xmax,bins,emin,emax);
  h8 = new TH2F("h8","",bins,emin,emax,bins,tmin,tmax);

  h5->SetXTitle("x-position (cm)");
  h6->SetXTitle("x-position (cm)");
  h7->SetXTitle("x-position (cm)");
  h8->SetXTitle("Energy (MeV)");
  h5->SetYTitle("y-position (cm)");
  h6->SetYTitle("time-of-flight (ns)");
  h7->SetYTitle("Energy (MeV)");
  h8->SetYTitle("time-of-flight (ns)");

  for(int i=5;i<9;i++){
    TString hname="h";
    hname+=i;
    TH2F * hInput=(TH2F *) gROOT->FindObject(hname.Data());
    hInput->SetMarkerStyle(mark_style);
    hInput->SetMarkerSize(mark_size);
    hInput->SetMarkerColor(4);
  }
  
  //cFit->Divide(2,2);
  cFit->cd(1);
  t2->Draw("y:x>>h5","","same");
  cFit->cd(2);
  //t2->Draw("t:x>>h6","","same");
  cFit->cd(3);
  t2->Draw("e:x>>h7","","same");
  cFit->cd(4);
  //t2->Draw("t:e>>h8","","same");
  h1->Write("");
  h2->Write("");
  h3->Write("");
  h4->Write("");
  h5->Write("");
  h6->Write("");
  h7->Write("");
  h8->Write("");

  cFit->cd(2);
  t2->Draw("t-pol2->GetParameter(2)*theta*theta-pol2->GetParameter(1)*theta:x>>h6","","same");
  h2->Draw("same");
    
  cFit->cd(4);
  t2->Draw("t-pol2->GetParameter(2)*theta*theta-pol2->GetParameter(1)*theta:e>>h8","","same");
}

  TTree *t3 = new TTree();
void mca2root3(TString fname="output3.mca")
{//assumes mca2root.C has been run
  

  //  t3->ReadFile(fname,"x/F:y/F:e/F:q/F:t/F");
  t3->ReadFile(fname,"x/F:y/F:e/F:q/F:t/F:theta/F:phi/F");
  t3->Write("");
  
  for(int i=9;i<13;i++){
    TString hname="h";
    hname+=i;
    if ((TH2F *) gROOT->FindObject(hname.Data())) {
      gROOT->FindObject(hname.Data())->Delete();  
      printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
    }
  }
  
  TH2F *h9;
  TH2F *h10;
  TH2F *h11;
  TH2F *h12;
  
  h9 = new TH2F("h9","",bins,xmin,xmax,bins,ymin,ymax);
  h10 = new TH2F("h10","",bins,xmin,xmax,bins,tmin,tmax);
  h11 = new TH2F("h11","",bins,xmin,xmax,bins,emin,emax);
  h12 = new TH2F("h12","",bins,emin,emax,bins,tmin,tmax);

  h9->SetXTitle("x-position (cm)");
  h10->SetXTitle("x-position (cm)");
  h11->SetXTitle("x-position (cm)");
  h12->SetXTitle("Energy (MeV)");
  h9->SetYTitle("y-position (cm)");
  h10->SetYTitle("time-of-flight (ns)");
  h11->SetYTitle("Energy (MeV)");
  h12->SetYTitle("time-of-flight (ns)");

  for(int i=9;i<13;i++){
    TString hname="h";
    hname+=i;
    TH2F * hInput=(TH2F *) gROOT->FindObject(hname.Data());
    hInput->SetMarkerStyle(mark_style);
    hInput->SetMarkerSize(mark_size);
    hInput->SetMarkerColor(3);
  }
  
  //cFit->Divide(2,2);
  cFit->cd(1);
  t3->Draw("y:x>>h9","","same");
  cFit->cd(2);
  t3->Draw("t:x>>h10","","same");
  cFit->cd(3);
  t3->Draw("e:x>>h11","","same");
  cFit->cd(4);
  t3->Draw("t:e>>h12","","same");
  h9->Write("");
  h10->Write("");
  h11->Write("");
  h12->Write("");

  cFit->cd(2);
  t3->Draw("t-pol2->GetParameter(2)*theta*theta-pol2->GetParameter(1)*theta:x>>h10","","same");
  h2->Draw("same");
  
  cFit->cd(4);
  t3->Draw("t-pol2->GetParameter(2)*theta*theta-pol2->GetParameter(1)*theta:e>>h12","","same");
}

  TTree *t4 = new TTree();
void mca2root4(TString fname="output4.mca")
{//assumes mca2root.C has been run
  
  //t4->ReadFile(fname,"x/F:y/F:e/F:q/F:t/F");
  t4->ReadFile(fname,"x/F:y/F:e/F:q/F:t/F:theta/F:phi/F");
  t4->Write("");
  
  for(int i=13;i<17;i++){
    TString hname="h";
    hname+=i;
    if ((TH2F *) gROOT->FindObject(hname.Data())) {
      gROOT->FindObject(hname.Data())->Delete();  
      printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
    }
  }
  
  TH2F *h13;
  TH2F *h14;
  TH2F *h15;
  TH2F *h16;
  
  h13 = new TH2F("h13","",bins,xmin,xmax,bins,ymin,ymax);
  h14 = new TH2F("h14","",bins,xmin,xmax,bins,tmin,tmax);
  h15 = new TH2F("h15","",bins,xmin,xmax,bins,emin,emax);
  h16 = new TH2F("h16","",bins,emin,emax,bins,tmin,tmax);

  h13->SetXTitle("x-position (cm)");
  h14->SetXTitle("x-position (cm)");
  h15->SetXTitle("x-position (cm)");
  h16->SetXTitle("Energy (MeV)");
  h13->SetYTitle("y-position (cm)");
  h14->SetYTitle("time-of-flight (ns)");
  h15->SetYTitle("Energy (MeV)");
  h16->SetYTitle("time-of-flight (ns)");

for(int i=13;i<17;i++){
    TString hname="h";
    hname+=i;
    TH2F * hInput=(TH2F *) gROOT->FindObject(hname.Data());
    hInput->SetMarkerStyle(mark_style);
    hInput->SetMarkerSize(mark_size);
    hInput->SetMarkerColor(6);
   
  }
  
  //cFit->Divide(2,2);
  cFit->cd(1);
  t4->Draw("y:x>>h13","","same");
  cFit->cd(2);
  t4->Draw("t:x>>h14","","same");
  cFit->cd(3);
  t4->Draw("e:x>>h15","","same");
  cFit->cd(4);
  t4->Draw("t:e>>h16","","same");

  h13->Write("");
  h14->Write("");
  h15->Write("");
  h16->Write("");
}

void fittest(TString fun="pol2",double sigma1=1)
{
  gRandom = new TRandom3();
  gRandom->SetSeed(0);

  hname="data"; 
  if ((TH2F *) gROOT->FindObject(hname)) {
    gROOT->FindObject(hname)->Delete();  
    printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
  }

  Float_t xmin=0;
  Float_t xmax=100;
  Float_t xbins=1000;
  
  TH2F * hist = new TH2F(hname, "hist", xbins, xmin, xmax,100, 0.0, 100.0);

  f1 = new TF1("f1",fun.Data(),xmin,xmax);

  Float_t x;
  
  for (int j = 0; j < 3*xbins; ++j) {
    x=gRandom->Uniform(xmin,xmax);
    for (int i = 0; i < 1e2; ++i)
      hist->Fill(x,gRandom->Gaus(f1(x), sigma1));
  }
  hist->Draw("col");

}

