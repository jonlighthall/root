/*----------------------PostScript "pretty-print" page width (97 columns)----------------------*/
/* Program: fit.cc
 *       Created  by Jon Lighthall
 *       Adapted from "linefit.cc" by Jack Winkelbauer
 *       Adapted from "peakfit.cc" by Alan Wuosmaa
 *       Adapted from "util.cc" by Alan Wuosmaa
 *       Adapted from "bkffit.cc" by Alan Wuosmaa
 * Purpose: 
 *       A package of utilities intended for use with HELIOS data analysis.  The utilities are
 *       divided into five groups: 1). Extension Utilities
 *                                 2). Plotting Utilities
 *                                 3). Fitting Utilities
 *                                 4). TSpectrum Utilities
 *                                 5). File Utilities
 * Requires:
 *       util.cc
 */
#include <iostream.h>
TH1F *hProj=0;
TH1F *hProf=0;
TH1F *hfit=0;
TH1F *hBkg=0;
TH1F *hResult=0;

TH2F *hInput=0;
TH2F *hOutput=0;
TH2F *h;

Float_t array[24][11];
Int_t det=0;
TString hname;

/* 1). Extension Utilities-----------------------
 * Extends the function of pre-existing utilities. 
 */
void add2(Char_t *histin1, Char_t *histin2, Char_t *histout=0, Float_t scale1=1.0, Float_t scale2=1.0)
{//Adds two 2D histograms.  If no output is given, a new histogram is made with copy2().
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  TH2F *hist1=(TH2F *) gROOT->FindObject(histin1);
  TH2F *hist2=(TH2F *) gROOT->FindObject(histin2);
  if(!gROOT->FindObject(histin1)) printf("Histogram \"%s\" not found.\n",histin1);
  if(!gROOT->FindObject(histin2)) printf("Histogram \"%s\" not found.\n",histin2);
  if(!histout){
    printf("No output histogram given.\n");
    hname=histin1; 
    hname+="_copy";     
    if ((TH2F *) gROOT->FindObject(hname)) {
      printf("Default output, %s, already exists",hname.Data());
      if(hname==histin2){
	printf(".  Creating new histogram...\n",hname.Data(),histin2);
	copy2(histin2,0);
      }
      else{
	printf("; it will be overwritten.\n",hname.Data());
    }
  }
    else{
      copy2(histin1,0);
    } 
    TH2F *hist3=(TH2F *) gROOT->FindObject(hname.Data());    
  }
  else{
     printf("Output histogram %s",histout);  
     if(!((TH2F *) gROOT->FindObject(histout))){
       printf(" does not exist.  Cloning %s...\n",histin1);
       hist1->Clone(histout);
     }
     else
       printf(" already exists; it will be overwritten.\n");
     TH2F *hist3=(TH2F *) gROOT->FindObject(histout);
  }  
  cFit->Clear();
  cFit->Divide(1,3);
  cFit->cd(1);
  hist1->Draw("COL");
  cFit->cd(2);
  hist2->Draw("COL");
  hist3->Add(hist1,hist2,scale1,scale2); 
  cFit->cd(3);
  hist3->Draw("COL");
  }

void dr(Char_t *histname,Float_t xmin=-999999.,Float_t xmax=999999.,Float_t ymin=-999999.,
	Float_t ymax=999999., Bool_t clear=1)
{//Extension to draw() and draw2().  
 //Accepts either 1-, 2- or 3-D histograms as input, then via the whatis() command, draws the
 //histogram.  Makes the cFit canvas, if it is not present, and clears it if it is.
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();    
  if(clear) cFit->Clear();
  if(gROOT->FindObject(histname)){
      if(whatis(histname,0)==1||whatis(histname,0)==2){
	draw(histname,"",xmin,xmax);//draw() takes the draw option as the second argument.
      }
      if(whatis(histname,0)==3||whatis(histname,0)==4)
	draw2(histname,xmin,xmax,ymin,ymax);
      if(whatis(histname,0)==5)
	gROOT->FindObject(histname)->Draw();
    }
  else
    printf("Histogram \"%s\" not recognized.\n",histname);
}

void fill1(Char_t *filename,Char_t *histname,Int_t reset=1)
{//Extension to fillhist2(), includes reset option.  
 //Fills a 1-dimensional histogram from a text file.
 //File is to be formatted as x-value, y-value, weight.
  Float_t x,y,w;
  TH1F *hist1=(TH1F *) gROOT->FindObject(histname);
  if (reset){
    printf("Resetting histogram \"%s\"\n",histname);
    hist1->Reset();
  }
  ifstream infile(filename);
  while (infile >> x) {
    //infile >> y;
    infile >> w;
    hist1->Fill(x,w);
  }
  hist1->Draw();
}

void fill2(Char_t *filename,Char_t *histname,Int_t reset=1)
{//Extension to fillhist2(), includes reset option.
 //Fills a 2-dimensional histogram from a text file.
 //File is to be formatted as x-value, y-value, weight.
  Float_t x,y,w;
  TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
  if (reset){
    printf("Resetting histogram \"%s\"\n",histname);
    hist2->Reset();
  }
  ifstream infile(filename);
  while (infile >> x) {
    infile >> y;
    infile >> w;
    hist2->Fill(x,y,w);
  }
  hist2->Draw("col2");
}

void opjx(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Int_t col=2)
{//"overlay ProjectionX" - extension to pjx()
  hInput=(TH2F *) gROOT->FindObject(histin);
  hname=histin;
  hname+="_px";
  if(maxpf==minpf){
    minpf=hInput->GetYaxis()->GetXmin();
    maxpf=hInput->GetYaxis()->GetXmax();
  }
  minpf=hInput->GetYaxis()->FindBin(minpf);
  maxpf=hInput->GetYaxis()->FindBin(maxpf);
  
  hInput->ProjectionX(hname,minpf,maxpf);

  hProj=(TH1F *) gROOT->FindObject(hname.Data());
  hProj->SetLineColor(col);
  hProj->Draw("same");
}

void pjxy(Char_t *histin)
{//plots the x- and y-projections of a 2D histogram
if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
cFit->Clear();
  cFit->Divide(1,2);
  cFit->cd(1);
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  //  hInput->Draw("COL2");

  hname=histin;
  hname+="_px"; 
  hInput->ProjectionX(hname);
  hProj=(TH1F *) gROOT->FindObject(hname.Data());
  hProj->Draw();
  
  cFit->cd(2);
  hname=histin;
  hname+="_py"; 
  hInput->ProjectionY(hname);
  hProj=(TH1F *) gROOT->FindObject(hname.Data());
  hProj->Draw();
}

void opjy(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Int_t col=2)
{//
  hInput=(TH2F *) gROOT->FindObject(histin);
  hname=histin;
  hname+="_py";
  if(maxpf==minpf){
    minpf=hInput->GetXaxis()->GetXmin();
    maxpf=hInput->GetXaxis()->GetXmax();
  }
  minpf=hInput->GetXaxis()->FindBin(minpf);
  maxpf=hInput->GetXaxis()->FindBin(maxpf);
  hInput->ProjectionY(hname,minpf,maxpf);

  hProj=(TH1F *) gROOT->FindObject(hname.Data());
  hProj->SetLineColor(col);
  hProj->Draw("same");
}

/* 2). Plotting Utilities-------------------------
 * Utilities for viewing, copying, shifting, and scaling histograms.
 */
void mkCanvas2(Char_t* cvname="cFit",Char_t *cvtitle="cFit",Int_t ww=632,Int_t wh=646)
{
  TCanvas * cFit=new TCanvas(cvname,cvtitle,0,0,ww,wh);
  if(!(cFit->GetShowEventStatus()))cFit->ToggleEventStatus();
  if(!(cFit->GetShowToolBar()))cFit->ToggleToolBar();
}

void doprint2(Char_t * cnvname="cFit", Char_t * filename="print.ps", Char_t * pr="f1-phaser")
{//prints to F-150
  const char cmd[255];
  char * fl = reinterpret_cast<char *>(filename);
  TCanvas *thecanvas=(TCanvas *)gROOT->FindObject(cnvname);
  thecanvas->SaveAs(filename);
  sprintf(cmd,"lpr -P %s %s",pr,fl);
  gSystem->Exec(cmd);
  cout<<fl<<" was sent to printer "<<pr<<endl;
}

void plotall(Char_t *histin,Float_t minX=0,Float_t maxX=0,Float_t minY=0,Float_t maxY=0,Int_t scale=1)
{//script to replace all of the macros in helios_plottools.cc
  Int_t col=0,row=0;
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2("cFit","cFit",1272,695);
  Int_t no=0,no1=0, no2=0; //number of histograms with given name
  for(Int_t i=1;i<25;++i){
    hname=histin;
    hname+=i;//had to change from "=hname+i" to "+=i" to work on ROOT 5.26
    if(gROOT->FindObject(hname.Data())) no++;
    //if(gROOT->FindObject(hname.Data())->InheritsFrom("TH1F")) no1++;
    //if(gROOT->FindObject(hname.Data())->InheritsFrom("TH2F")) no2++;
  }
  printf("Histograms with name: %d.  1D: %d.  2D: %d.\n",no,no1,no2);
  cFit->Clear();
  
  if(no!=0){
    printf("Plotting %2d Histograms...\n",no);  
    if(no>6)
      col=6;
    else
      col=ceil(no/2.);
    row=ceil(no/(Float_t)col);
    
    cFit->Divide(col,row);
    printf("Dividing Canvas as %d,%d\n",col,row);
    for(int i=1;i<(no+1);++i){
      cFit->cd(i);
      hname=histin;
      hname+=i;
      
      if(gROOT->FindObject(hname.Data())->InheritsFrom("TH2F")){

	hInput=(TH2F*)gROOT->FindObject(hname.Data());
	//	->InheritsFrom("TH1F"))
      
	if(maxX==minX){
	  minX=hInput->GetXaxis()->GetXmin();
	  maxX=hInput->GetXaxis()->GetXmax();
	}
	
	if(maxY==minY){
	  minY=hInput->GetYaxis()->GetXmin();
	  maxY=hInput->GetYaxis()->GetXmax();
	} 
	if(scale==1){
	  hInput->SetAxisRange(minX,maxX,"X");
	  hInput->SetAxisRange(minY,maxY,"Y");
	}
	else{
	  hInput->SetAxisRange(-1,-1,"X");
	  hInput->SetAxisRange(-1,-1,"Y");
	}
	hInput->Draw("COL2");
	   }
	   else{
	     hProj=(TH1F*)gROOT->FindObject(hname.Data());

	     if(maxX==minX){
	       minX=hProj->GetXaxis()->GetXmin();
	       maxX=hProj->GetXaxis()->GetXmax();
	     }
	
	     if(scale==1){
	       hProj->SetAxisRange(minX,maxX,"X");
	     }
	     else{
	       hProj->SetAxisRange(-1,-1,"X");
	     	     }
	     hProj->Draw("");

	     }
    }
}
  else{
    hname=histin;
    dr(hname);   
    //hInput=(TH2F*)gROOT->FindObject(hname.Data());
    //hInput->Draw("COL2");
  }
}
  


void plotallpjx(Char_t *histin,Float_t minY=0,Float_t maxY=0,Int_t scale=1)
{//note, mins and maxs not used.  consider matching input to pjx()
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2("cFit","cFit",1272,695);
  Int_t no=0;
  Int_t col=0,row=0;
  for(int i=1;i<25;++i){
    hname=histin;
    hname+=i;
    if((TH2F*)gROOT->FindObject(hname.Data()))
      no++;
  }
  printf("Plotting %2d Histogram Projections...\n",no);
  if(no>6)
    col=6;
  else
    col=ceil(no/2.);
  row=ceil(no/(Float_t)col);
  cFit->Clear();
  cFit->Divide(col,row);
  for(int i=1;i<(no+1);++i){
    cFit->cd(i);
    hname=histin;
    hname+=i;
    
    hInput=(TH2F*)gROOT->FindObject(hname.Data());
    if(maxY==minY){
      minY=hInput->GetYaxis()->GetXmin();
      maxY=hInput->GetYaxis()->GetXmax();
    }
    minY=hInput->GetYaxis()->FindBin(minY);
    maxY=hInput->GetYaxis()->FindBin(maxY);

    hname=hname+"_px";
    hInput->ProjectionX(hname,minY,maxY);
    
    hProj=(TH1F *) gROOT->FindObject(hname.Data());
    hProj->Draw();
  }
}

void plotallpjy(Char_t *histin,Float_t minX=0,Float_t maxX=0,Int_t scale=1)
{
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2("cFit","cFit",1272,695);
  Int_t no=0;
  Int_t col=0,row=0;
  for(int i=1;i<25;++i){
    hname=histin;
    hname+=i;
    if((TH2F*)gROOT->FindObject(hname.Data()))
      no++;
  }
  printf("Plotting %2d Histogram Projections...\n",no);
   if(no>6)
    col=6;
  else
    col=ceil(no/2.);
  row=ceil(no/(Float_t)col);
  cFit->Clear();
  cFit->Divide(col,row);
  for(int i=1;i<(no+1);++i){
    cFit->cd(i);
    hname=histin;
    hname+=i;
    
    hInput=(TH2F*)gROOT->FindObject(hname.Data());

    if(maxX==minX){
      minX=hInput->GetXaxis()->GetXmin();
      maxX=hInput->GetXaxis()->GetXmax();
    }
    minX=hInput->GetXaxis()->FindBin(minX);
    maxX=hInput->GetXaxis()->FindBin(maxX);

    hname=hname+"_py";
    //hInput->ProjectionY(hname,minX,maxX);
     hInput->ProjectionY(hname);

    hProj=(TH1F *) gROOT->FindObject(hname.Data());

    hProj->Draw();
  }
}

void setrange(Char_t *histin,Float_t minX=0,Float_t maxX=0,Float_t minY=0,Float_t maxY=0,Int_t scale=1)
{//doesn't work!
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  if(maxY==minY){
    minY=hInput->GetYaxis()->GetXmin();
    maxY=hInput->GetYaxis()->GetXmax();
    hInput->SetAxisRange(minY ,maxY ,"Y");
  }
  
  if(maxX==minX){
    minX=hInput->GetXaxis()->GetXmin();
    maxX=hInput->GetXaxis()->GetXmax();
    hInput->SetAxisRange(minX,maxX,"X");
  }
  
  if(scale==1){
    hInput->SetAxisRange(minX,maxX,"X");
    hInput->SetAxisRange(minY ,maxY ,"Y");
  }
  else{
    hInput->GetYaxis()->UnZoom();
    hInput->GetYaxis()->UnZoom();
  }
}

void setscale(Char_t *histin,Float_t minX=0,Float_t maxX=0,Float_t minY=0,Float_t maxY=0,Int_t scale=1){
  //this program doesn't work because the variable go out of scope!
 Float_t xmin,xmax; 
 Float_t ymin,ymax; 

 TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
 
 xmax=hInput->GetXaxis()->GetXmax();
 ymax=hInput->GetYaxis()->GetXmax();
 xmin=hInput->GetXaxis()->GetXmin();
 ymin=hInput->GetYaxis()->GetXmin();

 if(scale==1){
   hInput->SetAxisRange(minX,maxX,"X");
   hInput->SetAxisRange(minY,maxY,"Y");
 }
 else{
   hInput->SetAxisRange(xmin,xmax,"X");
   hInput->SetAxisRange(ymin,ymax,"Y");
 }
 
 hInput->Draw("COL2");
}

void mkhist(Char_t *histin="h", Int_t bins=3, Float_t size=10)
{//creates a small histogram to test copy2()
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();    
  hname=histin;
  if ((TH2F *) gROOT->FindObject(hname.Data())) {
    gROOT->FindObject(hname)->Delete();  
    printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
  } 

  h = new TH2F(hname.Data(),"Small Histogram",bins,0.,size,bins,0.,size);
  h->Fill(3,3,1);
  h->Fill(8,8,1);
  h->Fill(-5,5,1);//underfill
  h->Fill(12,8,1);//overfill
  h->Fill(8,8,1);//second entry in same bin
  //  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  h->Draw("colz");
}

void copy2(Char_t *histin, Float_t minz=0)
{//copies a 2D histogram with zero suppression
if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();    
 Float_t xmin,xmax; 
 Float_t ymin,ymax; 
 Int_t xbin;
 Int_t ybin,entry=0;
 Float_t x,y,z;
  
 TH2F * hInput=(TH2F *) gROOT->FindObject(histin);

 xbin=hInput->GetXaxis()->GetNbins();
 ybin=hInput->GetYaxis()->GetNbins();

 xmax=hInput->GetXaxis()->GetXmax();
 ymax=hInput->GetYaxis()->GetXmax();
 xmin=hInput->GetXaxis()->GetXmin();
 ymin=hInput->GetYaxis()->GetXmin();

 hname=histin;
 hname+="_copy"; 
 Char_t *htitle = hInput->GetTitle();
 printf("Output histogram is \"%s\"\n",hname.Data());

 if ((TH2F *) gROOT->FindObject(hname)) {
   gROOT->FindObject(hname)->Delete();  
   printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
 }

 // printf("Output histogram is constructed as:\n TH2F(\"%s\",\"%s\",%d,%1.0f,%1.0f,%d,%1.0f,%1.0f)\n",hname.Data(),htitle,xbin,xmin,xmax,ybin,ymin,ymax);
 TH2F * hOutput=new  TH2F(hname,htitle,xbin,xmin,xmax,ybin,ymin,ymax);

 for(int i=0;i<(xbin+2);i++){
   for(int j=0;j<(ybin+2);j++){
     //Note: The 0 bin contains the underflow, so the loop starts at 0;
     //      and the max_bin+1 contains overflow, so the loop terminates at max_bin+2
     x=hInput->GetXaxis()->GetBinCenter(i);
     y=hInput->GetYaxis()->GetBinCenter(j);
     z=hInput->GetBinContent(i,j);
     //printf("i=%2d, j=%2d, z=%2.0f \n",i,j,z);
     if(z!=0){
       if((Int_t)z-z)printf("Warning!  The content of bin (%d,%d) is not an integer! (%f)\n",i,j,z);
       for(int k=0;k<(z-minz);k++){//Each bin is filled with a for loop so the number of entries is the same in the copied histogram (for minz=0).  Otherwise, the number of entries is equal to the number of non-zero bins.
	 hOutput->Fill(x,y,1);
	 // entry++;
	 //	 printf("Entry %2d is %2f,%2f,%2.0f\n",entry,x,y,1);
       }
     }
   }
 }
 cFit->Clear();
 cFit->Divide(1,2);
 cFit->cd(1);
 hInput->Draw("colz");
 cFit->cd(2);
 hOutput->Draw("colz");
}

void shiftx2(Char_t *histin, Float_t shift=0, Int_t plot=2)
{//copies a 2D histogram with a given x-offset.
  if(plot>0)if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();    

 Float_t xmin,xmax; 
 Float_t ymin,ymax; 
 Int_t xbin;
 Int_t ybin,entry=0;
 Float_t x,y,z;
 if(gROOT->FindObject(histin))  
   TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
 else
    printf("Histogram \"%s\" not recognized\n",histin);

 printf("Notice: Overflow and underflow are neglected, so the\n        total number of entries may not match!\n");

 xbin=hInput->GetXaxis()->GetNbins();
 ybin=hInput->GetYaxis()->GetNbins();

 xmax=hInput->GetXaxis()->GetXmax();
 ymax=hInput->GetYaxis()->GetXmax();
 xmin=hInput->GetXaxis()->GetXmin();
 ymin=hInput->GetYaxis()->GetXmin();

 hname=histin;
 hname+="_shift"; 
 Char_t *htitle = hInput->GetTitle();
 printf("Output histogram is \"%s\"\n",hname.Data());

 if ((TH2F *) gROOT->FindObject(hname)) {
   gROOT->FindObject(hname)->Delete();  
   printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
 }
 // printf("Output histogram is constructed as:\n TH2F(\"%s\",\"%s\",%d,%1.0f,%1.0f,%d,%1.0f,%1.0f)\n",hname.Data(),htitle,xbin,xmin,xmax,ybin,ymin,ymax);
 TH2F * hOutput=new  TH2F(hname,htitle,xbin,xmin,xmax,ybin,ymin,ymax);

 for(int i=1;i<(xbin+1);i++){
   for(int j=1;j<(ybin+1);j++){
     //Note: The 0 bin contains the underflow, so the loop starts at 1;
     //      and the max_bin+1 contains overflow, so the loop terminates at max_bin+1
     x=hInput->GetXaxis()->GetBinCenter(i);
     y=hInput->GetYaxis()->GetBinCenter(j);
     z=hInput->GetBinContent(i,j);
     //printf("i=%2d, j=%2d, z=%2.0f \n",i,j,z);
     if(z!=0){
       for(int k=0;k<(z);k++){//Each bin is filled with a for loop so the number of entries is the same in the copied histogram (for minz=0).  Otherwise, the number of entries is equal to the number of non-zero bins.
	 hOutput->Fill(x+shift,y,1);
	 //	 entry++;
	 //	 printf("Entry %2d is %2f,%2f,%2.0f\n",entry,x,y,1);
       }
     }
   }
 }
 if(plot==2){
   cFit->Clear();
   cFit->Divide(1,2);
   cFit->cd(1);
   hInput->Draw("colz");
   cFit->cd(2);
 }
 if(plot==1)cFit->Clear();
 if(plot>0) hOutput->Draw("colz");
}

void slopex(Char_t *histin, Float_t slope=1, Float_t offset=0,
	    Bool_t scale=1, Float_t min=0, Float_t max=0)
{ //Copies and scales a 1D histogram with given slope and offset.
  //"scale" sets whether the bin size is scaled.
  //The axis range may be set manually with min and max.
if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();    

 Float_t xmin,xmax; 
 Int_t xbin_in,xbin_out;
 Float_t x,y,z;
 TH1F * hProj=(TH1F *) gROOT->FindObject(histin);
 Float_t xwidth_in;

 xbin_in=hProj->GetXaxis()->GetNbins();
 xbin_out=xbin_in;
 xmax=hProj->GetXaxis()->GetXmax();
 xmin=hProj->GetXaxis()->GetXmin(); 
 xwidth_in=hProj->GetBinWidth(0);

 if(min>=max){//i.e. no/bad range given
   if(slope<0){
     Float_t min=xmin;
     //  if(scale){
       xmin=(xmax-offset)/slope;
       xmax=(min-offset)/slope;
       // }
       // else{
       // xmin=xmax*-1;
       // xmax=min*-1;
       // }
   }
   else{
     // if(scale){
       xmax=(xmax-offset)/slope;
       xmin=(xmin-offset)/slope;
       // }
   }
 if(!scale)
   xbin_out=abs(xbin_in*slope);
 }
 else{
   if(scale)
     xbin_out=abs(xbin_in*slope*(max-min)/(xmax-xmin));
   xmin=min;
   xmax=max;
 }

 hname=histin;
 hname+="_slope"; 
 Char_t *htitle = hProj->GetTitle();
 if ((TH1F *) gROOT->FindObject(hname)) {
   gROOT->FindObject(hname)->Delete();  
   printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
 }
 printf("Output histogram is constructed as:\n TH1F(\"%s\",\"%s\",%d,%1.0f,%1.0f)\n",         hname.Data(),htitle,xbin_out,xmin,xmax);
 TH1F * hResult=new  TH1F(hname,htitle,xbin_out,xmin,xmax);

 for(int i=1;i<(xbin_in+1);i++){
     x=hProj->GetXaxis()->GetBinCenter(i);
     z=hProj->GetBinContent(i);
     if(z!=0){
       for(int k=0;k<(z);k++){
	 hResult->Fill((x-offset)/slope,1);
       }
     }
 }
 cFit->Clear();
 cFit->Divide(1,2);
 cFit->cd(1);
 hProj->Draw();
 cFit->cd(2);
 hResult->Draw();
 Float_t xwidth_out=hResult->GetBinWidth(0);
 if(offset<xwidth_out&&offset!=0)
   printf("offset smaller than bin size\n");

 if((fabs(xwidth_in-fabs(xwidth_out*slope))/xwidth_in)>0){
   printf("Warning: Bin width miss-match!\n");
   printf("Input bin width is  %f\n",xwidth_in);
   printf("Output bin width is %f",xwidth_out);
   printf(" (Scaled bin width is %f)\n",xwidth_out*slope);
 }
}

void slopexy(Char_t *histin, Float_t slopex=1, Float_t offsetx=0,
	     Float_t slopey=1, Float_t offsety=0)
{ //Copies and scales a 1D histogram with given slope and offset.
  //"scale" sets whether the bin size is scaled. (disabled)
  //The axis range may be set manually with min and max. (disabled)
  Bool_t scale=1; Float_t min=0; Float_t max=0;
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();    
  
  Float_t xmin,xmax; 
  Float_t ymin,ymax; 
  Int_t xbin_in,xbin_out;
  Int_t ybin_in,ybin_out;
  Int_t entry=0;
  Float_t x,y,z;
  Float_t xwidth_in;
  Float_t ywidth_in;
  
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  
  xbin_in=hInput->GetXaxis()->GetNbins();
  xbin_out=xbin_in;
  xmax=hInput->GetXaxis()->GetXmax();
  xmin=hInput->GetXaxis()->GetXmin(); 
  xwidth_in=hInput->GetXaxis()->GetBinWidth(0);
  
  ybin_in=hInput->GetYaxis()->GetNbins();
  ybin_out=ybin_in;
  ymax=hInput->GetYaxis()->GetXmax();
  ymin=hInput->GetYaxis()->GetXmin(); 
  ywidth_in=hInput->GetYaxis()->GetBinWidth(0);
  
  if(slopex<0){
    Float_t min=xmin;
    xmin=(xmax-offsetx)/slopex;
    xmax=(min-offsetx)/slopex;
  }
  else{
    xmax=(xmax-offsetx)/slopex;
    xmin=(xmin-offsetx)/slopex;
  }

  if(slopey<0){
    Float_t min=xmin;
    ymin=(ymax-offsety)/slopey;
    ymax=(min-offsety)/slopey;
  }
  else{
    ymax=(ymax-offsety)/slopey;
    ymin=(ymin-offsety)/slopey;
  }
    
  hname=histin;
  hname+="_slope"; 
  Char_t *htitle = hInput->GetTitle();
  if ((TH2F *) gROOT->FindObject(hname)) {
    gROOT->FindObject(hname)->Delete();  
    printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
  }
  printf("Output histogram is constructed as:\n TH1F(\"%s\",\"%s\",%d,%1.0f,%1.0f,%f,%f,%f)\n",         hname.Data(),htitle,xbin_out,xmin,xmax,ybin_out,ymin,ymax);
  TH2F * hOutput=new  TH2F(hname,htitle,xbin_out,xmin,xmax,ybin_out,ymin,ymax);
  
  for(int i=1;i<(xbin_in+1);i++){
    for(int j=0;j<(ybin_in+1);j++){
      x=hInput->GetXaxis()->GetBinCenter(i);
      y=hInput->GetYaxis()->GetBinCenter(j);
      z=hInput->GetBinContent(i,j);
      if(z!=0){
	for(int k=0;k<(z);k++){
	  hOutput->Fill((x-offsetx)/slopex,(y-offsety)/slopey,1);
	}
      }
    }
  }
  cFit->Clear();
  cFit->Divide(1,2);
  cFit->cd(1);
  hInput->Draw();
  cFit->cd(2);
  hOutput->Draw();
 
  //output warnings
 Float_t xwidth_out=hOutput->GetXaxis()->GetBinWidth(0);
  if(offsetx<xwidth_out&&offsetx!=0)
    printf("offset smaller than bin size\n");
  
  if((fabs(xwidth_in-fabs(xwidth_out*slopex))/xwidth_in)>0){
    printf("Warning: Bin width miss-match!\n");
    printf("Input bin width is  %f\n",xwidth_in);
    printf("Output bin width is %f",xwidth_out);
    printf(" (Scaled bin width is %f)\n",xwidth_out*slopex);
  }
}

void reflect(Char_t *histin, Float_t size=1)
{//reflects a 2d histogram about y=x and overlays it onto itself
if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();    
 Float_t xmin,xmax; 
 Float_t ymin,ymax; 
 Int_t xbin;
 Int_t ybin;
 Float_t x,y,z;
  
 TH2F * hInput=(TH2F *) gROOT->FindObject(histin);

 xbin=hInput->GetXaxis()->GetNbins();
 ybin=hInput->GetYaxis()->GetNbins();
 if(xbin!=ybin)
   printf("Warning: histogram not symmetric about y=x\n");
 xmax=hInput->GetXaxis()->GetXmax();
 ymax=hInput->GetYaxis()->GetXmax();
 xmin=hInput->GetXaxis()->GetXmin();
 ymin=hInput->GetYaxis()->GetXmin();

 hname=histin;
 hname+="_reflect"; 
 Char_t *htitle = hInput->GetTitle();
 printf("Output histogram is \"%s\"\n",hname.Data());

 if ((TH2F *) gROOT->FindObject(hname)) {
   gROOT->FindObject(hname)->Delete();  
   printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
 }

 // printf("Output histogram is constructed as:\n TH2F(\"%s\",\"%s\",%d,%1.0f,%1.0f,%d,%1.0f,%1.0f)\n",hname.Data(),htitle,xbin,xmin,xmax,ybin,ymin,ymax);
 TH2F * hOutput=new  TH2F(hname,htitle,ybin,ymin,ymax,xbin,xmin,xmax);

 for(int i=0;i<(xbin+2);i++){
   for(int j=0;j<(ybin+2);j++){
     //Note: The 0 bin contains the underflow, so the loop starts at 0;
     //      and the max_bin+1 contains overflow, so the loop terminates at max_bin+2
     x=hInput->GetXaxis()->GetBinCenter(i);
     y=hInput->GetYaxis()->GetBinCenter(j);
     z=hInput->GetBinContent(i,j);
     //printf("i=%2d, j=%2d, z=%2.0f \n",i,j,z);
     if(z!=0){
       if((Int_t)z-z)printf("Warning!  The content of bin (%d,%d) is not an integer! (%f)\n",i,j,z);
       for(int k=0;k<z;k++){//Each bin is filled with a for loop so the number of entries is the same in the copied histogram (for minz=0).  Otherwise, the number of entries is equal to the number of non-zero bins.
	 hOutput->Fill(y,x,1);
       }
     }
   }
 }
 cFit->Clear();
 if(size!=1){
   hInput ->SetMarkerStyle(21);
   hOutput->SetMarkerStyle(21);
 }
 else{
   hInput ->SetMarkerStyle(1);
   hOutput->SetMarkerStyle(1);
 }
 
 hInput->SetMarkerColor(1);
 hInput->SetMarkerSize(size);
 hInput->Draw("");

 hOutput->SetMarkerColor(2);
 hOutput->SetMarkerSize(size);
 hOutput->Draw("same");
}

void scale_slope(Char_t *histin,Float_t slope=1)
{
  if(slope>1){
    printf("Scale XN (x) by %f with slopexy(\"%s\",%f)\n",slope,histin,(1/slope));
    slopexy(histin,(1/slope));
  }
  else{
    printf("Scale XF (y) by %f with slopexy(\"%s\",1,0,%f)\n",1/slope,histin,slope);
    slopexy(histin,1,0,slope);
  }
  hname=histin;
  hname+="_slope";
  reflect(hname.Data());

 FILE * outfile;
 outfile=fopen("temp.lst","w");
 printf("The contents of \"temp.lst\" is: ");
 //in order to conform to a polynomial fit output, an arbitrary scaler offset is given.
 fprintf(outfile,"%d, %g\n",2048,-slope); 
 printf("2048, %g\n",slope);
 fclose(outfile);
}

void comb1(Char_t *histin1,Char_t *histin2 )
{//copies 1D histogram with given slope and offset
if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();    

 Float_t xmin,xmax; 
 Int_t xbin;
 Float_t xmin1,xmax1,xmin2,xmax2; 
 Int_t xbin1,xbin2;
 Float_t x,y,z;
 TH1F * hProj=(TH1F *) gROOT->FindObject(histin1);
 TH1F * hInput=(TH1F *) gROOT->FindObject(histin2);

 xbin1=hProj->GetXaxis()->GetNbins();
 xmax1=hProj->GetXaxis()->GetXmax();
 xmin1=hProj->GetXaxis()->GetXmin();
 
 xbin2=hInput->GetXaxis()->GetNbins();
 xmax2=hInput->GetXaxis()->GetXmax();
 xmin2=hInput->GetXaxis()->GetXmin();
 
 xmin=xmin1;
 if(xmin2<xmin)
   xmin=xmin2;
 xmax=xmax1;
 if(xmax2>xmax)
   xmax=xmax2;
 xbin=xbin1;
 if(xbin2<xbin)
   xbin=xbin2;
 
 hname=histin1;
 hname+="_sum"; 
 Char_t *htitle = hProj->GetTitle();
 printf("Output histogram is \"%s\"\n",hname.Data());

 if ((TH1F *) gROOT->FindObject(hname)) {
   gROOT->FindObject(hname)->Delete();  
   printf("Histogram \"%s\" already exists. Deleting old histogram.\n",hname.Data());
 }
  printf("Output histogram is constructed as:\n TH1F(\"%s\",\"%s\",%d,%1.0f,%1.0f)\n",hname.Data(),htitle,xbin,xmin,xmax);
 TH1F * hResult=new  TH1F(hname,htitle,xbin,xmin,xmax);

 for(int i=1;i<(xbin1+1);i++){
       //Note: The 0 bin contains the underflow, so the loop starts at 1;
     //      and the max_bin+1 contains overflow, so the loop terminates at max_bin+1
     x=hProj->GetXaxis()->GetBinCenter(i);
  
     z=hProj->GetBinContent(i);
         if(z!=0){
       for(int k=0;k<(z);k++){//Each bin is filled with a for loop so the number of entries is the same in the copied histogram (for minz=0).  Otherwise, the number of entries is equal to the number of non-zero bins.
	 hResult->Fill(x,1);
       }
     }
 }
 
 for(int i=1;i<(xbin2+1);i++){
       //Note: The 0 bin contains the underflow, so the loop starts at 1;
     //      and the max_bin+1 contains overflow, so the loop terminates at max_bin+1
     x=hInput->GetXaxis()->GetBinCenter(i);
       z=hInput->GetBinContent(i);
         if(z!=0){
       for(int k=0;k<(z);k++){//Each bin is filled with a for loop so the number of entries is the same in the copied histogram (for minz=0).  Otherwise, the number of entries is equal to the number of non-zero bins.
	 hResult->Fill(x,1);
       }
     }
 }
 
 cFit->Clear();
 cFit->Divide(1,2);
 cFit->cd(1);
 hProj->Draw();
 cFit->cd(2);
 hResult->Draw();
}

void shiftadd2(Char_t *histin1, Float_t shift1=0, Char_t *histin2, Float_t shift2=0)
{//adds two 2D histograms with given x-offsets
  TString name1=histin1;
  TString name2=histin2;
  shiftx2(histin1,shift1);
  shiftx2(histin2,shift2);
  name1+="_shift";
  name2+="_shift";
  add2(name1,name2);
  name1+="_copy";
  dr(name1);
}

/* 4). Fitting Utilities------------------------------
 */

void fitpfx(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Float_t minfit=0,Float_t maxfit=0,Int_t ord=1,Int_t scale=1)
{//same first three parameters as pfx(), next three same as pfit() (developed independently)
  Float_t cp=0 ;  
  Float_t a=0,b=0,c=0,d=0;  
  Float_t p0=0,p1=0,p2=0,p3=0,p4=0;
  Float_t fits[10];
  hname=histin;
  for(Int_t i=0;i<hname.Length();i++){//loop added by Jack
    TString temp="";
    temp=hname(i,hname.Length()-i);
    if(temp.IsFloat())
      {
	det=temp.Atoi();
	break;
      }
  }

  printf("Detector number %d is being read.\n",det);
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
  cFit->Clear();
  cFit->SetWindowPosition(0,0);
  cFit->Divide(1,2);
  cFit->cd(1);

  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  //  printf("Input histogram is %s\n",histin);    
  if(maxpf==minpf){
    minpf=hInput->GetYaxis()->GetXmin();
    maxpf=hInput->GetYaxis()->GetXmax();
    hInput->SetAxisRange(minpf ,maxpf ,"Y");
  }
  
  if(maxfit==minfit){
    minfit=hInput->GetXaxis()->GetXmin();
    maxfit=hInput->GetXaxis()->GetXmax();
    hInput->SetAxisRange(minfit,maxfit,"X");
  }
    
  if(scale==1){
    hInput->SetAxisRange(minfit,maxfit,"X");
    hInput->SetAxisRange(minpf ,maxpf ,"Y");
  }
  else{
    hInput->GetXaxis()->UnZoom();
    hInput->GetYaxis()->UnZoom();
  }

  hInput->Draw("COL2");
  cFit->cd(2); 
  a=minpf;
  b=maxpf; 
  printf("Projection Limits are %3.2f to %3.2f\n",minpf,maxpf);  
  minpf=hInput->GetYaxis()->FindBin(minpf);
  maxpf=hInput->GetYaxis()->FindBin(maxpf);
  printf("       Fit Limits are %3.2f to %3.2f\n",minfit,maxfit);

  hname=histin;
  hname+="_pfx"; 
  hInput->ProfileX(hname,minpf,maxpf);
 
  hProf=(TH1F *) gROOT->FindObject(hname.Data());
  cFit->cd(2);
  hProf->GetXaxis()->UnZoom();
  hProf->GetYaxis()->UnZoom();
  hProf->SetAxisRange(a,b,"Y");
  hProf->SetLineColor(2);
  hProf->Draw();
  
  if(ord>9||ord<1){
    printf("Polynomial fits only valid for orders 1-9.\n");
    ord=1;
  }
  hname="pol";
  hname+=ord;
  printf("Fit function is is \"%s\"\n",hname.Data()); 
 hProf->Fit(hname,"V","",minfit,maxfit);
 hProf->SetStats(kFALSE);
 switch(ord){
 case 1://original fitpfx()
    p0=hProf->GetFunction("pol1")->GetParameter(0);
    p1=hProf->GetFunction("pol1")->GetParameter(1);
    printf("p1 = %7.3f\n",p1);  
    break;
 case 2://copied from minfit() - used to find minimum (center) of hEX plots - then generalize to fit2pfx()
   p0=hProf->GetFunction("pol2")->GetParameter(0);//"c"
   p1=hProf->GetFunction("pol2")->GetParameter(1);//"b"
   p2=hProf->GetFunction("pol2")->GetParameter(2);//"a"
   cp=(-p1/(2*p2));
   Float_t zero1=(-p1-sqrt(p1*p1-(4*p2*p0)))/(2*p2);
   Float_t zero2=(-p1+sqrt(p1*p1-(4*p2*p0)))/(2*p2);
   printf("p1 = %7.3f p2=%7.3f\n",p1,p2);
   printf("Fit has critical point = %7.3f\n",p1,p2,cp);
   Float_t slope=0;
   if((p1*p1-(4*p2*p0))>0){
     printf("Fit has roots = %.1f, %.1f\n",zero1,zero2);
     TLine *line = new TLine(0,p0,zero1,0);
     line->SetLineStyle(2);
     line->SetLineWidth(2);
     line->Draw();
     slope=(0-p0)/(zero1-0);
     printf("Linear fit has slope %f\n",slope);
   }
   else{
     printf("Fit has imaginary root(s).\n");
     TLine *line = new TLine(minfit,minfit*minfit*p2+minfit*p1+p0,maxfit,maxfit*maxfit*p2+maxfit*p1+p0);
     line->SetLineStyle(2);
     line->SetLineWidth(2);
     line->Draw();
     slope=((maxfit*maxfit*p2+maxfit*p1+p0)-(minfit*minfit*p2+minfit*p1+p0))/(maxfit-minfit);
     printf("Linear fit has slope %f\n",slope);
   }
   if(slope<-1)
     printf("Scale XN (x) by %f with slopexy(\"%s\",%f)\n",-slope,histin,(-1/slope));
   else
     printf("Scale XF (y) by %f with slopexy(\"%s\",1,0,%f)\n",-1/slope,histin,-slope);
   
   FILE * outfile;
   outfile=fopen("slope.lst","w");
   printf("The contents of \"slope.lst\" is: ");
   fprintf(outfile,"%g\n",slope);
   printf("%g\n",slope);
   fclose(outfile);
   
   break;
 case 3:
   p0=hProf->GetFunction("pol3")->GetParameter(0);
   p1=hProf->GetFunction("pol3")->GetParameter(1);
   p2=hProf->GetFunction("pol3")->GetParameter(2);
   p3=hProf->GetFunction("pol3")->GetParameter(3);
   a=p3;
   b=p2;
   c=p1;
   d=p0;
   Float_t disc=18*a*b*c*d-4*b*b*b*d+b*b*c*c-27*a*a*d*d;
   if(disc<0){
     printf("Fit has one real root\n");
     Float_t zero=0;
     //     TF1 *fit = hProf->GetFuction("pol3");
     //     printf("Fit at zero is %f\n",fit(zero));
   }
   else
     printf("Fit discriminant = %f, yeilding all real roots.\n",disc);

   break;
 case 4:
   p0=hProf->GetFunction("pol4")->GetParameter(0);
   p1=hProf->GetFunction("pol4")->GetParameter(1);
   p2=hProf->GetFunction("pol4")->GetParameter(2);
   p3=hProf->GetFunction("pol4")->GetParameter(3);
   p4=hProf->GetFunction("pol4")->GetParameter(4);
   printf("p1 = %5.0f, p2 = %5.0f, p3 = %5.0f, p4 = %5.0f\ncp = %7.3f\n",p1,p2,p3,p4);
   break;
 default:
   break; 
 }
 
 FILE * outfile;
 outfile=fopen("temp.lst","w");
 printf("The contents of \"temp.lst\" are: ");
 for(Int_t i=0;i<=ord;i++){
   fprintf(outfile,"%g, ",hProf->GetFunction(hname)->GetParameter(i));
   printf("%g ",hProf->GetFunction(hname)->GetParameter(i));
 }
 fprintf(outfile,"\n");
 printf("\n");
 fclose(outfile);
}

void fitpfy(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Float_t minfit=0,Float_t maxfit=0,Int_t scale=1)
{
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  Float_t p0=0,p1=0,p2=0;
  
  hname=histin;
  for(Int_t i=0;i<hname.Length();i++){//loop added by Jack
    TString temp="";
    temp=hname(i,hname.Length()-i);
    if(temp.IsFloat())
      {
	det=temp.Atoi();
	break;
      }
  }

  cFit->Clear();
  cFit->Divide(1,2);
  cFit->cd(1);
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
    
  if(maxpf==minpf){
    minpf=hInput->GetXaxis()->GetXmin();
    maxpf=hInput->GetXaxis()->GetXmax();
  }
  
  if(maxfit==minfit){
    minfit=hInput->GetYaxis()->GetXmin();
    maxfit=hInput->GetYaxis()->GetXmax();
  }
    
  if(scale==1){
    hInput->SetAxisRange(minfit,maxfit,"Y");
    hInput->SetAxisRange(minpf ,maxpf ,"X");
  }
  else{
    hInput->SetAxisRange(-1,-1,"X");
    hInput->SetAxisRange(-1,-1,"Y");
  }

  hInput->Draw("COL2");
  cFit->cd(2); 
  a=minpf;
  b=maxpf;   
  printf("Projection Limits are %d to %d\n",minpf,maxpf);  
  minpf=hInput->GetXaxis()->FindBin(minpf);
  maxpf=hInput->GetXaxis()->FindBin(maxpf);
  printf("Fit Limits are %d to %d\n",minfit,maxfit);
 
  hname=histin;
  hname+="_pfy"; 
  
  hInput->ProfileY(hname,minpf,maxpf);
  hProf=(TH1F *) gROOT->FindObject(hname.Data());
  cFit->cd(2);
  hProf->GetXaxis()->UnZoom();
  hProf->GetYaxis()->UnZoom();
  hProf->SetAxisRange(a,b,"Y");
  hProf->SetLineColor(2);
  hProf->Draw();
 
  hProf->Fit("pol1","V","",minfit,maxfit);
  p0=hProf->GetFunction("pol1")->GetParameter(0);
  p1=hProf->GetFunction("pol1")->GetParameter(1);
  //  p2=hProf->GetFunction("pol2")->GetParameter(2);
  // cp=(-p1/(2*p2));
  printf("p1 = %7.3f\n",p1);
}

void fitallpjx(Char_t *histin,Float_t minpj=0,Float_t maxpj=0,Float_t minfit=0,Float_t maxfit=0,Int_t ord=6,Int_t writetofile=0)
{
 if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  Float_t a=0,b=0;
  for(int i=0;i<24;i++)
    for(int j=0;j<11;j++)
      array[i][j]=0;
  
  cFit->Clear();
  cFit->Divide(1,2);
  a=minpj;
  b=maxpj;

  for(int i=0;i<24;i++){
      array[i][0]=i+1;
      //      printf("input histogram is %s\n",histin); 
      hname=histin;
      hname=hname+(i+1);
      //      printf("input histogram is %s\n",hname.Data()); 
      
      hInput=(TH2F*)gROOT->FindObject(hname.Data());
      cFit->cd(1);
  
      if(maxpj==minpj){
      minpj=hInput->GetYaxis()->GetXmin();
      maxpj=hInput->GetYaxis()->GetXmax();
    }
    minpj=hInput->GetYaxis()->FindBin(minpj);
    maxpj=hInput->GetYaxis()->FindBin(maxpj);
  
    if(maxfit==minfit){
      minfit=hInput->GetXaxis()->GetXmin();
      maxfit=hInput->GetXaxis()->GetXmax();
    }

 hInput->Draw("COL2");
      
      cFit->cd(2); 
      hname+="_px"; 
      // printf("output histogram is %s\n",hname.Data()); 
      hInput->ProjectionX(hname,minpj,maxpj);
      hProj=(TH1F *) gROOT->FindObject(hname.Data());
      hProj->Draw();
     
      hname="pol";
      hname+=ord;
      //      printf("order is %s\n",hname.Data()); 
      hProj->Fit(hname,"Q","",minfit,maxfit);
       
      for (int j=0;j<ord+1;j++){
	array[i][j+1]=hProj->GetFunction(hname)->GetParameter(j);
	//	printf("function is %s, p%d = %f\n",hname.Data(),j,array[i][j+1]);
      }
      minpj=a;
      maxpj=b;
    }
  printcaldata();
 if(writetofile){
    createfile();
    printf("Fit parameters written to file \"calibration.cal\"\n");
  }
}

void pfitallpjx(Char_t *histin,Float_t minpj=0,Float_t maxpj=0,Float_t minfit=0,Float_t maxfit=0,Int_t ord=6,Int_t writetofile=0)
  //to be used with addpositions.C
{
if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
 TString fname;
  for(int i=0;i<24;i++)
    for(int j=0;j<11;j++)
      array[i][j]=0;//initializes array to zero
  
  cFit->Clear();
  cFit->Divide(1,2);

  hname=histin;
  hname=hname+1;
            
  hInput=(TH2F*)gROOT->FindObject(hname.Data());
  if(maxpj==minpj){
	minpj=hInput->GetYaxis()->GetXmin();
	maxpj=hInput->GetYaxis()->GetXmax();
      }
      minpj=hInput->GetYaxis()->FindBin(minpj);
      maxpj=hInput->GetYaxis()->FindBin(maxpj);
      
      if(maxfit==minfit){
	minfit=hInput->GetXaxis()->GetXmin();
	maxfit=hInput->GetXaxis()->GetXmax();
      }
      
      if(i==0){
	printf("Projection range is %f to %f\n",minpj,maxpj);
	printf("Fit range is %f %f\n",minfit,maxfit);
      }


for(int i=0;i<6;i++){
      array[i][0]=i+1;
      hname=histin;
      hname=hname+(i+1);
            
      hInput=(TH2F*)gROOT->FindObject(hname.Data());
      cFit->cd(1);
       hInput->Draw("COL2");
      
      cFit->cd(2); 
      hname+="_px"; 
      hInput->ProjectionX(hname,minpj,maxpj);
      hProj=(TH1F *) gROOT->FindObject(hname.Data());
      hProj->Draw();
      fname="pol";
      fname+=ord;
      if(i==0)printf("Fit function is \"%s\"\n",fname.Data()); 
      hProj->Fit(fname,"Q","",minfit,maxfit);
       
      for (int j=0;j<ord+1;j++){
	array[i][j+1]=hProj->GetFunction(fname)->GetParameter(j);
	//	printf("function is %s, p%d = %f\n",fname.Data(),j,array[i][j+1]);
      }
    }
  printcaldata();
  if(writetofile){
    createfile();
    printf("Fit parameters written to file \"calibration.cal\"\n");
  }
  else
        printf("Fit parameters NOT written to file.\n");
}

void timefit(Char_t *histin,Float_t minE=1,Float_t maxE=12,Int_t minpf=0,Int_t maxpf=1200)
{//prototype of fitpqpfy
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
  Float_t cp=0 ;
  Int_t i=0;
  Float_t p0=0,p1=0,p2=0;
  
  cFit->Clear();
  cFit->Divide(1,2);
 
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  cFit->cd(1);
  hInput->SetAxisRange(minpf,maxpf,"X");
  hInput->SetAxisRange(minE,maxE,"Y");
  hInput->Draw("COL2");
  cFit->cd(2); 
    
  hname=histin;
  hname+="_pfy"; 
  minpf=hInput->GetXaxis()->FindBin(minpf);
  maxpf=hInput->GetXaxis()->FindBin(maxpf);
  printf("Projection Limits are t = %d to %d\n",minpf,maxpf);
  //minpf=(Int_t)((256/1200.0)*minpf);
  //maxpf=(Int_t)((256/1200.0)*maxpf);
  hInput->ProfileY(hname,minpf,maxpf);
  hProf=(TH1F *) gROOT->FindObject(hname.Data());
  cFit->cd(2);

  hProf->Draw();
 
  printf("Fit range is %6.3f to %6.3f\n",minE,maxE);
  hProf->Fit("pol2","V","",minE,maxE);
  p0=hProf->GetFunction("pol2")->GetParameter(0);
  p1=hProf->GetFunction("pol2")->GetParameter(1);
  p2=hProf->GetFunction("pol2")->GetParameter(2);
  cp=(-p1/(2*p2));
  printf("p0 = %7.3f, p1 = %7.3f, p2 = %7.3f, critical point at z = %5.3f\n",p0,p1,p2,cp);
  
    while((fabs((maxE-cp)/cp)>.03)&&i<200){

    if(cp<maxE){
      maxE-=(maxE-cp)/10;
      printf("Fit range is %6.3f MeV to %6.3f MeV, ",minE,maxE);
         
      hProf->Fit("pol2","Q","",minE,maxE);
      p0=hProf->GetFunction("pol2")->GetParameter(0);
      p1=hProf->GetFunction("pol2")->GetParameter(1);
      p2=hProf->GetFunction("pol2")->GetParameter(2);
      cp=(-p1/(2*p2));
      
      printf("Critical Point at E = %5.3f\n          p0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",cp,p0,p1,p2);
    }
     else{
       maxE+=(cp-maxE)/10;
       printf("Fit range is %6.3f MeV to %6.3f MeV, \n",minE,maxE);
     
       hProf->Fit("pol2","Q","",minE,maxE);
       p0=hProf->GetFunction("pol2")->GetParameter(0);
       p1=hProf->GetFunction("pol2")->GetParameter(1);
       p2=hProf->GetFunction("pol2")->GetParameter(2);
       cp=(-p1/(2*p2));
      
 printf("Critical Point at E = %5.3f\n           p0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",cp,p0,p1,p2);
     }
    i++;
    printf("Iter. %2d: ",i);
    if(i==100)
      printf("Loop terminated.  Could not converge.\n");
  }
}
void timefit2(Char_t *histin,Float_t minE=1,Float_t maxE=24)
{
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
  Float_t cp=0 ;
  Int_t i=0;
  Float_t p0=0,p1=0,p2=0;
  
  cFit->Clear();
  cFit->Divide(1,2);
 
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  cFit->cd(1);
  hInput->Draw("COL2");
  cFit->cd(2); 
   
  hInput->ProfileY();
  hname=histin;
  hname+="_pfy";
  hProf=(TH1F *) gROOT->FindObject(hname.Data());
  cFit->cd(2);
  hProf->Draw();
 
  printf("Fit range is %6.3f to %6.3f\n",minE,maxE);
  hProf->Fit("pol2","V","",minE,maxE);
  p0=hProf->GetFunction("pol2")->GetParameter(0);
  p1=hProf->GetFunction("pol2")->GetParameter(1);
  p2=hProf->GetFunction("pol2")->GetParameter(2);
  cp=(-p1/(2*p2));
  printf("p0 = %7.3f, p1 = %7.3f, p2 = %7.3f, critical point at z = %5.3f\n",p0,p1,p2,cp);
  
    while((fabs((((maxE+minE)/2)-cp)/cp)>.03)&&i<200){

    if(cp<((maxE+minE)/2)){
      maxE-=(((maxE+minE)/2)-cp)/10;
      printf("Fit range is %6.3f MeV to %6.3f MeV, ",minE,maxE);
         
      hProf->Fit("pol2","Q","",minE,maxE);
      p0=hProf->GetFunction("pol2")->GetParameter(0);
      p1=hProf->GetFunction("pol2")->GetParameter(1);
      p2=hProf->GetFunction("pol2")->GetParameter(2);
      cp=(-p1/(2*p2));
      
      printf("Critical Point at E = %5.3f\n          p0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",cp,p0,p1,p2);
    }
     else{
       maxE+=(((maxE+minE)/2)-maxE)/10;
       printf("Fit range is %6.3f MeV to %6.3f MeV, \n",minE,maxE);
     
       hProf->Fit("pol2","Q","",minE,maxE);
       p0=hProf->GetFunction("pol2")->GetParameter(0);
       p1=hProf->GetFunction("pol2")->GetParameter(1);
       p2=hProf->GetFunction("pol2")->GetParameter(2);
       cp=(-p1/(2*p2));
      
 printf("Critical Point at E = %5.3f\n           p0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",cp,p0,p1,p2);
     }
    i++;
    printf("Iter. %2d: ",i);
    if(i==100)
      printf("Loop terminated.  Could not converge.\n");
  }
}
void timefitg(Char_t *histin,Float_t minT=100,Float_t maxT=1200)
{//redundant
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
  cFit->Clear();
  cFit->Divide(1,2);
 
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  cFit->cd(1);
  hInput->Draw("COL2");
  cFit->cd(2); 
   
  hInput->ProjectionX();
  hname=histin;
  hname+="_px";
  hProj=(TH1F *) gROOT->FindObject(hname.Data());
  cFit->cd(2);
  hProj->Draw();
 
  hProj->Fit("gaus","V","",minT,maxT);
}

void fitpqpfy(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Float_t minfit=0,Float_t maxfit=0,Int_t scale=1)
{//"Fit piece-wise quadratic, ProfileY" - generalized from timefit()
if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
Float_t cp=0 ;
  Int_t i=0;
  Float_t p0=0,p1=0,p2=0;
  Float_t a=0,b=0,c=0;
  Float_t tol=.03;

  cFit->Clear();
  cFit->Divide(1,2);
 
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  cFit->cd(1);

   if(maxpf==minpf){
    minpf=hInput->GetXaxis()->GetXmin();
    maxpf=hInput->GetXaxis()->GetXmax();
  }
  
  if(maxfit==minfit){
    minfit=hInput->GetYaxis()->GetXmin();
    maxfit=hInput->GetYaxis()->GetXmax();
  }
    
  if(scale==1){
    hInput->SetAxisRange(minfit,maxfit,"Y");
    hInput->SetAxisRange(minpf ,maxpf ,"X");
  }
  else{
    hInput->SetAxisRange(-1,-1,"X");
    hInput->SetAxisRange(-1,-1,"Y");
  }

  hInput->Draw("COL2");
  cFit->cd(2); 
    
  hname=histin;
  hname+="_pfy"; 
  a=minpf;
  b=maxpf;
  c=maxfit;

  minpf=hInput->GetXaxis()->FindBin(minpf);
  maxpf=hInput->GetXaxis()->FindBin(maxpf);
  printf("Projection Limits are t = %d to %d\n",a,b);
  hInput->ProfileY(hname,minpf,maxpf);
  hProf=(TH1F *) gROOT->FindObject(hname.Data());
  cFit->cd(2);
hProf->GetXaxis()->UnZoom();
  hProf->GetYaxis()->UnZoom();
  hProf->SetAxisRange(a,b,"Y");
  hProf->Draw();
 
  printf("Fit range is %6.3f to %6.3f\n",minfit,maxfit);
  hProf->Fit("pol2","V","",minfit,maxfit);
  p0=hProf->GetFunction("pol2")->GetParameter(0);
  p1=hProf->GetFunction("pol2")->GetParameter(1);
  p2=hProf->GetFunction("pol2")->GetParameter(2);
  cp=(-p1/(2*p2));
  printf("p0 = %7.3f, p1 = %7.3f, p2 = %7.3f, critical point at z = %5.3f\n",p0,p1,p2,cp);
  
    while((fabs((maxfit-cp)/cp)>tol)&&i<200){

    if(cp<maxfit){
      maxfit-=(maxfit-cp)/10;
      //      printf("Fit range is %6.3f to %6.3f, ",minfit,maxfit);
      hProf->Fit("pol2","Q","",minfit,maxfit);
      p0=hProf->GetFunction("pol2")->GetParameter(0);
      p1=hProf->GetFunction("pol2")->GetParameter(1);
      p2=hProf->GetFunction("pol2")->GetParameter(2);
      cp=(-p1/(2*p2));
      //      printf("Critical Point at E = %5.3f\n          p0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",cp,p0,p1,p2);
      printf("-");
    }
     else{
       maxfit+=(cp-maxfit)/10;
       //       printf("Fit range is %6.3f to %6.3f, ",minfit,maxfit);
       hProf->Fit("pol2","Q","",minfit,maxfit);
       p0=hProf->GetFunction("pol2")->GetParameter(0);
       p1=hProf->GetFunction("pol2")->GetParameter(1);
       p2=hProf->GetFunction("pol2")->GetParameter(2);
       cp=(-p1/(2*p2));
       // printf("Critical Point at E = %5.3f\n           p0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",cp,p0,p1,p2);
       printf("+");
     }
    i++;
    // printf("Iter. %3d: ",i);
    if(i==200){
      printf("\nLoop terminated.  Could not converge.\n");
      tol+=.01;
      maxfit=c;
      hProf->Fit("pol2","Q","",minfit,maxfit);
       p1=hProf->GetFunction("pol2")->GetParameter(1);
       p2=hProf->GetFunction("pol2")->GetParameter(2);
       cp=(-p1/(2*p2));

      printf("Re-initiating loop with tolerance set to %2.0f%%\n",tol*100);

     i=0;
    }
  }
    printf("Loop Exited at Iteration %3d.\nEndpoint tolerance is %2.0f%%\nFit range is %6.3f to %6.3f\nCritical Point is %5.3f\np0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",i,tol*100,minfit,maxfit,cp,p0,p1,p2);
    printf("Emax = %5.2f, p1 = %4.0f, p2 = %7.2f\n",cp,p1,p2);
}
void fitpqpfx(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Float_t minfit=0,Float_t maxfit=0,Int_t scale=1)
{//"Fit piece-wise quadratic, ProfileX"
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2(); 
  Float_t cp=0 ;
  Int_t i=0;
  Float_t p0=0,p1=0,p2=0;
  Float_t a=0,b=0,c=0;
  Float_t tol=.03;

  cFit->Clear();
  cFit->Divide(1,2);
 
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  cFit->cd(1);
  //if no limits are given, the entire range is projected and fitted
  if(maxpf==minpf){
    minpf=hInput->GetYaxis()->GetXmin();
    maxpf=hInput->GetYaxis()->GetXmax();
  }
  if(maxfit==minfit){
    minfit=hInput->GetXaxis()->GetXmin();
    maxfit=hInput->GetXaxis()->GetXmax();
  }
  
  if(scale==1){
    hInput->SetAxisRange(minfit,maxfit,"X");
    hInput->SetAxisRange(minpf ,maxpf ,"Y");
  }
  else{
    hInput->SetAxisRange(-1,-1,"X");
    hInput->SetAxisRange(-1,-1,"Y");
  }
 
  hInput->Draw("COL2");
  cFit->cd(2); 
    
  hname=histin;
  hname+="_pfx"; 
  a=minpf;
  b=maxpf;
  c=maxfit;

  minpf=hInput->GetYaxis()->FindBin(minpf);
  maxpf=hInput->GetYaxis()->FindBin(maxpf);
  printf("Projection Limits are t = %f to %f\n",a,b);
  hInput->ProfileX(hname,minpf,maxpf);
  hProf=(TH1F *) gROOT->FindObject(hname.Data());
  cFit->cd(2);
hProf->GetXaxis()->UnZoom();
  hProf->GetYaxis()->UnZoom();
  hProf->SetAxisRange(a,b,"Y");
  hProf->Draw();
 
  printf("Fit range is %6.3f to %6.3f\n",minfit,maxfit);
  hProf->Fit("pol2","V","",minfit,maxfit);
  p0=hProf->GetFunction("pol2")->GetParameter(0);
  p1=hProf->GetFunction("pol2")->GetParameter(1);
  p2=hProf->GetFunction("pol2")->GetParameter(2);
  cp=(-p1/(2*p2));
  printf("p0 = %7.3f, p1 = %7.3f, p2 = %7.3f, critical point at z = %5.3f\n",p0,p1,p2,cp);
  
  while((fabs((maxfit-cp)/cp)>tol)&&i<200){

    if(cp<maxfit){
      maxfit-=(maxfit-cp)/10;
      //      printf("Fit range is %6.3f to %6.3f, ",minfit,maxfit);
         
      hProf->Fit("pol2","Q","",minfit,maxfit);
      p0=hProf->GetFunction("pol2")->GetParameter(0);
      p1=hProf->GetFunction("pol2")->GetParameter(1);
      p2=hProf->GetFunction("pol2")->GetParameter(2);
      cp=(-p1/(2*p2));
      printf("-");
      //      printf("Critical Point at E = %5.3f\n          p0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",cp,p0,p1,p2);
    }
    else{
      maxfit+=(cp-maxfit)/10;
      //       printf("Fit range is %6.3f to %6.3f, ",minfit,maxfit);
      hProf->Fit("pol2","Q","",minfit,maxfit);
      p0=hProf->GetFunction("pol2")->GetParameter(0);
      p1=hProf->GetFunction("pol2")->GetParameter(1);
      p2=hProf->GetFunction("pol2")->GetParameter(2);
      cp=(-p1/(2*p2));
      // printf("Critical Point at E = %5.3f\n           p0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",cp,p0,p1,p2);
      printf("+");
    }
    i++;
    //    printf("Iter. %3d: ",i);
    if(i==200){
      printf("Loop terminated.  Could not converge.\n");
      tol+=.01;
      maxfit=c;
      hProf->Fit("pol2","Q","",minfit,maxfit);
      p1=hProf->GetFunction("pol2")->GetParameter(1);
      p2=hProf->GetFunction("pol2")->GetParameter(2);
      cp=(-p1/(2*p2));
      
      printf("Re-initiating loop with tolerance set to %2.0f%%\n",tol*100); cp=-1;

      i=0;
    }
  }
  printf("
Loop Exited at Iteration %3d.\nEndpoint tolerance is %2.0f%%\nFit range is %6.3f to %6.3f\nCritical Point is %5.3f\np0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",i,tol*100,minfit,maxfit,cp,p0,p1,p2);
}

void fitpqRpfx(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Float_t minfit=0,Float_t maxfit=0,Int_t scale=1)
{"Fit piece-wise quadratic, ProfileX, from right"
if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
Float_t cp=0 ;
  Int_t i=0;
  Float_t p0=0,p1=0,p2=0;
  Float_t a=0,b=0,c=0;
  Float_t tol=.03;

  cFit->Clear();
  cFit->Divide(1,2);
 
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  cFit->cd(1);
  //if no limits are given, the entire range is projected and fitted
  if(maxpf==minpf){
    minpf=hInput->GetYaxis()->GetXmin();
    maxpf=hInput->GetYaxis()->GetXmax();
  }
  if(maxfit==minfit){
    minfit=hInput->GetXaxis()->GetXmin();
    maxfit=hInput->GetXaxis()->GetXmax();
  }
 
  if(scale==1){
    hInput->SetAxisRange(minfit,maxfit,"X");
    hInput->SetAxisRange(minpf ,maxpf ,"Y");
  }
  else{
    hInput->SetAxisRange(-1,-1,"X");
    hInput->SetAxisRange(-1,-1,"Y");
  }
 
  hInput->Draw("COL2");
  cFit->cd(2); 
    
  hname=histin;
  hname+="_pfx"; 
  a=minpf;
  b=maxpf;
  c=minfit;

  minpf=hInput->GetYaxis()->FindBin(minpf);
  maxpf=hInput->GetYaxis()->FindBin(maxpf);
  printf("Projection Limits are t = %d to %d\n",a,b);
  hInput->ProfileX(hname,minpf,maxpf);
  hProf=(TH1F *) gROOT->FindObject(hname.Data());
  cFit->cd(2);
hProf->GetXaxis()->UnZoom();
  hProf->GetYaxis()->UnZoom();
  hProf->SetAxisRange(a,b,"Y");
  hProf->Draw();
 
  printf("Fit range is %6.3f to %6.3f\n",minfit,maxfit);
  hProf->Fit("pol2","V","",minfit,maxfit);
  p0=hProf->GetFunction("pol2")->GetParameter(0);
  p1=hProf->GetFunction("pol2")->GetParameter(1);
  p2=hProf->GetFunction("pol2")->GetParameter(2);
  cp=(-p1/(2*p2));
  printf("p0 = %7.3f, p1 = %7.3f, p2 = %7.3f, critical point at z = %5.3f\n",p0,p1,p2,cp);
  
    while((fabs((minfit-cp)/cp)>tol)&&i<200){

    if(cp<minfit){
      minfit-=(minfit-cp)/10;
      //      printf("Fit range is %6.3f to %6.3f, ",minfit,maxfit);
         
      hProf->Fit("pol2","Q","",minfit,maxfit);
      p0=hProf->GetFunction("pol2")->GetParameter(0);
      p1=hProf->GetFunction("pol2")->GetParameter(1);
      p2=hProf->GetFunction("pol2")->GetParameter(2);
      cp=(-p1/(2*p2));
            printf("-"); 
      //      printf("Critical Point at X = %5.3f\n          p0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",cp,p0,p1,p2);
    }
     else{
       minfit+=(cp-minfit)/10;
       //       printf("Fit range is %6.3f to %6.3f, ",minfit,maxfit);
     
       hProf->Fit("pol2","Q","",minfit,maxfit);
       p0=hProf->GetFunction("pol2")->GetParameter(0);
       p1=hProf->GetFunction("pol2")->GetParameter(1);
       p2=hProf->GetFunction("pol2")->GetParameter(2);
       cp=(-p1/(2*p2));
             printf("+");
       // printf("Critical Point at X = %5.3f\n           p0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",cp,p0,p1,p2);
     }
    i++;
    //    printf("Iter. %3d: ",i);
        if(i==200){
      printf("\nLoop terminated.  Could not converge.\n");
      tol+=.01;
      minfit=c;
      hProf->Fit("pol2","Q","",minfit,maxfit);
      p1=hProf->GetFunction("pol2")->GetParameter(1);
      p2=hProf->GetFunction("pol2")->GetParameter(2);
      cp=(-p1/(2*p2));

      printf("Reiniating loop with tolerance set to %2.0f%%\n",tol*100);
      i=0;
	}
  }
printf("\nLoop Exited at Iteration %3d.\nEndpoint tolerance is %2.0f%%\nFit range is %6.3f to %6.3f\nCritical Point is %5.3f\np0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",i,tol*100,minfit,maxfit,cp,p0,p1,p2);
}

void getline(void)
{//shows information about a TLine
  Float_t x1,x2,y1,y2;
  Float_t slope=0,b,theta;
  x1=TLine->GetX1();
  y1=TLine->GetY1();
  x2=TLine->GetX2();
  y2=TLine->GetY2();
  slope=(y2-y1)/(x2-x1);
  b=y1-slope*x1;
  theta=atan(slope)/(4.0*atan(1.0))*180;
  if(fabs(slope)>.01) 
    printf("Slope of TLine is %6.6f\n",slope);
 else
   printf("Slope of TLine is %6.3e\n",slope);
   printf("Intercept of TLine is %6.3f\n",b);
   printf("TLine angle is %6.3f\n",theta);
}

void getpave()
{//shows the position of a TPave
Float_t a,b,c,d,e,f,g,h;
 a=TPave->GetX1();
 b=TPave->GetY1();
 c=TPave->GetX2();
 d=TPave->GetY2();
 e=TPave->GetX1NDC();
 f=TPave->GetY1NDC();
 g=TPave->GetX2NDC();
 h=TPave->GetY2NDC();
 printf("%.1f,%.1f,%.1f,%.1f,\"br\"\n",a,b,c,d);
 printf("%.3f,%.3f,%.3f,%.3f,\"NDC\"\n",e,f,g,h);
 FILE * savefile;
 savefile=fopen("pave.lst","w");
 fprintf(savefile,"%.1f,%.1f,%.1f,%.1f,\"br\"\n",a,b,c,d);
 fprintf(savefile,"%.3f,%.3f,%.3f,%.3f,\"NDC\"\n",e,f,g,h);
 fclose (savefile);
}

/* 4). Spectrum Utilities------------------------
 * 1-D peak search          : peakfit(), peakfitx(), peakfity() 
 * 1-D background estimation: bkgfit2()
 * 1-D smoothing            : smooth()
 */
void peakfit(Char_t *histin, Char_t *filename, Float_t resolution=2, Double_t sigma=3, Double_t threshold=0.05, Char_t *option="")
{//Program by AHW.  Modified to run in fit.cc and in "modern" version of ROOT.
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();//added
  hname=histin;
  for(Int_t i=0;i<hname.Length();i++){//loop added by Jack
    TString temp="";
    temp=hname(i,hname.Length()-i);
    if(temp.IsFloat())
      {
	det=temp.Atoi();
	break;
      }
  }

  cFit->Clear();
  cFit->Divide(1,2);
  Float_t *positions;//moved * before variable name
  Float_t energies[100];
  Float_t slope,offset,width;
  Float_t ein;
  Float_t max=0,min=0;//added
  Int_t npeaks,nlist=0;
  TSpectrum *spectrum=new TSpectrum();
  ifstream listfile(filename);
  while(listfile>>ein) {
    energies[nlist]=ein;
    if(nlist==0)min=energies[nlist];//added
    if(energies[nlist]>max)max=energies[nlist];//added
    if(energies[nlist]<min)min=energies[nlist];//added
    cout << "Energy "<<nlist<< "= "<<energies[nlist]<<endl;
    nlist++;
  }
if(nlist<=1){
    printf("WARNING: Only %d energy found in file \"%s\"\n         Check file to ensure there is a carriage return on the last line!\n",nlist,filename);
  }
  TH1F *hProj=(TH1F *) gROOT->FindObject(histin);//hInput changed to hProj from here on.
  if(gROOT->FindObject("hPeakFit"))hPeakFit->Delete();//added
  hfit=new TH1F("hPeakFit","hPeakFit",1024,min-(max-min)/4,max+(max-min)/4);//added
  if((min-(max-min)/4)<0)printf("Notice: \"%s\" contains negatives value(s).\n        All zero-content bins are shown.\n",hfit->GetTitle());//added
  //  TH1F *hfit =(TH1F *) gROOT->FindObject("hPeakFit");//needed?
  hfit->Reset();
  cFit->cd(1);
  spectrum->SetResolution(resolution);
  spectrum->Search(hProj,sigma,option,threshold);
  positions=spectrum->GetPositionX();//in ROOT 5.26 this array is ordered by peak height!
  npeaks=spectrum->GetNPeaks();
  cout << "Found "<<npeaks<<" peaks in spectrum."<<endl;
  if(gROOT->GetVersionInt()<52400)
    printf("Deprecated (online) ROOT version.\n");
  else{
    printf("Modern (offline) ROOT version.\n");
    Float_t sorted[100];
    sorted[0]=positions[0];//initializes sorted array with a valid position
    for(Int_t i=0;i<npeaks;i++){
      if((positions[i])<(sorted[0])){
	sorted[0]=positions[i];//locates minimum
      }
    }
    for(Int_t i=1;i<npeaks;i++){
      sorted[i]=hProj->GetXaxis()->GetXmax();//initializes array element with non-minimum
      for(Int_t j=0;j<npeaks;j++){
	if(((positions[j])<(sorted[i]))&&((positions[j])>(sorted[i-1])))
	  sorted[i]=positions[j];//locates next-smallest position
      }
    }
        positions=sorted;  
  }//end version IF 
  for (Int_t i=0; i<npeaks; i++){
    cout<<"Peak " <<i<<" found at channel "<<positions[i]<<endl;
    //  printf("                        %f\n",sorted[i]);
  }
  hProj->Draw();
  for (Int_t i=0; i<npeaks; i++) {
    hfit->Fill(energies[i],positions[i]);
  }
  hfit->SetAxisRange(min,max);//added to set fit range to peak range
  hfit->Fit("pol1","Q");
  hfit->GetXaxis()->UnZoom();
  hfit->SetStats(kFALSE);//
  slope=hfit->GetFunction("pol1")->GetParameter(1);
  offset=hfit->GetFunction("pol1")->GetParameter(0);
  hProj->Fit("gaus","QW","",positions[npeaks-1]-25,positions[npeaks-1]+25);
  width=hProj->GetFunction("gaus")->GetParameter(2);
  cout<<"Fit parameters are:  Slope= "<<slope<<" offset= "<<offset<<" sigma(peak "<<npeaks-1<<")="<<width<<endl;
  printf("Fit parameters are: Slope = %3.3f, Offset = %3.3f\n",slope,offset);
  hfit->SetMarkerStyle(2);
  hfit->SetMarkerColor(2);
  hfit->SetMarkerSize(3);
  cFit->cd(2);
  hfit->Draw("P");
  delete spectrum;
FILE * outfile;
 outfile=fopen("temp.lst","w");
 fprintf(outfile,"%g, %g\n",slope,offset);
 fclose(outfile);
}

void peakfitx(Char_t *histin, Char_t *filename, Float_t resolution=2, Double_t sigma=3, Double_t threshold=0.05, Char_t *option="")
{//extension of peakfit() - takes a 2D histogram as input
 hname=histin;
  for(Int_t i=0;i<hname.Length();i++){//loop added by Jack
    TString temp="";
    temp=hname(i,hname.Length()-i);
    if(temp.IsFloat())
      {
	det=temp.Atoi();
	break;
      }
  }
  Float_t * positions;
  Float_t energies[100];
  Float_t slope,offset,width;
  Float_t ein; 
  Float_t max=0,min=0;
  Float_t a=0,b=0;
  Int_t npeaks,nlist=0;
  TSpectrum *spectrum=new TSpectrum();
  ifstream listfile(filename);

  while(listfile>>ein) {
    energies[nlist]=ein;
    if(nlist==0)min=energies[nlist];
    if(energies[nlist]>max)max=energies[nlist];
    if(energies[nlist]<min)min=energies[nlist];
    cout << "Energy "<<nlist<< "= "<<energies[nlist]<<endl;
    nlist++;
  }
  if(nlist<=1){
    printf("WARNING: Only %d energy found in file \"%s\"\n         Check file to ensure there is a carriage return on the last line!\n",nlist,filename);
  }
  
  if(gROOT->FindObject("hPeakFit"))hPeakFit->Delete();//added
  hfit=new TH1F("hPeakFit","hPeakFit",1024,min-(max-min)/4,max+(max-min)/4);//added
  if((min-(max-min)/4)<0)printf("Notice: \"%s\" contains negatives value(s).\n        All zero-content bins are shown.\n",hfit->GetTitle());//added
  //  TH1F *hfit =(TH1F *) gROOT->FindObject("hPeakFit");  //needed?
  hfit->Reset();

  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  cFit->Clear();
  cFit->Divide(1,3);

  cFit->cd(1);
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin); 
  hInput->Draw("COL");
 
  cFit->cd(2);
  hname=histin;
  hname+="_px";
  hInput->ProjectionX(hname);
  hProj=(TH1F *) gROOT->FindObject(hname.Data());
  hProj->Draw(); 
  a=hProj->GetXaxis()->GetXmin();
  b=hProj->GetXaxis()->GetXmax();

  spectrum->SetResolution(resolution);
  spectrum->Search(hProj,sigma,option,threshold);
  positions=spectrum->GetPositionX();
  for(Int_t i=0;i<15;i++)
  npeaks=spectrum->GetNPeaks();
  cout << "Found "<<npeaks<<" peaks in spectrum."<<endl;
  if(gROOT->GetVersionInt()<52400)
    printf("Deprecated (online) ROOT version.\n");
  else{
    printf("Modern (offline) ROOT version.\n");
    Float_t sorted[10];
    sorted[0]=positions[0];//initializes sorted array with a valid position
    for(Int_t i=0;i<npeaks;i++){
      if((positions[i])<(sorted[0])){
	sorted[0]=positions[i];//locates minimum
      }
    }
    for(Int_t i=1;i<npeaks;i++){
      sorted[i]=hProj->GetXaxis()->GetXmax();//initializes array element with non-minimum
      for(Int_t j=0;j<npeaks;j++){
	if(((positions[j])<(sorted[i]))&&((positions[j])>(sorted[i-1])))
	  sorted[i]=positions[j];//locates next-smallest position
      }
    }
    positions=sorted;  
  }//end version IF 
  for (Int_t i=0; i<npeaks; i++){
      cout<<"Peak " <<i<<" found at channel "<<positions[i]<<endl;
  }
  for (Int_t i=0; i<npeaks; i++) {
    hfit->Fill(energies[i],positions[i]);
  }
  hfit->SetAxisRange(min-(max-min)/2,max+(max-min)/2);
  hfit->Fit("pol1","Q");
  slope=hfit->GetFunction("pol1")->GetParameter(1);
  offset=hfit->GetFunction("pol1")->GetParameter(0);
  hProj->Fit("gaus","Q","",positions[npeaks-1]-(b-a)/15,positions[npeaks-1]+(b-a)/15);
  width=hProj->GetFunction("gaus")->GetParameter(2);
  cout<<"Fit parameters are:  Slope= "<<slope<<" offset= "<<offset<<" sigma(peak "<<npeaks-1<<")="<<width<<endl;
  printf("Fit parameters are: Slope = %3.3f, Offset = %3.3f\n",slope,offset);
  hfit->SetMarkerStyle(2);
  hfit->SetMarkerColor(2);
  hfit->SetMarkerSize(3);
  cFit->cd(3);
  hfit->Draw("P");
  delete spectrum;
FILE * outfile;
 outfile=fopen("temp.lst","w");
 fprintf(outfile,"%g, %g\n",slope,offset);
 fclose(outfile);
}

void peakfity(Char_t *histin, Char_t *filename, Float_t resolution=2, Double_t sigma=3, Double_t threshold=0.05, Char_t *option="")
{//extension of peakfit() - takes a 2D histogram as input
hname=histin;
  for(Int_t i=0;i<hname.Length();i++){//loop added by Jack
    TString temp="";
    temp=hname(i,hname.Length()-i);
    if(temp.IsFloat())
      {
	det=temp.Atoi();
	break;
      }
  }

  Float_t * positions;
  Float_t energies[10];
  Float_t slope,offset,width;
  Float_t ein;
  Float_t max=0,min=0;
  Float_t a=0,b=0;
  Int_t npeaks,nlist=0;
  TSpectrum *spectrum=new TSpectrum();

  ifstream listfile(filename);

  while(listfile>>ein) {
    energies[nlist]=ein;
    if(nlist==0)min=energies[nlist];
    if(energies[nlist]>max)max=energies[nlist];
    if(energies[nlist]<min)min=energies[nlist];
    cout << "Energy "<<nlist<< "= "<<energies[nlist]<<endl;
    nlist++;
  }
  if(nlist<=1){
    printf("WARNING: Only %d energy found in file \"%s\"\n         Check file to ensure there is a carriage return on the last line!\n",nlist,filename);
  }
  
  if(gROOT->FindObject("hPeakFit"))hPeakFit->Delete();//added
  hfit=new TH1F("hPeakFit","hPeakFit",1024,min-(max-min)/4,max+(max-min)/4);//added
  if((min-(max-min)/4)<0)printf("Notice: \"%s\" contains negatives value(s).\n        All zero-content bins are shown.\n",hfit->GetTitle());//added
  //  TH1F *hfit =(TH1F *) gROOT->FindObject("hPeakFit");
  //hfit =(TH1F *) gROOT->FindObject("hPeakFit");
  hfit->Reset();
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  cFit->Clear();
  cFit->Divide(1,3);
  
  cFit->cd(1);
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  hInput->Draw("COL");

  cFit->cd(2);
  hname=histin;
  hname+="_py"; 
  hInput->ProjectionY(hname);
  hProj=(TH1F *) gROOT->FindObject(hname.Data());
  hProj->Draw();
  a=hProj->GetXaxis()->GetXmin();
  b=hProj->GetXaxis()->GetXmax();
 
  spectrum->SetResolution(resolution);
  spectrum->Search(hProj,sigma,option,threshold);
  positions=spectrum->GetPositionX();
  npeaks=spectrum->GetNPeaks();
  cout << "Found "<<npeaks<<" peaks in spectrum."<<endl;
  if(gROOT->GetVersionInt()<52400)
    printf("Deprecated (online) ROOT version.\n");
  else{
    printf("Modern (offline) ROOT version.\n");
    Float_t sorted[10];
    sorted[0]=positions[0];//initializes sorted array with a valid position
    for(Int_t i=0;i<npeaks;i++){
      if((positions[i])<(sorted[0])){
	sorted[0]=positions[i];//locates minimum
      }
    }
    for(Int_t i=1;i<npeaks;i++){
      sorted[i]=hProj->GetXaxis()->GetXmax();//initializes array element with non-minimum
      for(Int_t j=0;j<npeaks;j++){
	if(((positions[j])<(sorted[i]))&&((positions[j])>(sorted[i-1])))
	  sorted[i]=positions[j];//locates next-smallest position
      }
    }
    positions=sorted;  
  }//end version IF  
  for (Int_t i=0; i<npeaks; i++){
      cout<<"Peak " <<i<<" found at channel "<<positions[i]<<endl;
    }
  for (Int_t i=0; i<npeaks; i++) {
    hfit->Fill(energies[i],positions[i]);
  }
  hProj->Fit("gaus","","",positions[npeaks-1]-25,positions[npeaks-1]+25);
  width=hProj->GetFunction("gaus")->GetParameter(2);

 cFit->cd(3);
 hfit->SetAxisRange(min-(max-min)/2,max+(max-min)/2);
 hfit->SetMarkerStyle(2);
 hfit->SetMarkerColor(2);
 hfit->SetMarkerSize(3);
 hfit->Fit("pol1","","P");
 slope=hfit->GetFunction("pol1")->GetParameter(1);
 offset=hfit->GetFunction("pol1")->GetParameter(0);
 printf("Fit parameters are: Slope = %3.3f, Offset = %3.3f\n",slope,offset);
 // cout<<"Fit parameters are:  Slope= "<<slope<<" offset= "<<offset<<" sigma(peak "<<npeaks-1<<")="<<width<<endl;
 FILE * outfile;
 outfile=fopen("temp.lst","w");
 fprintf(outfile,"%g, %g\n",slope,offset);
 fclose(outfile);

}

void bkgfit2(Char_t *histin, Float_t resolution=2, Double_t sigma=3, Double_t threshold=.05, Int_t niter=20, Char_t *option="")
{//modified to run in fit.cc and automatically copy histin.
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();//added
  if(!((TH1F *) gROOT->FindObject("hPeakFit"))) hfit=new TH1F("hPeakFit","hPeakFit",1000,0,100);  //added
  hfit=(TH1F *) gROOT->FindObject("hPeakFit");//added
 
  cFit->Clear();//cFitCanvas renamed to CFit
  cFit->Divide(1,2);
  
  Float_t *positions;//* moved to start to avoid automatic variable warning
  Float_t energies[10];
  Float_t slope,offset,width;
  Float_t ein;
  Int_t npeaks,nlist=0;
  float result[2048];
  Int_t nbins;
  TSpectrum *spectrum=new TSpectrum();
 
  //The non-histogram-declaration lines in the following block were added to copy the input histogram automatically.
  TH1F * hProj=(TH1F *) gROOT->FindObject(histin);
  hname=histin;
  hname+="_bkg"; 
  hProj->Clone(hname);
  TString histbkg=hname.Data();
  TH1F *hBkg=(TH1F *) gROOT->FindObject(hname.Data());
  hname=histin;
  hname+="_out"; 
  hProj->Clone(hname);
  TString histout=hname.Data();
  TH1F *hRresult=(TH1F *) gROOT->FindObject(hname.Data());

  hfit->Reset();
  cFit->cd(1);
  spectrum->SetResolution(resolution);
  spectrum->Search(hProj,sigma,option,threshold);
  positions=spectrum->GetPositionX();
  npeaks=spectrum->GetNPeaks();
  cout << "Found "<<npeaks<<" peaks in spectrum."<<endl;
  for (Int_t i=0; i<npeaks; i++){
      cout<<"Peak " <<i<<" found at channel "<<positions[i]<<endl;
    }
  cFit->cd(1);
  hProj->Draw();
  nbins=hProj->GetNbinsX();
  for (Int_t i=0; i<nbins; i++){
      result[i]=hProj->GetBinContent(i);
    }
  spectrum->Background(result,nbins,niter,1,1,kFALSE,1,kFALSE);
  for (Int_t i=0; i<nbins; i++){
      hBkg->SetBinContent(i,result[i]);
    }
  hBkg->SetLineColor(2);
  hBkg->Draw("same");
  cFit->cd(2);
  subtract(histin,histbkg,histout);
 }

void smooth(Char_t *histin,Int_t averWindow=3)
{
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
  TH1F * hProj=(TH1F *) gROOT->FindObject(histin);
  hname=histin;
  hname+="_smooth";
  hProj->Clone(hname);
  TH1F *hResult=(TH1F *) gROOT->FindObject(hname.Data());

  Float_t nbins;
  nbins=hProj->GetXaxis()->GetNbins();
 
  Float_t * source = new float[nbins];
  for(int i=0;i<nbins;i++){
    source[i]=hProj->GetBinContent(i+1);
  }
      
  TSpectrum *s = new TSpectrum();
  s->SmoothMarkov(source,nbins,averWindow);  //3, 7, 10
  for (int i = 0; i < nbins; i++) hResult->SetBinContent(i + 1,source[i]);
  cFit->Clear();//cFitCanvas renamed to CFit
  cFit->Divide(1,2);
  cFit->cd(1);
  hProj->Draw();
  cFit->cd(2);
  hResult->Draw("L");
 }

/* 5). File Utilities---------------------------------
 */

void printcaldata(void)
{
  Float_t p0av=0;
  Int_t entries=0;
  printf("Contents of array are: \n");
  for(int i=0;i<24;i++){
      for(int j=0;j<11;j++){
	  //cout<<array[i][j]<<" ";
	  if(j==0)
	    printf("%2d ",array[i][j]);
	  else
	    printf("%7.0f ",array[i][j]);
	}
      cout<<endl;
      if(array[i][0]){
	  // printf("Det %2d, p0 is %f, p0av is %f\n",i+1,array[i][1],p0av);
	  p0av+=array[i][1];
	  entries++;
	}
    }
  p0av=p0av/entries;
  printf("p0 average is %5.1f for %d entries\n",p0av,entries);
}

void createfile(Int_t numbered=0)
{//writes array[i][j] to file or creates blank file
  ofstream outfile("calibration.cal");
  for(int i=0;i<24;i++){
    outfile<<i+1<<" ";
    for(int j=0;j<11;j++){
      switch(numbered){
      case 0:
	outfile<<array[i][j]<<" ";
	break;
      case 1:
	outfile<<0<<" ";
	break;
      case 2:
	outfile<<j<<" ";
	break;
      default:break;  
      }
    }
    outfile<<endl;
  }
}

void readfile(Char_t *filename="test.cal", Bool_t showtest=1)
{
  Float_t param[24][50]; 
  Int_t errorline=-1;
  Int_t size=sizeof(param[0])/sizeof(param[0][0]);
  printf("param array size is: [24][%d].\n",size); 
  Bool_t fit=kFALSE;
  for(Int_t i=0;i<24;i++)
    for(Int_t j=0;j<size;j++)
      param[i][j]=0;//initializes all array elements to zero
 
  FILE * infile;
  Int_t k=1;
  while(!fit&&k<=(size)){
      infile = fopen (filename,"r");
      for(Int_t i=0;i<24;i++){
	for(Int_t j=0;j<k;j++){
	  fscanf(infile,"%f",&param[i][j]);
	}
      }
      fclose(infile);
      if(showtest) printf("Testing array length %d:\n",k);

      if (param[0][0]==1)fit=kTRUE;
      else fit=kFALSE;
      for(Int_t i=0;i<24;i++){
	fit=(fit&&(param[i][0]==(i+1)));	
	if(fit){
	  if(showtest) printf("%2.0f ",param[i][0]);
	  for(Int_t j=1;j<k;j++){
	    if(showtest) printf("%7.2f ",param[i][j]);
	    }
	 if(showtest) printf("\n");
	  if((i+1)>errorline)errorline=i+1;
	}
      }
           
     if(fit)printf("File \"%s\" has %d elements per line.\n",filename,k);
     
    k++;
}
  if(fit){   
  FILE * outfile;
  outfile=fopen("output.txt","w");
  printf("The contents of \"%s\" are:\n",filename);
  for(Int_t i=0;i<24;i++){
    fprintf(outfile,"%2.0f ",param[i][0]);
    printf("%2.0f ",param[i][0]);
    for(Int_t j=1;j<k-1;j++){
      fprintf(outfile,"%g ",param[i][j]);
      printf("%7.2f ",param[i][j]);
    }
    fprintf(outfile,"\n");
    printf("\n");
  }
  fclose(outfile);
  }
  else
    printf("File \"%s\" has more than %d elements per line, or there is an error on line %d.\n",filename,size,errorline);
}

void truncfile(Char_t *filename="source.cal",Char_t *output="shortsource.cal",Int_t trunc=6)
{
  Float_t param[24][50]; 
  Int_t errorline=-1;
  Int_t size=sizeof(param[0])/sizeof(param[0][0]);
  printf("param array size is: [24][%d].\n",size); 
  Bool_t fit=kFALSE;
  for(Int_t i=0;i<24;i++)
    for(Int_t j=0;j<size;j++)
      param[i][j]=0;//initializes all array elements to zero
 
  FILE * infile;
  Int_t k=2;
  while(!fit&&k<=(size)){
      infile = fopen (filename,"r");
      for(Int_t i=0;i<24;i++){
	for(Int_t j=0;j<k;j++){
	  fscanf(infile,"%f",&param[i][j]);
	}
      }
      fclose(infile);
      printf("Testing array length %2d:\n",k);
      if (param[0][0]==1)fit=kTRUE;
      else fit=kFALSE;
      for(Int_t i=0;i<24;i++){
	fit=(fit&&(param[i][0]==(i+1)));	
	if(fit){
	  printf("%2.0f ",param[i][0]);
	  for(Int_t j=1;j<k;j++){
	    printf("%7.3f ",param[i][j]);
	  }
	  printf("\n");
	  if((i+1)>errorline)errorline=i+1;
	}
      }
      if(fit)printf("File \"%s\" has %d elements per line.\n",filename,k);
      k++;
    }
  if(fit){   
    FILE * outfile;
    outfile=fopen(output,"w");
    printf("The new contents of \"%s\" are:\n",output);
    for(Int_t i=0;i<24;i++){
      fprintf(outfile,"%2.0f ",param[i][0]);
      printf("%2.0f ",param[i][0]);
      for(Int_t j=1;j<trunc;j++){
	fprintf(outfile,"%g ",param[i][j]);
	printf("%7.3f ",param[i][j]);
      }
      fprintf(outfile,"\n");
      printf("\n");
    }
    fclose(outfile);
  }
  else
    printf("File \"%s\" has more than %d elements per line, or there is an error on line %d.\n
File \"%s\" not written.\n",filename,size,errorline,output);
}

void writetemp(Char_t *calfile="calibration.cal", Int_t detno=-1, Char_t *tempfile="temp.lst",Int_t index=0,Bool_t overwrite=1,Bool_t showtest=0)
{
  if(detno==-1)
    detno=det;
  Float_t param[24][50]; 
  Int_t errorline=-1;
  Int_t size=sizeof(param[0])/sizeof(param[0][0]);
  if(showtest)  printf("param array size is: [24][%d].\n",size); 
  Bool_t fit=kFALSE;
  for(Int_t i=0;i<24;i++)
    for(Int_t j=0;j<size;j++)
      param[i][j]=0;//initializes all array elements to zero
  
  FILE * infile;
  Int_t k=1;
  
  while(!fit&&k<=(size)){
    infile = fopen (calfile,"r");
    if(infile!=NULL){
      for(Int_t i=0;i<24;i++){
	for(Int_t j=0;j<k;j++){
	  fscanf(infile,"%f",&param[i][j]);
	}
      }
      fclose(infile);
    }
    else{
      printf("File is NULL!\n");
      //for(Int_t l=0;l<size;l++){
      //	Fit[l]=0;
    }
      
      if(showtest) printf("Testing array length %d:\n",k);
      
      if(param[0][0]==1)
	fit=kTRUE;
      else
	fit=kFALSE;
      
      for(Int_t i=0;i<24;i++){
	fit=(fit&&(param[i][0]==(i+1)));	
	if(fit){
	  if(showtest) printf("%2.0f ",param[i][0]);
	  for(Int_t j=1;j<k;j++){
	    if(showtest) printf("%7.2f ",param[i][j]);
	    }
	 if(showtest) printf("\n");
	  if((i+1)>errorline)errorline=i+1;
	}
      }
           
     if(fit)printf("File \"%s\" has %d elements per line.\n",calfile,k);
     
    k++;
  }
    
  if(fit){   
    Float_t Fit[10];
    Int_t size=sizeof(Fit)/sizeof(Fit[0]); 
    infile = fopen (tempfile,"r");
    if(infile!=NULL){
      for(Int_t l=0;l<size;l++){
	if(fscanf(infile,"%f, ",&Fit[l])==-1)
	  Fit[l]=0;
	if(Fit[l])//added
	  param[detno-1][l+1+index]=Fit[l];//note use of index
      }
      fclose(infile);
    }
    else
      printf("File is NULL!\n");
     
    FILE * outfile;
    if(overwrite)
      outfile=fopen(calfile,"w");
    else{
      outfile=fopen("output.txt","w");
      printf("Output file is \"output.txt\"\n");
    }
    printf("The contents of \"%s\" are approximately:\n",calfile);
    for(Int_t i=0;i<24;i++){
      fprintf(outfile,"%2.0f ",param[i][0]);
      printf("%2.0f ",param[i][0]);
    for(Int_t j=1;j<k-1;j++){
      fprintf(outfile,"%g ",param[i][j]);
      printf("%5.0f ",param[i][j]);
    }
    fprintf(outfile,"\n");
    printf("\n");
  }
  fclose(outfile);
  }
  else
    printf("File \"%s\" has more than %d elements per line, or there is an error on line %d.\n",calfile,size,errorline);
  
}

