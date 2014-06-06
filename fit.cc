/*
Line fitting program.  Fits a 2-D histogram with a line.
*/
#include <iostream.h>
TH1F *hProj=0;
TH1F *hProf=0;
TH2F *hInput=0;
TH2F *hOutput=0;
TH1F *hfit=0;
Float_t array[24][11];
Float_t slope=0;
Float_t offset=0;
Int_t det=0;
TString hname;
//TFile *_file0;
//TFile *_file1;
//TFile *_file2 = new TFile("100_cut.root");

void plotsim(Char_t* hname="test.dat"){
  filltree(hname.Data(),"tree");
  analyze("tree");
}

void mkCanvas2(Char_t* cvname="cFit",Char_t *cvtitle="cFit",Int_t ww=640,Int_t wh=646){
  TCanvas * cFit=new TCanvas(cvname,cvtitle,0,0,ww,wh);
  if(!(cFit->GetShowEventStatus()))cFit->ToggleEventStatus();
  if(!(cFit->GetShowToolBar()))cFit->ToggleToolBar();
}

void plotall(Char_t *histin,Float_t minX=0,Float_t maxX=0,Float_t minY=0,Float_t maxY=0,Int_t scale=1){
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2("cFit","cFit",1272,695);
  Int_t no=0;
  for(int i=1;i<25;++i){
    hname=histin;
    hname=hname+i;
    if((TH2F*)gROOT->FindObject(hname.Data()))
      no++;
  }
  printf("Plotting %2d Histograms...\n",no);
  cFit->Clear();
  cFit->Divide(6,(no/6));
  for(int i=1;i<(no+1);++i){
    cFit->cd(i);
    hname=histin;
    hname=hname+i;
    
    hInput=(TH2F*)gROOT->FindObject(hname.Data());
     
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
    hInput->Draw();
  }
}

void plotallpjx(Char_t *histin,Float_t minX=0,Float_t maxX=0,Float_t minY=0,Float_t maxY=0,Int_t scale=1){
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2("cFit","cFit",1272,695);
  Int_t no=0;
  for(int i=1;i<25;++i){
    hname=histin;
    hname=hname+i;
    if((TH2F*)gROOT->FindObject(hname.Data()))
      no++;
  }
  printf("Plotting %2d Histograms...\n",no);
  cFit->Clear();
  cFit->Divide(6,(no/6));
  for(int i=1;i<(no+1);++i){
    cFit->cd(i);
    hname=histin;
    hname=hname+i;
    
    hInput=(TH2F*)gROOT->FindObject(hname.Data());
    /*    
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
    */
    hInput->ProjectionX();
    hname=hname+"_px";
    hProj=(TH1F *) gROOT->FindObject(hname.Data());
    hProj->Draw();
  }
}

void plotallpjy(Char_t *histin,Float_t minX=0,Float_t maxX=0,Float_t minY=0,Float_t maxY=0,Int_t scale=1){
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2("cFit","cFit",1272,695);
  Int_t no=0;
  for(int i=1;i<25;++i){
    hname=histin;
    hname=hname+i;
    if((TH2F*)gROOT->FindObject(hname.Data()))
      no++;
  }
  printf("Plotting %2d Histograms...\n",no);
  cFit->Clear();
  cFit->Divide(6,(no/6));
  for(int i=1;i<(no+1);++i){
    cFit->cd(i);
    hname=histin;
    hname=hname+i;
    
    hInput=(TH2F*)gROOT->FindObject(hname.Data());
    /*    
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
    */
    hInput->ProjectionY();
    hname=hname+"_py";
    hProj=(TH1F *) gROOT->FindObject(hname.Data());
    hProj->Draw();
  }
}

void setrange(Char_t *histin,Float_t minX=0,Float_t maxX=0,Float_t minY=0,Float_t maxY=0,Int_t scale=1){
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

void opjy(Char_t *histin,Int_t col=2)
{
  hInput=(TH2F *) gROOT->FindObject(histin);
  hInput->ProjectionY();
  hname=histin;
  hname+="_py";

  hProj=(TH1F *) gROOT->FindObject(hname.Data());
  hProj->SetLineColor(col);
  hProj->Draw("same");
}

void opjx(Char_t *histin,Int_t col=2)
{
  hInput=(TH2F *) gROOT->FindObject(histin);
  hInput->ProjectionX();
  hname=histin;
  hname+="_px";

  hProj=(TH1F *) gROOT->FindObject(hname.Data());
  hProj->SetLineColor(col);
  hProj->Draw("same");
}



void fitpfx(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Float_t minfit=0,Float_t maxfit=0,Int_t scale=1)
{
  Float_t a=0,b=0;  
  Float_t p0=0,p1=0,p2=0;
 
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();  
  cFit->Clear();
  cFit->SetWindowPosition(0,0);
  cFit->Divide(1,2);
  cFit->cd(1);

  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
  printf("input histogram is %s",histin);    
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
    hInput->GetYaxis()->UnZoom();
    hInput->GetYaxis()->UnZoom();

  }

  hInput->Draw("COL2");
  cFit->cd(2); 
  a=minpf;
  b=maxpf; 
  printf("Projection Limits are %3.2f to %3.2f\n",minpf,maxpf);  
  minpf=hInput->GetYaxis()->FindBin(minpf);
  maxpf=hInput->GetYaxis()->FindBin(maxpf);
  printf("Fit Limits are %3.2f to %3.2f\n",minfit,maxfit);
 
  hname=histin;
  hname+="_pfx"; 
  
  hInput->ProfileX(hname,minpf,maxpf);
  hProf=(TH1F *) gROOT->FindObject(hname.Data());
  cFit->cd(2);
  hProf->SetAxisRange(a,b,"Y");
  hProf->Draw();
 
  hProf->Fit("pol1","V","",minfit,maxfit);
  p0=hProf->GetFunction("pol1")->GetParameter(0);
  p1=hProf->GetFunction("pol1")->GetParameter(1);
 
  printf("p1 = %7.3f\n",p1);

}

void fit2pfx(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Float_t minfit=0,Float_t maxfit=0,Int_t scale=1)
{
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2(); 
Float_t cp=0 ;
   Float_t p0=0,p1=0,p2=0;
   Float_t a=0,b=0;  
  cFit->Clear();
 cFit->SetWindowPosition(0,0); 
 cFit->Divide(1,2);
  cFit->cd(1);
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
    
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
     a=minpf;
  b=maxpf;  
  printf("Projection Limits are %d to %d\n",minpf,maxpf);  
  minpf=hInput->GetYaxis()->FindBin(minpf);
  maxpf=hInput->GetYaxis()->FindBin(maxpf);
  printf("Fit Limits are %d to %d\n",minfit,maxfit);
 
  hname=histin;
  hname+="_pfx"; 
  
  hInput->ProfileX(hname,minpf,maxpf);
  hProf=(TH1F *) gROOT->FindObject(hname.Data());
  cFit->cd(2);
   hProf->SetAxisRange(a,b,"Y");
  hProf->Draw();
 
  hProf->Fit("pol2","V","",minfit,maxfit);
  p0=hProf->GetFunction("pol2")->GetParameter(0);
  p1=hProf->GetFunction("pol2")->GetParameter(1);
  p2=hProf->GetFunction("pol2")->GetParameter(2);
  cp=(-p1/(2*p2));
  printf("p1 = %7.3f p2=%7.3f\ncp = %7.3f\n",p1,p2,cp);
}
void fit4pfx(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Float_t minfit=0,Float_t maxfit=0,Int_t scale=1)
{
   if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
Float_t cp=0 ;
   Float_t p0=0,p1=0,p2=0,p3=0,p4=0;
   Float_t a=0,b=0;  
  cFit->Clear();
  cFit->SetWindowPosition(0,0);
  cFit->Divide(1,2);
  cFit->cd(1);
  TH2F * hInput=(TH2F *) gROOT->FindObject(histin);
    
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
  a=minpf;
  b=maxpf;  
  printf("Projection Limits are %d to %d\n",minpf,maxpf);  
  minpf=hInput->GetYaxis()->FindBin(minpf);
  maxpf=hInput->GetYaxis()->FindBin(maxpf);
  printf("Fit Limits are %d to %d\n",minfit,maxfit);
 
  hname=histin;
  hname+="_pfx"; 
  
  hInput->ProfileX(hname,minpf,maxpf);
  hProf=(TH1F *) gROOT->FindObject(hname.Data());
  cFit->cd(2);
   hProf->SetAxisRange(a,b,"Y");
  hProf->Draw();
 
  hProf->Fit("pol4","V","",minfit,maxfit);
  p0=hProf->GetFunction("pol4")->GetParameter(0);
  p1=hProf->GetFunction("pol4")->GetParameter(1);
  p2=hProf->GetFunction("pol4")->GetParameter(2);
  p3=hProf->GetFunction("pol4")->GetParameter(3);
  p4=hProf->GetFunction("pol4")->GetParameter(4);
  printf("p1 = %5.0f, p2 = %5.0f, p3 = %5.0f, p4 = %5.0f\ncp = %7.3f\n",p1,p2,p3,p4);
}

void fitpfy(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Float_t minfit=0,Float_t maxfit=0,Int_t scale=1)
{
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
Float_t p0=0,p1=0,p2=0;
  
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
  hProf->SetAxisRange(a,b,"Y");
  hProf->Draw();
 
  hProf->Fit("pol1","V","",minfit,maxfit);
  p0=hProf->GetFunction("pol1")->GetParameter(0);
  p1=hProf->GetFunction("pol1")->GetParameter(1);
  //  p2=hProf->GetFunction("pol2")->GetParameter(2);
  // cp=(-p1/(2*p2));
  printf("p1 = %7.3f\n",p1);
}
void fitall(void)
{
 
 //this function fills the array with the calibration data.
  for(int i=0;i<24;i++)
    {
      array[i][0]=i+1;//just the detector number
      hname="hESum";
      hname+=(i+1);
      cout<<"hname = "<<hname<<endl;
      det=array[i][0];array[i][3]=slope;array[i][4]=offset;
    }
  linefit("hESum1");det=array[0][0];array[0][3]=slope;array[0][4]=offset;
  linefit("hESum2");det=array[1][0];array[1][3]=slope;array[1][4]=offset;
  linefit("hESum3");det=array[2][0];array[2][3]=slope;array[2][4]=offset;
  linefit("hESum4");det=array[3][0];array[3][3]=slope;array[3][4]=offset;
  linefit("hESum5");det=array[4][0];array[4][3]=slope;array[4][4]=offset;
  linefit("hESum6");det=array[5][0];array[5][3]=slope;array[5][4]=offset;
  linefit("hESum7");det=array[6][0];array[6][3]=slope;array[6][4]=offset;
  linefit("hESum8");det=array[7][0];array[7][3]=slope;array[7][4]=offset;
  linefit("hESum9");det=array[8][0];array[8][3]=slope;array[8][4]=offset;
  linefit("hESum10");det=array[9][0];array[9][3]=slope;array[9][4]=offset;
  linefit("hESum11");det=array[10][0];array[10][3]=slope;array[10][4]=offset;
  linefit("hESum12");det=array[11][0];array[11][3]=slope;array[11][4]=offset;
  linefit("hESum13");det=array[12][0];array[12][3]=slope;array[12][4]=offset;
  linefit("hESum14");det=array[13][0];array[13][3]=slope;array[13][4]=offset;
  linefit("hESum15");det=array[14][0];array[14][3]=slope;array[14][4]=offset;
  linefit("hESum16");det=array[15][0];array[15][3]=slope;array[15][4]=offset;
  linefit("hESum17");det=array[16][0];array[16][3]=slope;array[16][4]=offset;
  linefit("hESum18");det=array[17][0];array[17][3]=slope;array[17][4]=offset;
  linefit("hESum19");det=array[18][0];array[18][3]=slope;array[18][4]=offset;
  linefit("hESum20");det=array[19][0];array[19][3]=slope;array[19][4]=offset;
  linefit("hESum21");det=array[20][0];array[20][3]=slope;array[20][4]=offset;
  linefit("hESum22");det=array[21][0];array[21][3]=slope;array[21][4]=offset;
  linefit("hESum23");det=array[22][0];array[22][3]=slope;array[22][4]=offset;
  linefit("hESum24");det=array[23][0];array[23][3]=slope;array[23][4]=offset;
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
	printf("Projectin range is %f to %f\n",minpj,maxpj);
	printf("Fit range is %f %f\n",minfit,maxfit);
      }


for(int i=0;i<6;i++)
    {
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
      if(i==0)printf("Fit funcion is \"%s\"\n",fname.Data()); 
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

void createfile(void)
{
  ofstream outfile("calibration.cal");
  for(int i=0;i<24;i++)
    {
      for(int j=0;j<11;j++)
	{
	  outfile<<array[i][j]<<" ";
	}
      outfile<<endl;
    }
}

void printcaldata(void)
{
  Float_t p0av=0;
  Int_t entries=0;
  printf("Contents of array are: \n");
  for(int i=0;i<24;i++)
    {
      for(int j=0;j<11;j++)
	{
	  //cout<<array[i][j]<<" ";
	  if(j==0)
	    printf("%2d ",array[i][j]);
	  else
	    printf("%7.0f ",array[i][j]);
	}
      cout<<endl;
      if(array[i][0])
	{
	  // printf("Det %2d, p0 is %f, p0av is %f\n",i+1,array[i][1],p0av);
	  p0av+=array[i][1];
	  entries++;
	}
    }
  p0av=p0av/entries;
  printf("p0 average is %5.1f for %d entries\n",p0av,entries);
}

void minfit(Char_t *histin, Int_t minpf=0, Int_t maxpf=1200)
{
 if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();
 cFit->Clear();
  cFit->Divide(1,2);
  hInput=(TH2F *) gROOT->FindObject(histin);
  hInput->SetAxisRange(minpf,maxpf,"Y"); 
  /*
  hInput->ProjectionX();
  hname=histin;
  hname+="_px";
  hProj=(TH1F *) gROOT->FindObject(hname.Data());
  */
  cFit->cd(1);
  hInput->Draw("COL2");

  minpf=hInput->GetYaxis()->FindBin(minpf);
  maxpf=hInput->GetYaxis()->FindBin(maxpf);
  printf("Projection Limits are t = %d to %d",minpf,maxpf);
  hname=histin;
  hname+="_pfx";

  hInput->ProfileX(hname,minpf,maxpf);

  hProj=(TH1F *) gROOT->FindObject(hname.Data());

  cFit->cd(2);
  hProj->Draw();

  hProj->Fit("pol2","","",0,1);
  Float_t p2=hProj->GetFunction("pol2")->GetParameter(2);
  Float_t p1=hProj->GetFunction("pol2")->GetParameter(1);
  Float_t p0=hProj->GetFunction("pol2")->GetParameter(0);
}

void timefit(Char_t *histin,Float_t minE=1,Float_t maxE=12,Int_t minpf=0,Int_t maxpf=1200)
{
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

void fitpqpfy(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Float_t minfit=0,Float_t maxfit=0,Int_t scale=1)
{
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

      printf("Reiniating loop with tolerance set to %2.0f%%\n",tol*100);

     i=0;
    }
  }
    printf("Loop Exited at Iteration %3d.\nEndpoint tolerance is %2.0f%%\nFit range is %6.3f to %6.3f\nCritical Point is %5.3f\np0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",i,tol*100,minfit,maxfit,cp,p0,p1,p2);
    printf("Emax = %5.2f, p1 = %4.0f, p2 = %7.2f\n",cp,p1,p2);
}
void fitpqpfx(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Float_t minfit=0,Float_t maxfit=0,Int_t scale=1)
{
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

      printf("Reiniating loop with tolerance set to %2.0f%%\n",tol*100); cp=-1;

      i=0;
     }
  }
printf("
Loop Exited at Iteration %3d.\nEndpoint tolerance is %2.0f%%\nFit range is %6.3f to %6.3f\nCritical Point is %5.3f\np0 = %7.3f, p1 = %7.3f, p2 = %7.3f\n",i,tol*100,minfit,maxfit,cp,p0,p1,p2);
}
void fitpqRpfx(Char_t *histin,Float_t minpf=0,Float_t maxpf=0,Float_t minfit=0,Float_t maxfit=0,Int_t scale=1)
{
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
{
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

void getline(void){
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

void pjxy(Char_t *histin)
{
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

void peakfit(Char_t *histin, Char_t *filename, Float_t resolution=2, Double_t sigma=3, Double_t threshold=0.05, Char_t *option="")
{
  if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2(); 
  cFit->Clear();
  cFit->Divide(1,2);
  
  Float_t *positions;
  Float_t energies[10];
  Float_t slope,offset,width;
  Float_t ein;
  Float_t max=0,min=0;
  Float_t a=0,b=0;
  Int_t npeaks,nlist=0;
  TSpectrum *spectrum=new TSpectrum();
 if(!((TH1F *) gROOT->FindObject("hPeakFit"))) hfit=new TH1F("hPeakFit","hPeakFit",1000,0,100);
  ifstream listfile(filename);
  while(listfile>>ein) {
    energies[nlist]=ein;
    if(nlist==0)min=energies[nlist];
    if(energies[nlist]>max)max=energies[nlist];
    if(energies[nlist]<min)min=energies[nlist];
    cout << "Energy "<<nlist<< "= "<<energies[nlist]<<endl;
    nlist++;
  }
  TH1F *hProj=(TH1F *) gROOT->FindObject(histin);
  TH1F *hfit =(TH1F *) gROOT->FindObject("hPeakFit");
  hfit->Reset();
  cFit->cd(1);
  spectrum->SetResolution(resolution);
  spectrum->Search(hProj,sigma,option,threshold);
  positions=spectrum->GetPositionX();
  npeaks=spectrum->GetNPeaks();
  cout << "Found "<<npeaks<<" peaks in spectrum."<<endl;
  for (Int_t i=0; i<npeaks; i++)
    {
      cout<<"Peak " <<i<<" found at channel "<<positions[i]<<endl;
    }
  hProj->Draw();
  for (Int_t i=0; i<npeaks; i++) {
    hfit->Fill(energies[i],positions[i]);
  }
  hfit->SetAxisRange(min-(max-min)/2,max+(max-min)/2);
  hfit->Fit("pol1","Q");
  slope=hfit->GetFunction("pol1")->GetParameter(1);
  offset=hfit->GetFunction("pol1")->GetParameter(0);
  hProj->Fit("gaus","QW","",positions[npeaks-1]-25,positions[npeaks-1]+25);
  width=hProj->GetFunction("gaus")->GetParameter(2);
   cout<<"Fit parameters are:  Slope= "<<slope<<" offset= "<<offset<<" sigma(peak "<<npeaks-1<<")="<<width<<endl;
  printf("Fit parameters are: Slope = %3.3f, Ofset = %3.3f\n",slope,offset);
  hfit->SetMarkerStyle(2);
  hfit->SetMarkerColor(2);
  hfit->SetMarkerSize(3);
  cFit->cd(2);
  hfit->Draw("P");
  delete spectrum;
}

void peakfitx(Char_t *histin, Char_t *filename, Float_t resolution=2, Double_t sigma=3, Double_t threshold=0.05, Char_t *option="")
{
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
  if(!((TH1F *) gROOT->FindObject("hPeakFit"))) hfit=new TH1F("hPeakFit","hPeakFit",1000,0,100);  
  TH1F *hfit=(TH1F *) gROOT->FindObject("hPeakFit");
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
  for (Int_t i=0; i<npeaks; i++)
    {
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
  printf("Fit parameters are: Slope = %3.3f, Ofset = %3.3f\n",slope,offset);
  hfit->SetMarkerStyle(2);
  hfit->SetMarkerColor(2);
  hfit->SetMarkerSize(3);
  cFit->cd(3);
  hfit->Draw("P");
  delete spectrum;
}

void peakfity(Char_t *histin, Char_t *filename, Float_t resolution=2, Double_t sigma=3, Double_t threshold=0.05, Char_t *option="")
{
  Float_t * positions;
  Float_t energies[10];
  Float_t slope,offset,width;
  Float_t ein;
  Float_t max=0,min=0;
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

  TH1F *hfit=(TH1F *) gROOT->FindObject("hPeakFit");
  if(!((TH1F *) gROOT->FindObject("hPeakFit"))) hfit=new TH1F("hPeakFit","hPeakFit",1000,-100,100);
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
 
  spectrum->SetResolution(resolution);
  spectrum->Search(hProj,sigma,option,threshold);
  positions=spectrum->GetPositionX();
  npeaks=spectrum->GetNPeaks();
  cout << "Found "<<npeaks<<" peaks in spectrum."<<endl;
  for (Int_t i=0; i<npeaks; i++)
    {
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
 printf("Fit parameters are: Slope = %3.3f, Ofset = %3.3f\n",slope,offset);
 cout<<"Fit parameters are:  Slope= "<<slope<<" offset= "<<offset<<" sigma(peak "<<npeaks-1<<")="<<width<<endl;
}

void loadfile(Char_t *histin="")
{
  hname="100";
  hname+=histin; 
  hname+=".root"; 
  printf("Input File is %s\n",hname.Data());
  TFile *_file0 = new TFile(hname.Data());
}

void readfile(Char_t *filename="test.cal")
{
  Float_t param[24][50];
  Int_t size=sizeof(param[0])/sizeof(param[0][0]);
  Bool_t fit=kFALSE;
  for(Int_t i=0;i<24;i++)
    for(Int_t j=0;j<size;j++)
      param[i][j]=0;//initializes all array elements to zero
 
  FILE * infile;
  
  Int_t k=1;
  while(!fit&&k<=(size))
{
    infile = fopen (filename,"r");
    for(Int_t i=0;i<24;i++){
      for(Int_t j=0;j<k;j++){
        fscanf(infile,"%f",&param[i][j]);
      }
    }
  
    fclose(infile);
    for(Int_t i=0;i<24;i++){
            printf("%2.0f ",param[i][0]);
      if (param[0][0]==1)fit=kTRUE;
      else fit=kFALSE;
      fit=(fit&&(param[i][0]==(i+1)));
          for(Int_t j=1;j<k;j++){
              printf("%7.3f ",param[i][j]);
        }
            printf("\n");
	  
    }
    if(fit)printf("File has %d elements per line.\n",k);
     else printf("Array length %d does not fit\n",k);
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
      fprintf(outfile,"%7.3f ",param[i][j]);
      printf("%7.3f ",param[i][j]);
    }
    fprintf(outfile,"\n");
    printf("\n");
  }
  fclose(outfile);
  }
  else
    printf("File is larger than array!\n");
}

void merge(Char_t *histin)
{
  _file0->cd();
  hname=histin;
   TH2F * hInput=(TH2F *) gROOT->FindObject(hname.Data());
  hname+="_100"; 
  hInput->Clone(hname.Data());
  hOutput=(TH2F *) gROOT->FindObject(hname.Data());
  hOutput->SetDirectory(home);
  
  _file1->cd();
  hname=histin;
  TH2F *   hInput=(TH2F *) gROOT->FindObject(hname.Data());
  hname+="_350"; 
  hInput->Clone(hname.Data());
  hOutput=(TH2F *) gROOT->FindObject(hname.Data());
  hOutput->SetDirectory(home);
  
  _file2->cd();
  hname=histin;
  TH2F * hInput=(TH2F *) gROOT->FindObject(hname.Data());
  hname+="_500"; 
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
  add2(hname+"_100",hname+"_350",hname+"_all");
  add2(hname+"_500",hname+"_all",hname+"_all");
}

void loadfiles(Char_t *histin="")
{
  TDirectory *home=gDirectory;
  hname="100";
  hname+=histin; 
  hname+=".root"; 
  printf("Input File is %s\n",hname.Data());
  TFile *_file0 = new TFile(hname.Data());
  hname="350";
  hname+=histin; 
  hname+=".root"; 
  printf("Input File is %s\n",hname.Data());
  TFile *_file1 = new TFile(hname.Data());
  hname="500";
  hname+=histin; 
  hname+=".root"; 
  printf("Input File is %s\n",hname.Data());
  TFile *_file2 = new TFile(hname.Data());
}

//The following files have been copied from my personal modifications to helios_plottools.cc and util.cc

void doprint2(Char_t * cnvname="cFit", Char_t * filename="print.ps", Char_t * pr="f1-phaser")
{
  const char cmd[255];
  char * fl = reinterpret_cast<char *>(filename);
  TCanvas *thecanvas=(TCanvas *)gROOT->FindObject(cnvname);
  thecanvas->SaveAs(filename);
  sprintf(cmd,"lpr -P %s %s",pr,fl);
  gSystem->Exec(cmd);
  cout<<fl<<" was sent to printer "<<pr<<endl;
}

void add2(Char_t *histin1, Char_t *histin2, Char_t *histout, Float_t scale1=1.0, Float_t scale2=1.0)
{
  TH2F *hist1=(TH2F *) gROOT->FindObject(histin1);
  TH2F *hist2=(TH2F *) gROOT->FindObject(histin2);
  TH2F *hist3=(TH2F *) gROOT->FindObject(histout);
  hist3->Add(hist1,hist2,scale1,scale2);
  hist3->Draw("COL");
}

void copy2(Char_t *histin, Float_t minz=0)
{
if(!((TCanvas *) gROOT->FindObject("cFit"))) mkCanvas2();    

 Float_t xmin,xmax; 
 Float_t ymin,ymax; 
 Int_t xbin;
 Int_t ybin;
 Float_t x,y,z;
 TString hname;
  
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

 if ((TH2F *) gROOT->FindObject(hname)) gROOT->FindObject(hname)->Delete();  
 TH2F * hOutput=new  TH2F(hname,htitle,xbin,xmin,xmax,ybin,ymin,ymax);

 for(int i=0;i<xbin;i++){
   for(int j=0;j<ybin;j++){
     x=hInput->GetXaxis()->GetBinCenter(i);
     y=hInput->GetYaxis()->GetBinCenter(j);
     z=hInput->GetBinContent(i,j);
     if(z>minz)
       hOutput->Fill(x,y,z-minz);
   }
 }
 
 cFit->Clear();
 cFit->Divide(1,2);
 cFit->cd(1);
 hInput->Draw("col");
 cFit->cd(2);
 hOutput->Draw("col");
}
