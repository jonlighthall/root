//#include <iostream.h>
//#include <iomanip.h>
//#include <fstream.h>
//#include <stdio.h>
//#include <Type.h>

void add(Char_t *histin1, Char_t *histin2, Char_t *histout, Float_t scale1=1.0, Float_t scale2=1.0)
{
  TH1F *hist1=(TH1F *) gROOT->FindObject(histin1);
  TH1F *hist2=(TH1F *) gROOT->FindObject(histin2);
  TH1F *hist3=(TH1F *) gROOT->FindObject(histout);
  hist3->Add(hist1,hist2,scale1,scale2);
  hist3->Draw();
}

void add2(Char_t *histin1, Char_t *histin2, Char_t *histout, Float_t scale1=1.0, Float_t scale2=1.0)
{//added?
  TH2F *hist1=(TH2F *) gROOT->FindObject(histin1);
  TH2F *hist2=(TH2F *) gROOT->FindObject(histin2);
  TH2F *hist3=(TH2F *) gROOT->FindObject(histout);
  hist3->Add(hist1,hist2,scale1,scale2);
  hist3->Draw();
}

void newcanvas(Char_t *name)
{
  TCanvas *canvasname=new TCanvas(name,name,1);
  canvasname->cd();
}

void browse()
{
  TBrowser *b=new TBrowser();
}

void cutg(Char_t *histname,Char_t * newcutname, Int_t graphit=1, Int_t npts=0, Char_t *xvar="", Char_t *yvar="")
{
  if (graphit==1) {
      
    if (!(TCanvas *) gROOT->GetListOfCanvases()->FindObject("c1"))
      TCanvas *c1=new TCanvas("c1","c1");
    c1->cd();
    if ((TCutG *) gROOT->GetListOfSpecials()->FindObject("CUTG"))
      CUTG->Delete();
    if ((TCutG *) gROOT->GetListOfSpecials()->FindObject(newcutname)) {
      (TCutG *) gROOT->GetListOfSpecials()->FindObject(newcutname)->Delete();
    }
    TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
    hist2->Draw("col2");
    c1->EditorBar();
    cout << "Click on \"Graphical Cut\" to create cut"<<endl;
    cout << "When finished don't forget to setname(\"CUTG\",\""<<newcutname<<"\")"<<endl;
  }
  else {
    if (npts==0) {
      cout << "enter number of points in graphical cut: ";
      cin >> npts;
    }
    TCutG *newcut=new TCutG(newcutname,npts);
    Double_t x,y;
    for (Int_t i=0; i<npts-1; i++) {
      cout <<"enter x, y for point "<<i<<", one per line:"<<endl;
      cin>>x>>y;
      newcut->SetPoint(i,x,y);
    }
    newcut->GetPoint(0,x,y);
    newcut->SetPoint(npts-1,x,y);
    newcut->Draw();
    newcut->SetVarX(xvar);
    newcut->SetVarY(yvar);
    //      newcut->Print();
    newcut->SetName(newcutname);      
  }
}   

void dir(void)
{
  gROOT->GetListOfFiles()->Print();
}

void divide(Char_t *histin1, Char_t *histin2, Char_t *histout, Int_t ierr=0)
{
  // ierr=-1: set all resultant errors to 0
  // ierr=0: no recalculation of errors
  // ierr=1: use statistical errors on histogram 1
  // ierr=2: use statistical errors on histogram 2
  // ierr=3: combine statistical errors from histograms 1 and 2
  // ierr=4: use existing errors from histogram 1
  // ierr=5: use existing errors from histogram 2
  // ierr=6: combine existing errors from histograms 1 and 2

  Float_t dy1,dy2,dy3,y1,y2,y3;
  TH1F *hist1=(TH1F *) gROOT->FindObject(histin1);
  TH1F *hist2=(TH1F *) gROOT->FindObject(histin2);
  TH1F *hist3=(TH1F *) gROOT->FindObject(histout);
  hist3->Divide(hist1,hist2);
  if (ierr!=0) {
    // recalculate errors assuming statistical or existing errors
    for (Int_t i=0; i<hist1->GetNbinsX(); i++) {
      y1 = hist1->GetBinContent(i);
      y2 = hist2->GetBinContent(i);
      y3 = hist3->GetBinContent(i);
      dy1= hist1->GetBinError(i);
      dy2= hist2->GetBinError(i);
      switch (ierr) {
      case -1:
	dy3 = 0;
	break;
      case 1:
	if (y1!=0) {
	  dy3=TMath::Sqrt(1/y1)*y3;
	} else {
	  dy3=0;
	}
	break;
      case 2:
	if (y2 !=0) {
	  dy3=TMath::Sqrt(1/y2)*y3;
	} else {
	  dy3=0.;
	}
	break;
      case 3:
	if (y1!=0 && y2 !=0) {
	  dy3=TMath::Sqrt(1/y1 + 1/y2)*y3;
	} else if (y1!=0 && y2==0) {
	  dy3=TMath::Sqrt(1/y1) * y3;
	} else if (y1==0 && y2 !=0) {
	  dy3=TMath::Sqrt(1/y2) * y3;
	}
	break;
      case 4:
	if (y1!=0) {
	  dy3=(dy1/y1) * y3;
	} else {
	  dy3=0.;
	}
	break;
      case 5:
	if (y2 !=0) {
	  dy3=(dy2/y2) * y3;
	} else {
	  dy3=0.;
	}
	break;
      case 6:
	//	cout << i << " "<<y1<<" "<<dy1<<" "<<y2<<" "<<dy2<<" "<<y3<<endl;
	if (y1!=0 && y2 !=0) {
	  dy3=TMath::Sqrt((dy1*dy1)/(y1*y1) + (dy2*dy2)/(y2*y2))*y3;
	} else if (y1!=0 && y2==0) {
	  dy3=(dy1/y1) * y3;
	} else if (y1==0 && y2 !=0) {
	  dy3=(dy2/y2) * y3;
	}
	break;
      }
      hist3->SetBinError(i,dy3);
    }
    hist3->Draw("E");
  } else {		     
    hist3->Draw("HIST");
  }
}

void draw(Char_t *histname,Char_t *dopt="",Float_t xmin=-999999.,Float_t xmax=999999.)
{
  Int_t minbin,maxbin;
  Int_t lowbin,highbin;
  Float_t xlow,xhigh;
   
  TH1F *hist1=(TH1F *) gROOT->FindObject(histname);
  lowbin=hist1->GetXaxis()->GetFirst();
  highbin=hist1->GetXaxis()->GetLast();
  if (xmin==-999999. && xmax==999999.) {
    hist1->GetXaxis()->UnZoom();//added
    hist1->GetXaxis()->SetRange(lowbin,highbin);
  } else if (xmax==999999. && xmin !=-999999.) {
    minbin=hist1->FindBin(xmin);
    hist1->GetXaxis()->SetRange(minbin,highbin);
  } else {
    minbin=hist1->FindBin(xmin);
    maxbin=hist1->FindBin(xmax);
    hist1->GetXaxis()->SetRange(minbin,maxbin);
  }
  hist1->Draw(dopt);
}

void draw2(Char_t *histname,Float_t xmin=-999999.,Float_t xmax=999999.,
	   Float_t ymin=-999999.,Float_t ymax=999999.)
{
  Int_t minxbin,maxxbin;
  Int_t lowxbin,highxbin;
  Float_t xlow,xhigh;
  Int_t minybin,maxybin;
  Int_t lowybin,highybin;
  Float_t ylow,yhigh;
   
  TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
   
  lowxbin=hist2->GetXaxis()->GetFirst();
  highxbin=hist2->GetXaxis()->GetLast();
  lowybin=hist2->GetYaxis()->GetFirst();
  highybin=hist2->GetYaxis()->GetLast();
   
  if (xmin==-999999. && xmax==999999.) {
    hist2->GetXaxis()->UnZoom();//added
    hist2->GetXaxis()->SetRange(lowxbin,highxbin);
  } else if (xmax==999999. && xmin !=-999999.) {
    minxbin=hist2->GetXaxis()->FindBin(xmin);
    hist2->GetXaxis()->SetRange(minxbin,highxbin);
  } else {
    minxbin=hist2->GetXaxis()->FindBin(xmin);
    maxxbin=hist2->GetXaxis()->FindBin(xmax);
    hist2->GetXaxis()->SetRange(minxbin,maxxbin);
  }
  if (ymin==-999999. && ymax==999999.) {
    hist2->GetYaxis()->UnZoom();//added    
    hist2->GetYaxis()->SetRange(lowybin,highybin);
  } else if (ymax==999999. && ymin !=-999999.) {
    minybin=hist2->GetYaxis()->FindBin(ymin);
    hist2->GetYaxis()->SetRange(minybin,highybin);
  } else {
    minybin=hist2->GetYaxis()->FindBin(ymin);
    maxybin=hist2->GetYaxis()->FindBin(ymax);
    hist2->GetYaxis()->SetRange(minybin,maxybin);
  }
  hist2->Draw("col2");
}

void drawline(Char_t *filename, Char_t *linename="gline",Bool_t showpoints=1,Int_t linestyle=1,Int_t linewidth=1)
{//new
  //Int_t maxpoints=100;//not used
  Int_t point=0;
  Float_t x,y;
  //Int_t npts=0;//not used

  TGraph *graph = new TGraph(2);
  graph->SetName(linename);
  graph->SetTitle(linename);
  graph->SetFillColor(1);
  graph->SetLineStyle(linestyle);
  graph->SetLineWidth(linewidth);

  ifstream infile(filename);
  while (infile >> x) {
    infile>>y;
    graph->SetPoint(point,x,y);
    if(showpoints)
      cout<<"reading point "<<point<<endl;
    point++;
  }
  //   point=point/2;
  cout<<"points = "<<point<<endl;
  
  /*   for (Int_t i=point; i<100; i++) {
       graph->RemovePoint(i);
       npts=graph->GetN();
       cout<<"removed point "<<i<<" npoints = "<<npts<<endl;
       }
  */
  if(showpoints)  
    graph->Print();
  //  cout<<"number of points= "<<npts<<endl;
  graph->Draw("l");
}


void dump(Char_t *histname, Char_t *filename, Int_t ierr=0, Int_t zsupp=1)
{
  //   gDirectory->pwd();
  Int_t whatsit=0;
  whatsit=whatis(histname);
   
  switch(whatsit) {
  case 1:
    TH1F *hist1=(TH1F *) gROOT->FindObject(histname);
    break;
  case 2:
    TH1D *hist1=(TH1D *) gROOT->FindObject(histname);
    break;
  case 3:
    TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
    break;
  case 4:
    TH2D *hist2=(TH2D *) gROOT->FindObject(histname);
    break;
  case 5:
    TH3F *hist3=(TH3F *) gROOT->FindObject(histname);
    break;
  case 6:
    TProfile *hist1=(TProfile *) gROOT->FindObject(histname);
    break;
  }
   
  ofstream outfile(filename);
   
  if (whatsit==1 || whatsit==2 || whatsit==6) {//1D object
      
    for (Int_t ibin=0; ibin < hist1->GetNbinsX()+1;ibin++) {//includes underflow bin
      if (ierr==0) {
	if ((zsupp==1 && hist1->GetBinContent(ibin)!=0) ||
	    zsupp==0) {
	  outfile<< hist1->GetBinCenter(ibin) << " " <<
	    hist1->GetBinContent(ibin)<<endl;
	}
      } else {
	if ((zsupp==1 && hist1->GetBinContent(ibin)!=0) ||
	    zsupp==0) {
	  outfile<<hist1->GetBinCenter(ibin)<<" "<<
	    hist1->GetBinContent(ibin)<<" "<<
	    hist1->GetBinError(ibin)<<endl;
	}
      }
    }
  }
   
  if (whatsit==3 ||  whatsit==4) {//2D object
    for (Int_t ixbin=0; ixbin < hist2->GetNbinsX();ixbin++) {
      for (Int_t iybin=0; iybin<hist2->GetNbinsY();iybin++) {
	    
	Float_t content=hist2->GetBinContent(ixbin,iybin);
	Float_t xbin=hist2->GetXaxis()->GetBinCenter(ixbin);
	Float_t ybin=hist2->GetYaxis()->GetBinCenter(iybin);
	if (ierr==0) {
	  if ((zsupp==1 && content!=0) ||
	      zsupp==0) {
	    outfile<< xbin << " " << ybin <<" "<<content<<endl;
	  }
	} else {
	  if ((zsupp==1 && content!=0) ||
	      zsupp==0) {
	    outfile<<xbin<<" "<< ybin<<" "<<content<<endl;
	  }
	}
      }
    }
  }
  outfile.close();
}

void fillhist0(Char_t *filename, Char_t *histname)
{
  Float_t x;
  TH1F *hist1=(TH1F *) gROOT->FindObject(histname);
  hist1->Reset();
  ifstream infile(filename);
  while (infile >> x) {
    hist1->Fill(x);
  }
  hist1->Draw();
}

void fillhist(Char_t *filename, Char_t *histname, Int_t opt1=0, Int_t ierr=0)
{
  // if opt1=0 use x as channel to fill, else use bin number
  // if ierr=1 also read in errors and set bin errors

  Float_t x,y,dy,xval;
  Int_t xbin=1;
  TH1F *hist1=(TH1F *) gROOT->FindObject(histname);
  hist1->Reset();
  ifstream infile(filename);
  switch (ierr) {
  case 0:
    while (infile >> x) {
      infile >> y;
      xval=x;
      if (opt1==1) {xval=(Float_t) xbin;}
      hist1->Fill(xval,y);
      xbin++;
    }
    hist1->Draw();
    break;
  case 1:
    while (infile >> x) {
      infile>>y;
      infile>>dy;
      xval=x;
      if (opt1==1) {xval=(Float_t) xbin;}
      hist1->Fill(xval,y);
      hist1->SetBinError(xbin,dy);
      xbin++;
    }
    hist1->Draw("E");
    break;
  }
}

void fillhist2(Char_t *filename,Char_t *histname)
{
  Float_t x,y,z;
  TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
  hist2->Reset();
  ifstream infile(filename);
  while (infile >> x) {
    infile >> y;
    infile >> z;
    hist2->Fill(x,y,z);
  }
  hist2->Draw("col2");
}

void fillgraph(Char_t *filename, Char_t *graphname, Int_t npts, Int_t ierr=0)
{
  gr=new TGraphErrors(npts);
  Double_t x,y,dx,dy;
  ifstream infile(filename);
  Int_t i=0;
  while (infile >> x) {
    infile >> y;
    infile >> dy;
    gr->SetPoint(i,x,y);
    gr->SetPointError(i,dx,dy);
    i++;
  }
  gr->SetName(graphname);
  gr->SetTitle("graph test");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("AP");
}

void getfile(Char_t *filename)
{
  TFile *f=new TFile(filename);
}

void setdef(Char_t *dirname)
{
  gDirectory->cd(dirname);
  cout << "Current directory is " <<gDirectory->pwd()<<endl;
}

void h1(Char_t *histname, Char_t *title, Int_t nchan=100,Float_t low=0,Float_t high=100)
{
  TH1F *newhist=new TH1F(histname,title,nchan,low,high);
  newhist->SetName(histname);
}

void h2(Char_t *histname, Char_t *title, Int_t nxchan=100,Float_t lowx=0,Float_t highx=100,
	Int_t nychan,Float_t lowy=0,Float_t highy=100)
{
  TH2F *newhist=new TH2F(histname,title,nxchan,lowx,highx,nychan,lowy,highy);
  newhist->SetName(histname);
}

void h3(Char_t *histname, Char_t *title, Int_t nxchan=100,Float_t lowx=0,Float_t highx=100,
	Int_t nychan,Float_t lowy=0,Float_t highy=100, Int_t nzchan=100, 
	Float_t lowz=0, Float_t highz=100)
{
  TH3F *newhist=new TH3F(histname,title,nxchan,lowx,highx,nychan,lowy,highy,
			 nzchan,lowz,highz);
  newhist->SetName(histname);
}

void listall(Int_t option=0)
{
  switch (option) {
      
  case 0:
    cout << endl<<"All histograms in memory:"<<endl<<endl;
    gDirectory->GetList()->ls();
    cout << endl<<"All Special objects in memory:"<<endl<<endl;;
    gROOT->GetListOfSpecials()->ls();
    cout << endl<<"All Canvases in memory:"<<endl<<endl;
    gROOT->GetListOfCanvases()->ls();
    break;
      
  case 1:
    cout << endl<<"All histograms in memory:"<<endl<<endl;
    gDirectory->GetList()->ls();
    break;
      
  case 2:
    cout << endl<<"All Special objects in memory:"<<endl<<endl;;
    gROOT->GetListOfSpecials()->ls();
    break;
      
  case 3:
    cout << endl<<"All Canvases in memory:"<<endl<<endl;
    gROOT->GetListOfCanvases()->ls();
    break;
  }
      
}

void pj3x(Char_t *histname, Float_t ylow=-999999., Float_t yhi=999999., Float_t zlow=-999999., Float_t zhi=999999.)
{
  Int_t ylobin,yhibin,zlobin,zhibin;
  TH3F *hist3=(TH3F *) gROOT->FindObject(histname);
  Int_t nch=strlen(histname+3);
  char * projname=new char[nch];
  sprintf(projname,"%s_%s",histname,"x"); 

  if (ylow==-999999. && yhi=999999) {
    ylow=hist3->GetYaxis()->GetXmin();
    yhi=hist3->GetYaxis()->GetXmax();
  }
  if (zlow==-999999. && zhi=999999) {
    zlow=hist3->GetZaxis()->GetXmin();
    zhi=hist3->GetZaxis()->GetXmax();
  }

  ylobin=hist3->GetYaxis()->FindBin(ylow);
  yhibin=hist3->GetYaxis()->FindBin(yhi);
  zlobin=hist3->GetZaxis()->FindBin(zlow);
  zhibin=hist3->GetZaxis()->FindBin(zhi);

  hist3->GetYaxis()->SetRange(ylobin,yhibin);
  hist3->GetZaxis()->SetRange(zlobin,zhibin);
  hist3->GetXaxis()->SetRange(0,hist3->GetNbinsX());
  hist3->Project3D("x")->Draw();
}

void pj3y(Char_t *histname, Float_t xlow, Float_t xhi, Float_t zlow, Float_t zhi)
{
  TH3F *hist3=(TH3F *) gROOT->FindObject(histname);
  hist3->GetXaxis()->SetRange(xlow,xhi);
  hist3->GetZaxis()->SetRange(zlow,zhi);
  hist3->Project3D("y");
}

void pj3z(Char_t *histname, Float_t xlow, Float_t xhi, Float_t ylow, Float_t yhi)
{
  TH3F *hist3=(TH3F *) gROOT->FindObject(histname);
  hist3->GetXaxis()->SetRange(xlow,xhi);
  hist3->GetYaxis()->SetRange(ylow,yhi);
  hist3->Project3D("z");
}

void pj3xy(Char_t *histname, Float_t zmin=-999999., Float_t zmax=999999.)
{
  Int_t zbin1,zbin2;
  TH3F *hist3=(TH3F *) gROOT->FindObject(histname);
  // search for an existing projection and delete it if it exists:
  Int_t nch=strlen(histname+4);
  char * projname=new char[nch];
  sprintf(projname,"%s_%s",histname,"xy"); 

  if (zmin==-999999. && zmax==999999.) {
    zmin=hist3->GetZaxis()->GetXmin();
    zmax=hist3->GetZaxis()->GetXmax();
  } else if (zmax==999999. && zmin !=-999999.) {
    zmax=zmin;
  } 
     
  zbin1=hist3->GetZaxis()->FindBin(zmin);
  zbin2=hist3->GetZaxis()->FindBin(zmax);
   
  hist3->GetZaxis()->SetRange(zbin1,zbin2);
     
  hist3->GetXaxis()->SetRange(0,hist3->GetNbinsX());
  hist3->GetYaxis()->SetRange(0,hist3->GetNbinsY());
   
  hist3->Project3D("xy")->Draw("col2");
}

void pj3yx(Char_t *histname, Float_t zmin=-999999., Float_t zmax=999999.)
{
  Int_t zbin1,zbin2;
  TH3F *hist3=(TH3F *) gROOT->FindObject(histname);
  // search for an existing projection and delete it if it exists:
  Int_t nch=strlen(histname+4);
  char * projname=new char[nch];
  sprintf(projname,"%s_%s",histname,"xy"); 
  //   if ((TH2D *) gROOT->FindObject(projname)) 
  //     (TH2D *) gROOT->FindObject(projname)->Delete();     
   
  if (zmin==-999999. && zmax==999999.) {
    zmin=hist3->GetZaxis()->GetXmin();
    zmax=hist3->GetZaxis()->GetXmax();
  } else if (zmax==999999. && zmin !=-999999.) {
    zmax=zmin;
  } 
     
  zbin1=hist3->GetZaxis()->FindBin(zmin);
  zbin2=hist3->GetZaxis()->FindBin(zmax);
   
  hist3->GetZaxis()->SetRange(zbin1,zbin2);
   
  //   hist3->GetYaxis()->SetRange(hist3->GetYaxis()->GetXmin(),
  //			       hist3->GetYaxis()->GetXmax());
 
  //   hist3->GetXaxis()->SetRange(hist3->GetXaxis()->GetXmin(),
  //			       hist3->GetXaxis()->GetXmax());
  
  hist3->GetXaxis()->SetRange(0,hist3->GetNbinsX());
  hist3->GetYaxis()->SetRange(0,hist3->GetNbinsY());
   
  hist3->Project3D("yx")->Draw("col2");
}

void pjxwin(Char_t *histname, Char_t *cutname)
{
  TCutG *tempcut;
  Int_t whatsit=whatis(histname,0);
  if (whatsit != 3 && whatsit !=4) {
    cout << "histogram is not a 2D or is not found!"<<endl;
    return;
  }
  switch(whatsit) {
  case 3:
    TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
    break ;
  case 4:
    TH2D *hist2=(TH2D *) gROOT->FindObject(histname);
    break;
  }
   
  if (!(TCutG *) gROOT->GetListOfSpecials()->FindObject(cutname)) {
    cout << "Graphical cut "<<cutname<<" is not found!"<<endl;
    return;
  } else {
    tempcut = (TCutG *) gROOT->GetListOfSpecials()->FindObject(cutname);
  }
   
  if ((TH1D *) gROOT->FindObject("xwproj")) xwproj->Delete();
  TH1D * xwproj=new TH1D("xwproj","X projection",hist2->GetNbinsX(),
			 hist2->GetXaxis()->GetXmin(),
			 hist2->GetXaxis()->GetXmax());
   
  for (Int_t xbin=0; xbin < hist2->GetNbinsX();xbin++) {
    for (Int_t ybin=0; ybin < hist2->GetNbinsY(); ybin++) {
      if (tempcut->IsInside(hist2->GetXaxis()->GetBinCenter(xbin),
			    hist2->GetYaxis()->GetBinCenter(ybin))) {
	xwproj->Fill(hist2->GetXaxis()->GetBinCenter(xbin),
		     hist2->GetBinContent(xbin,ybin));
      }
    }
  }
  xwproj->Draw();
}	 

void pjywin(Char_t *histname, Char_t *cutname)
{
  TCutG *tempcut;
  Int_t whatsit=whatis(histname,0);
  if (whatsit != 3 && whatsit !=4) {
    cout << "histogram is not a 2D or is not found!"<<endl;
    return;
  }
  switch(whatsit) {
  case 3:
    TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
    break ;
  case 4:
    TH2D *hist2=(TH2D *) gROOT->FindObject(histname);
    break;
  }
   
  if (!(TCutG *) gROOT->GetListOfSpecials()->FindObject(cutname)) {
    cout << "Graphical cut "<<cutname<<" is not found!"<<endl;
    return;
  } else {
    tempcut = (TCutG *) gROOT->GetListOfSpecials()->FindObject(cutname);
  }
   
  if ((TH1D *) gROOT->FindObject("ywproj")) ywproj->Delete();
  TH1D * ywproj=new TH1D("ywproj","Y projection",hist2->GetNbinsY(),
			 hist2->GetYaxis()->GetXmin(),
			 hist2->GetYaxis()->GetXmax());
   
  for (Int_t xbin=0; xbin < hist2->GetNbinsX();xbin++) {
    for (Int_t ybin=0; ybin < hist2->GetNbinsY(); ybin++) {
      if (tempcut->IsInside(hist2->GetXaxis()->GetBinCenter(xbin),
			    hist2->GetYaxis()->GetBinCenter(ybin))) {
	ywproj->Fill(hist2->GetYaxis()->GetBinCenter(ybin),
		     hist2->GetBinContent(xbin,ybin));
      }
    }
  }
  ywproj->Draw();
}	 

void pj2win(Char_t *histname, Char_t *cutname)
{
  TCutG *tempcut;
  Int_t whatsit=whatis(histname,0);
  if (whatsit != 3 && whatsit !=4) {
    cout << "histogram is not a 2D or is not found!"<<endl;
    return;
  }
  switch(whatsit) {
  case 3:
    TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
    break ;
  case 4:
    TH2D *hist2=(TH2D *) gROOT->FindObject(histname);
    break;
  }
   
  if (!(TCutG *) gROOT->GetListOfSpecials()->FindObject(cutname)) {
    cout << "Graphical cut "<<cutname<<" is not found!"<<endl;
    return;
  } else {
    tempcut = (TCutG *) gROOT->GetListOfSpecials()->FindObject(cutname);
  }
   
  if ((TH2F *) gROOT->FindObject("xywproj")) xywproj->Delete();
  TH2F * xywproj=new TH2F("xywproj","XY copy",hist2->GetNbinsX(),
			  hist2->GetXaxis()->GetXmin(),
			  hist2->GetXaxis()->GetXmax(),
			  hist2->GetNbinsY(),
			  hist2->GetYaxis()->GetXmin(),
			  hist2->GetYaxis()->GetXmax());
   
  for (Int_t xbin=0; xbin < hist2->GetNbinsX();xbin++) {
    for (Int_t ybin=0; ybin < hist2->GetNbinsY(); ybin++) {
      if (tempcut->IsInside(hist2->GetXaxis()->GetBinCenter(xbin),
			    hist2->GetYaxis()->GetBinCenter(ybin))) {
	xywproj->Fill(hist2->GetXaxis()->GetBinCenter(xbin),
		      hist2->GetYaxis()->GetBinCenter(ybin),
		      hist2->GetBinContent(xbin,ybin));
      }
    }
  }
  xywproj->Draw("col2");
}	 

void pjx(Char_t *histname, Float_t ymin=-999999., Float_t ymax=999999.)
{
  Int_t ybin1,ybin2;
  if ((TH1D *) gROOT->FindObject("xproj")) xproj->Delete();
  TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
   
  if (ymin==-999999. && ymax==999999.) {
    ybin2=hist2->GetNbinsY();
    ybin1=0;
  } else if (ymax==999999. && ymin !=-999999.) {
    ybin1=hist2->GetYaxis()->FindBin(ymin);
    ybin2=ybin1;
  } else {
    ybin1=hist2->GetYaxis()->FindBin(ymin);
    ybin2=hist2->GetYaxis()->FindBin(ymax);
  }
   
  TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
  hist2->ProjectionX("xproj",ybin1,ybin2)->Draw();
}

void pjy1(Char_t *histname, Float_t x=-999999.) {
  pjy(histname, x, x);
}

void pjy(Char_t *histname, Float_t xmin=-999999., Float_t xmax=999999.)
{
  Int_t xbin1,xbin2;
  if ((TH1D *) gROOT->FindObject("yproj")) yproj->Delete();

  TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
  if (xmin==-999999. && xmax==999999.) {
    xbin2=hist2->GetNbinsX();
    xbin1=0;
  } else if (xmax==999999. && xmin !=-999999.) {
    xbin1=hist2->GetXaxis()->FindBin(xmin);
    xbin2=xbin1;
  } else {
    xbin1=hist2->GetXaxis()->FindBin(xmin);
    xbin2=hist2->GetXaxis()->FindBin(xmax);
  }
  hist2->ProjectionY("yproj",xbin1,xbin2)->Draw();
}

void savecanvas(Char_t *cnvname, Char_t * filename)
{
  TCanvas *thecanvas=(TCanvas *) gROOT->FindObject(cnvname);
  thecanvas->SaveAs(filename);
}

void doprint(Char_t * cnvname="cc", Char_t * filename="print.ps", Char_t * pr="m0-epson")
{
  const char cmd[255];
  char * fl = reinterpret_cast<char *>(filename);
  TCanvas *thecanvas=(TCanvas *)gROOT->FindObject(cnvname);
  thecanvas->SaveAs(filename);
  sprintf(cmd,"lpr -P%s %s",pr,fl);
  gSystem->Exec(cmd);
  cout<<fl<<" was sent to printer "<<pr<<endl;
}

void gfit(Char_t *histname, Float_t xmin=-999999., Float_t xmax=999999)
{
  Float_t sigma=0, width=0; 
  TH1F *hist1=(TH1F*) gROOT->FindObject(histname);
  if (xmin==-999999. && xmax==999999.) {
    xmin=hist1->GetXaxis()->GetXmin();
    xmax=hist1->GetXaxis()->GetXmax();
  } else if (xmax==999999. && xmin !=-999999.) {
    xmax=xmin;
  } 
  hist1->Fit("gaus","W","",xmin,xmax);
  sigma=hist1->GetFunction("gaus")->GetParameter(2);
  width=sigma*2.35482;
  printf("Width of peak is %f or %f FWHM\n",sigma,width);
  printf("Width of peak is %f ns or %f FWHM ns, mean %f ns\n",sigma/5.,width/5.,hist1->GetFunction("gaus")->GetParameter(1)/5.);
  printf("Width of peak is %f mm or %f FWHM mm, mean %f mm\n",sigma/5./2.5,width/5./2.5,hist1->GetFunction("gaus")->GetParameter(1)/5./2.5);
}

void pfit(Char_t *histname, Float_t xmin=-999999., Float_t xmax=999999,Int_t order=1)
{
  TH1F *hist1=(TH1F*) gROOT->FindObject(histname);
  if (xmin==-999999. && xmax==999999.) {
    xmin=hist1->GetXaxis()->GetXmin();
    xmax=hist1->GetXaxis()->GetXmax();
  } else if (xmax==999999. && xmin !=-999999.) {
    xmax=xmin;
  } 
  if (order<10) {
    switch(order) {
    case 0:
      cout <<" Not meaningful!!"<<endl;
      break;
    case 1:
      hist1->Fit("pol1","W","",xmin,xmax);
      break;
    case 2:
      hist1->Fit("pol2","W","",xmin,xmax);
      break;
    case 3:
      hist1->Fit("pol3","W","",xmin,xmax);
      break;
    case 4:
      hist1->Fit("pol4","W","",xmin,xmax);
      break;
    case 5:
      hist1->Fit("pol5","W","",xmin,xmax);
      break;
    case 6:
      hist1->Fit("pol6","W","",xmin,xmax);
      break;
    case 7:
      hist1->Fit("pol7","W","",xmin,xmax);
      break;
    case 8:
      hist1->Fit("pol8","W","",xmin,xmax);
      break;
    case 9:
      hist1->Fit("pol9","W","",xmin,xmax);
      break;
    }
  } else {
    cout << "Only works to 9th order."<<endl;
  }
}
 
void pfx(Char_t *histname,Float_t ymin=-999999.,Float_t ymax=999999.)
{
  if ((TProfile *) gROOT->FindObject("xprof")) xprof->Delete();
  Int_t yminbin,ymaxbin;
  TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
  if (ymin==-999999.) {
    yminbin=0;
  } else {
    yminbin=hist2->GetYaxis()->FindBin(ymin);
  }
  if (ymax==999999.) {
    ymaxbin=hist2->GetNbinsY();
  } else {
    ymaxbin=hist2->GetYaxis()->FindBin(ymax);
  }
  hist2->ProfileX("xprof",yminbin,ymaxbin)->Draw();
}

void pfy(Char_t *histname,Float_t xmin=0, Float_t xmax=0)
{
  if ((TProfile *) gROOT->FindObject("yprof")) yprof->Delete();
  Int_t xminbin,xmaxbin;
  TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
  if (xmin==0) {
    xminbin=0;
  } else {
    xminbin=hist2->GetXaxis()->FindBin(xmin);
  }
  if (xmax==0) {
    xmaxbin=hist2->GetNbinsX();
  } else {
    xmaxbin=hist2->GetXaxis()->FindBin(xmax);
  }
  hist2->ProfileY("yprof",xminbin,xmaxbin)->Draw();

}

void drawcut(Char_t *cutname) 
{ if( c1 == 0 ) {
    cout<<"please display the histogram first.\n";
    return;
  }
  char path[124]=gROOT->GetPath();
  readcuts("cut_file.root");
  TCutG *cut1 = (TCutG*) gROOT->FindObject(cutname);
  if( cut1 == 0 ) { 
    cout<<"Sorry, "<<cutname<<" does not exist.\n";
  } else {
    c1->cd();
    cut1->SetLineColor(2);
    cut1->SetLineWidth(2);
    cut1->Draw();
  }
  gROOT->cd(path);
}

void readcuts(Char_t *filename)
{
  char path[124]=gROOT->GetPath();//added
  cout<<path<<endl;
  TFile *fin=new TFile(filename);
  TList *speclist=(TList*) fin->GetListOfSpecials();
  speclist->Read();
  //the following 3 lines are copied from a similar block of code...
  fin->ReadAll();
  fin->Purge();
  fin->ls("-d");
  fin->Close();
  gROOT->cd(path);
}

void deletecuts(Char_t *filename, Char_t *cut_name) 
{
  char path[124]=gROOT->GetPath();
  cout<<path<<endl;
  TFile *fin=new TFile(filename,"update");
  fin->ReadAll();
  TString cutname=TString(cut_name)+";*";
  cout<<"Do you want to delete "<<cutname<<" ? (y/n)"<<endl;
  char opt;
  cin>>opt;
  if( opt == 'y' ) {
    fin->Delete( cutname.Data() );
    cout<<cutname<<" is deleted."<<endl;
  }
  fin->Close();
  gROOT->cd(path);
}

void setname(Char_t *oldname, Char_t *newname)
{
  TNamed *obj=(TNamed *) gROOT->FindObject(oldname);
  obj->SetName(newname);
}

void saveall(Char_t *filename)
{
  TList *hlist=gDirectory->GetList();
  TList *clist=gROOT->GetListOfCanvases();
  TList *speclist=gROOT->GetListOfSpecials();
  TFile *fout=new TFile(filename,"recreate");
  fout->cd();
  hlist->Write();
  clist->Write();
  speclist->Write();
  fout->Close();
  cout << "Histograms, canvases and specials written to file "<<filename<<endl;
}

void savecuts(Char_t *filename)
{
  TList *speclist=gROOT->GetListOfSpecials();
  TFile *fout=new TFile(filename,"recreate");
  fout->cd();
  speclist->Write();
  fout->Purge();//added
  fout->ls();
  fout->Close();
  cout << "Special objects written to file "<<filename<<endl;
}

void listcuts()
{    
  TList *speclist=gROOT->GetListOfSpecials();
  cout<<"To be implemented."<<endl;
}

void setcutg(Char_t *cutname="CUTG", Char_t *xvar="",Char_t *yvar="")
{
  TCutG *gcut=(TCutG *) gROOT->GetListOfSpecials()->FindObject(cutname);
  gcut->SetVarX(xvar);
  gcut->SetVarY(yvar);
}

void setplain(void)
{
  //gROOT->SetStyle("Plain");
  //gStyle->SetPalette(1,0);
  printf("Default style being used.\n");
}

void show(Char_t *histname)
{
  Int_t whatsit=whatis(histname);
  switch(whatsit) {
  case 1:
    TH1F *hist1=(TH1F *) gROOT->FindObject(histname);
    break;
  case 2:
    TH1D *hist1=(TH1D *) gROOT->FindObject(histname);
    break;
  case 3:
    TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
    break;
  case 4:
    TH2D *hist2=(TH2D *) gROOT->FindObject(histname);
    break;
  case 5:
    TH3F *hist3=(TH3F *) gROOT->FindObject(histname);
    break;
  case 6:
    TProfile *hist1=(TProfile *) gROOT->FindObject(histname);
    break;
  }
   
  if (whatsit==1||whatsit==2) {
    cout << "1D histogram "<<histname<<" : "<<endl;
    cout << hist1->GetNbinsX()<<" bins from " <<hist1->GetXaxis()->GetXmin()
	 << " to " << hist1->GetXaxis()->GetXmax()<<endl;
  }
  if (whatsit==3 || whatsit==4) {
    cout <<"2D histogram "<<histname<<" : "<<endl;
    cout << "X axis: "<<hist2->GetNbinsX()<<" bins from " <<
      hist2->GetXaxis()->GetXmin()<<" to "<<hist2->GetXaxis()->GetXmax()
	 <<endl;
    cout << "Y axis: "<<hist2->GetNbinsY()<<" bins from " <<
      hist2->GetYaxis()->GetXmin()<<" to "<<hist2->GetYaxis()->GetXmax()
	 <<endl;
  }
}

void subtract(Char_t *histin1, Char_t *histin2, Char_t *histout, Float_t scale1=1.0, Float_t scale2=-1.0)
{//updated to use scale, omitts ierr

  TH1F *hist1=(TH1F *) gROOT->FindObject(histin1);
  TH1F *hist2=(TH1F *) gROOT->FindObject(histin2);
  TH1F *hist3=(TH1F *) gROOT->FindObject(histout);
  hist3->Add(hist1,hist2,scale1,scale2);
  hist3->Draw();
}

void subtract2(Char_t *histin1, Char_t *histin2, Char_t *histout, Float_t scale1=1.0, Float_t scale2=-1.0)
{//added

  TH2F *hist1=(TH2F *) gROOT->FindObject(histin1);
  TH2F *hist2=(TH2F *) gROOT->FindObject(histin2);
  TH2F *hist3=(TH2F *) gROOT->FindObject(histout);
  hist3->Add(hist1,hist2,scale1,scale2);
  hist3->Draw();
}

void sum(Char_t *histname,Float_t xmin=-999999.,Float_t xmax=999999.,Float_t ymin=-999999.,Float_t ymax=999999.)
{
  Int_t minbin,maxbin;
  Int_t lowbin,highbin;
  Float_t sum,centroid,RMS,background=0.;
  Float_t xlow,xhigh;
  Bool_t usebkg=kTRUE;
  if(ymin==-999999. || ymax==999999.) usebkg=kFALSE;
   
  TH1F *hist1=(TH1F *) gROOT->FindObject(histname);
  lowbin=hist1->GetXaxis()->GetFirst();
  highbin=hist1->GetXaxis()->GetLast();
  if (xmin==-999999. && xmax==999999.) {
    //      maxbin=hist1->GetNbinsX();
    //      minbin=0;
    minbin=lowbin;
    maxbin=highbin;
  } else if (xmax==999999. && xmin !=-999999.) {
    minbin=hist1->FindBin(xmin);
    maxbin=minbin;
  } else {
    minbin=hist1->FindBin(xmin);
    maxbin=hist1->FindBin(xmax);
  }
  xlow=hist1->GetXaxis()->GetBinCenter(minbin);
  xhigh=hist1->GetXaxis()->GetBinCenter(maxbin);
  hist1->GetXaxis()->SetRange(minbin,maxbin);
  centroid=hist1->GetMean();
  RMS=hist1->GetRMS();
  if (usebkg) background=(maxbin-minbin+1)*0.5*(ymax+ymin);
  sum=hist1->Integral(minbin,maxbin)-background;
  printf("%.2f - %.2f : Sum=%.1f Bkg= %.1f Centroid=%.3f RMS=%.3f\n",xlow,xhigh,sum,background,centroid,RMS); 
  hist1->GetXaxis()->SetRange(lowbin,highbin);
}

void sum2(Char_t *histname,Float_t xmin=-999999.,Float_t xmax=999999.,Float_t ymin=-999999.,Float_t ymax=999999.)
{
  Int_t minxbin,maxxbin;
  Int_t minybin,maxybin;
  Int_t lowxbin,highxbin;
  Int_t lowybin,highybin;
   
  Float_t sum,centroidx,RMSx,centroidy,RMSy;
  Float_t xlow,xhigh;
  Float_t ylow,yhigh;
   
  TH2F *hist2=(TH2F *) gROOT->FindObject(histname);
  lowxbin=hist2->GetXaxis()->GetFirst();
  highxbin=hist2->GetXaxis()->GetLast();
  lowybin=hist2->GetYaxis()->GetFirst();
  highybin=hist2->GetYaxis()->GetLast();
   
  if (xmin==-999999. && xmax==999999.) {
    //      maxxbin=hist2->GetNbinsX();
    //      minxbin=0;
    minxbin=lowxbin;
    maxxbin=highxbin;
  } else if (xmax==999999. && xmin !=-999999.) {
    minxbin=hist2->GetXaxis()->FindBin(xmin);
    maxxbin=minxbin;
  } else {
    minxbin=hist2->GetXaxis()->FindBin(xmin);
    maxxbin=hist2->GetXaxis()->FindBin(xmax);
  }

  if (ymin==-999999. && ymax==999999.) {
    //      maxybin=hist2->GetNbinsY();
    //      minybin=0;
    minybin=lowybin;
    maxybin=highybin;
  } else if (ymax==999999. && ymin !=-999999.) {
    minybin=hist2->GetYaxis()->FindBin(ymin);
    maxybin=minybin;
  } else {
    minybin=hist2->GetYaxis()->FindBin(ymin);
    maxybin=hist2->GetYaxis()->FindBin(ymax);
  }

   
  xlow=hist2->GetXaxis()->GetBinCenter(minxbin);
  xhigh=hist2->GetXaxis()->GetBinCenter(maxxbin);
  ylow=hist2->GetYaxis()->GetBinCenter(minybin);
  yhigh=hist2->GetYaxis()->GetBinCenter(maxybin);
   
  hist2->GetXaxis()->SetRange(minxbin,maxxbin);
  hist2->GetYaxis()->SetRange(minybin,maxybin);
   
  //   cout << "X axis bins: " << minxbin << " to " << maxxbin<<endl;
  //   cout << "Y axis bins: " << minybin << " to " << maxybin<<endl;

  sum=hist2->Integral(minxbin,maxxbin,minybin,maxybin);
  centroidx=hist2->GetMean(1);
  centroidy=hist2->GetMean(2);
  RMSx=hist2->GetRMS(1);
  RMSy=hist2->GetRMS(2);
  printf("%7.2f - %7.2f : Sum=%9.1f XCentroid=%.3f RMS=%.3f\n",xlow,xhigh,sum,centroidx,RMSx); 
  printf("%7.2f - %7.2f : Sum=%9.1f YCentroid=%.3f RMS=%.3f\n",ylow,yhigh,sum,centroidy,RMSy); 
  hist2->GetXaxis()->SetRange(lowxbin,highxbin);
  hist2->GetYaxis()->SetRange(lowybin,highybin);
}


Int_t whatis(Char_t *hname,Int_t verbose=1)
{
  Int_t returnvalue=0;
   
  //   gDirectory->pwd();
  //   if (gROOT->FindObject(hname)) {
  //      cout << "found it somewhere."<<endl;
  //   }

  if (gROOT->FindObject(hname)->InheritsFrom("TH1F")) {
    if (verbose==1) cout << "histogram " <<hname<<" is a TH1F"<<endl;
    returnvalue=1;
  }
  if (gROOT->FindObject(hname)->InheritsFrom("TH1D")) {
    if (verbose==1) cout << "histogram " <<hname<<" is a TH1D"<<endl;
    returnvalue=2;
  }
  if (gROOT->FindObject(hname)->InheritsFrom("TH2F")) {
    if (verbose==1) cout << "histogram "<< hname <<" is a TH2F"<<endl;
    returnvalue=3;
  }
  if (gROOT->FindObject(hname)->InheritsFrom("TH2D")) {
    if (verbose==1) cout << "histogram "<<hname<<"is a TH2D"<<endl;
    returnvalue=4;
  }
  if (gROOT->FindObject(hname)->InheritsFrom("TH3F")) {
    if (verbose==1) cout << "histogram "<<hname<<"is a TH3F"<<endl;
    returnvalue=5;
  }
  if (gROOT->FindObject(hname)->InheritsFrom("TProfile")) {
    if (verbose==1) cout << "histogram "<<hname<<" is a TProfile"<<endl;
    returnvalue=6;
  }
  return returnvalue;
}

void where(void)
{
  cout << "Current directory is "<<gDirectory->pwd()<<endl;
}

void help_util(void)
{
   
  Char_t junk;
  cout <<endl<< "Histogram utility options:"<<endl<<endl;
  cout << "cutg(\"histname\",\"cutname\"):"<<endl<<
    " Create a graphical cut (2D condition) named \"cutname\" using the"<<endl<<
    " histogram \"histname\".  The default name is CUTG and you will be"<<endl<<
    " prompted to rename the cut."<<endl;
  cout << "draw(\"histname\"[,xmin][,xmax]):"<<endl<<
    "  Draw 1D histogram histname between xmin and xmax"<<endl;
  cout << "draw2(\"histname\"[,xmin][,xmax][,ymin][,ymax]):"<<endl<<
    "  Draw 2D histogram histname between xmin and xmax,ymin and ymax"<<endl;
  cout <<"fillhist(\"histname\",\"filename\")"<<endl<<
    "  Fill the 1D histogram \"histname\" with the contents of the file \"filename\"."<<endl;
  cout << "gfit(\"histname\"[,xmin][,xmax]):"<<endl<<
    "  Fit histogram histname to a Gaussian between xmin and xmax."<<endl<<
    "  If xmin and xmax are not given the fit is over the entire"<<endl<<
    "  histogram range."<<endl;
  cout << "listall():"<<endl<<
    "  List all histograms, cuts, and canvases stored in memory."<<endl;
  cout << endl<<"Hit any character followed by return for more..."<<endl;
  cin >> junk;

  cout << "pfit(\"histname\"[,xmin][,xmax][,order]):"<<endl<<
    "  Fit histogram histname to an nth order polynomial, n given by order,"<<endl<< 
    "  between xmin and xmax."<<endl<<
    "  If xmin and xmax are not given the fit is over the entire"<<endl<<
    "  histogram range."<<endl;
  cout << "pfx(\"histname\"[,ymin][,ymax]):"<<endl<<
    "  Profile histogram histname onto X axis. If ymax is not given the"<<endl<<
    "  profile is for a single bin at ymin. If neither ymax nor ymin"<<endl<<
    "  are given the profile is for the entire y range."<<endl;     
  cout << "pfy(\"histname\"[,xmin][,xmax]):"<<endl<<
    "  Profile histogram histname onto Y axis. If xmax is not given the"<<endl<<
    "  profile is for a single bin at xmin. If neither xmax nor xmin"<<endl<<
    "  are given the profile is for the entire x range."<<endl;        
  cout << "pjx(\"histname\"[,ymin][,ymax]):"<<endl<<
    "  Project histogram histname onto X axis. If ymax is not given the"<<endl<<
    "  projection is for a single bin at ymin. If neither ymax nor ymin"<<endl<<
    "  are given the projection is for the entire y range."<<endl;     
  cout << "pj2win(\"histname\",\"cutname\"):"<<endl<<
    "  Copy the contents of the 2D histograms \"histname\" that are within"<<endl<<
    "  the 2D window \"cutname\" into a new histogram."<<endl;
  cout << endl<<"Hit any character followed by return for more..."<<endl;
  cin >> junk;
		 
  cout << "pjxwin(\"histname\",\"cutname\"):"<<endl<<
    "  Project the contents of \"histname\" that are inside the 2D cut called"<<endl<<
    "  \"cutname\" onto the X axis."<<endl;
  cout << "pjy(\"histname\"[,xmin][,xmax]):"<<endl<<
    "  Project histogram histname onto Y axis. If xmax is not given the"<<endl<<
    "  projection is for a single bin at xmin. If neither xmax nor xmin"<<endl<<
    "  are given the projection is for the entire x range."<<endl;
  cout << "pjywin(\"histname\",\"cutname\"):"<<endl<<
    "  Project the contents of \"histname\" that are inside the 2D cut called"<<endl<<
    "  \"cutname\" onto the Y axis."<<endl;
  cout << "saveall(\"filename\"):"<<endl<<
    "  Saves all current histograms and canvases in a file called \"filename\". "<<endl<<
    "  (You should use .root as an extension.)"<<endl;
  cout << "setplain():"<<endl<<
    "  Sets the display to a nice set of defaults."<<endl;
  cout << "sum(\"histname\"[,xmin][,xmax]):"<<endl<<
    "  Sum histogram histname from xmin to xmax.  If xmax is not given the"<<endl<<
    "  sum is for a single bin at xmin.  If neither xmin nor xmax are"<<endl<<
    "  given the sum is for the entire histogram."<<endl;
  cout << endl<<"Hit any character followed by return for more..."<<endl;
  cin >> junk;
  cout << "whatis(\"histname\"):"<<endl<<
    "  Tells you what kind of histogram \"histname\" is."<<endl;
  cout << "zap(\"histname\"):"<<endl<<
    "  Clears the contents of histogram \"histname\"."<<endl;
  cout << "zapall():"<<endl<<
    "  Clears the contents of all histograms in memory."<<endl;
   
  cout <<endl<<"REMEMBER, all histogram and file names must be enclosed"<<endl
       <<" in double quotation marks!"<<endl;
} 
