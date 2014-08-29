/*---------------------------------------------------------------------------------------------*/
//Copied from fit.cc for use with Scarlet at Argonne

void hdr(Char_t *histname,Float_t xmin=-999999.,Float_t xmax=999999.,Float_t ymin=-999999.,
	 Float_t ymax=999999., Bool_t clear=1)
{//Extension to dr() for online use.
  hup();
  dr(histname,xmin,xmax,ymin,ymax,1);
}

void hlowstat(Char_t *histname, Int_t style=20, Int_t size=1, Int_t color=1)
{//Extension to lowstat() for online use.
  hup();
  lowstat(histname,style,size,color);
}

void hplotall(Char_t *histin,Char_t *suffix="", Bool_t log=0,Float_t minX=0,Float_t maxX=0,Float_t minY=0,Float_t maxY=0,Int_t scale=1)
{//online version of plotall()
  hup();
  plotall(histin,suffix,log,minX,maxX,minY,maxY,scale);

}

void hplotalllow(Char_t *histin, Char_t *suffix="", Int_t style=7, Int_t size=1, Int_t color=1)
{//online version of plotalllow()
  hup();
  plotalllow(histin,suffix,style,size,color);
}

/*---------------------------------------------------------------------------------------------*/
//Copied from util_new.cc for use with Scarlet at Argonne
void d1(Char_t *histname)
{
  hup();
  if( c1==0 ) mkCanvas();
  else c1->cd();
  draw(histname);
}

void d2(Char_t *histname)
{
  hup();
  if( c1==0 ) mkCanvas();
  else c1->cd();
  draw2(histname);
}

void zap(Char_t *histname)
{
  switch(whatis(histname)) {
  case 0:
    cout << "object " << histname << " was not found."<<endl;
    break;
  case 1:
    TH1F *hist=(TH1F*) gROOT->FindObject(histname);
    break;
  case 2:
    TH1D *hist=(TH1D*) gROOT->FindObject(histname);
    break;
  case 3:
    TH2F *hist=(TH2F*) gROOT->FindObject(histname);
    break;
  case 4:
    TH2D *hist=(TH2D*) gROOT->FindObject(histname);
    break;
  case 5:
    TH3F *hist=(TH3F *) gROOT->FindObject(histname);
    break;
  case 6:
    TProfile *hist=(TProfile *) gROOT->FindObject(histname);
    break;
  }
  hist->Reset();
}

void zapall()
{
  TList *hlist=gDirectory->GetList();
  TIterator *hiterator=hlist->MakeIterator();
  TH1 *htemp;
  while (htemp= (TH1 *) hiterator->Next()) {
    if (htemp->InheritsFrom("TH1")) htemp->Reset();
  }
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

void doprint2(Char_t * cnvname="cFit", Char_t * filename="print.ps", Char_t * pr="f1-phaser")
{//prints to F-150, adapted from plot_tools.cc
  const char cmd[255];
  char * fl = reinterpret_cast<char *>(filename);
  TCanvas *thecanvas=(TCanvas *)gROOT->FindObject(cnvname);
  thecanvas->SaveAs(filename);
  sprintf(cmd,"lpr -P %s %s",pr,fl);
  gSystem->Exec(cmd);
  cout<<fl<<" was sent to printer "<<pr<<endl;
}

void doprint3(Int_t printer_no=0)
{
  switch(printer_no){
  case 0://same as dorprint2()
    printf("Sending canvas \"cFit\" to Room F154 Phaser...\n");
    doprint2("cFit","print.ps","f1-phaser");
    break;
  case 1:
    printf("Sending canvas \"cFit\" to ATLAS Data Room HP LaserJet...\n");
    doprint2("cFit","print.ps","m0-clj3700");
    break;
  defualt:
    printf("No printer assigned to printer_n0=%d\n",printer_no);
    // break;
  }
}
