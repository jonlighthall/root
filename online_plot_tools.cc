


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


