{
  printf("Using GitHub version of rootlogon.C.\n");
  //Set macro location in .rootrc

  //------------------------------------
  //Historic settings
  //------------------------------------
  //gROOT->LoadMacro("/net/bavarians/utilities/util.cc");
  //gROOT->LoadMacro("/net/bavarians/utilities/dsp.cc");
  //gROOT->LoadMacro("/net/helios/B12/helios_b12_plottools.cc");

  //------------------------------------
  //Load general macros and settings
  //------------------------------------
  gROOT->LoadMacro("fit.cc");  
  //gROOT->LoadMacro("scripts/online_plot_tools.cc");
  setdisplay(); //calls function in fit.cc
  settime(true);
  sethome();

  //------------------------------------
  //Load specfic macros and settings
  //------------------------------------
  gROOT->LoadMacro("load_and_plot.cc");
  //gROOT->LoadMacro("temp.cc");
  gStyle->SetStatH(0.05);//reduces the size of the stat box

  //------------------------------------
  //Load HELIOS simulation  
  //------------------------------------
  //gROOT->LoadMacro("C:/root/macros/PlotSim.cxx");  
  
  //------------------------------------
  //Load PGAC simulation  
  //------------------------------------
  //gROOT->LoadMacro("generator.cc");
  //setbeam(0,0);
  //setangles(0,4.5,0,360);
  //setresm(0.6);
  //source();

  //------------------------------------
  //Load ANASEN trees
  //------------------------------------
  //printf("loading tree_structure.h...\n");
  //gROOT->LoadMacro("/home/lighthall/anasen/analysis_software/tree_structure.h");
  //gROOT->ProcessLine(".L /home/lighthall/anasen/analysis_software/tree_structure.h");
  Double_t radius=3.846284509;
  Double_t gold_pos=27.7495;
  TF1 *wfit = new TF1("wfit","[0]*sin(fmod(6-x+24,24)*TMath::TwoPi()/24+[1])+[2]",0,24);
  Double_t spacer[8]={27.7495,
		       22.8988,
		       16.8789,
		       15.6083,
		       12.4741,
		       7.4268,
		       61.5495,
		      -2.8505};
}
