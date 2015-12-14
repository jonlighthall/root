{
  //------------------------------------
  //Load general macros and settings
  //------------------------------------
  //gROOT->LoadMacro("C:/root/macros/fit.cc");
  gROOT->LoadMacro("scripts/fit.cc");  
  //gROOT->LoadMacro("/home/emma/offline/lighthall/rootana/analysis/scripts/fit.cc");
  //gROOT->LoadMacro("scripts/online_plot_tools.cc");
  setdisplay(); //calls function in fit.cc
  settime(true);
  sethome();

  //------------------------------------
  //Load specfic macros and settings
  //------------------------------------
  //gROOT->LoadMacro("C:/root/macros/load_and_plot.cc");
  gROOT->LoadMacro("scripts/load_and_plot.cc");
  //gROOT->LoadMacro("/home/emma/offline/lighthall/rootana/analysis/scripts/load_and_plot.cc");
  gStyle->SetStatH(0.05);//reduces the size of the stat box

  //------------------------------------
  //Load HELIOS simulation  
  //------------------------------------
  //gROOT->LoadMacro("C:/root/macros/PlotSim.cxx");  
  
  //------------------------------------
  //Load PGAC simulation  
  //------------------------------------
  gROOT->LoadMacro("scripts/generator.cc");
  //gROOT->LoadMacro("/home/emma/offline/lighthall/rootana/analysis/scripts/generator.cc");
  //setbeam(0,0);
  //setangles(0,4.5,0,360);
  //setresm(0.6);
  //source();
}
