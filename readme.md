# ROOT Macros

This repository contains ROOT macros developed for data analysis.

## Installation
Use either of the following commands to download a copy of the git repository.
`git clone https://github.com/jonlighthall/root.git` (HTTPS)
`git clone git@github.com:jonlighthall/root.git` (SSH)

## Loading
Set the macro path in '.rootrc' to include this directory.
````
# Path used to find macros.
# Paths are different for Unix and Windows. The example shows the defaults
# for all ROOT applications for either Unix or Windows.
Unix.*.Root.MacroPath:      .:~/root/scripts
WinNT.*.Root.MacroPath:     .;C:\Users\jonli\OneDrive\Documents\.cygwin_home\root.git\scripts;
````
Alternatively, the macros can be loaded via the `rootlogon.C`.

## Contents
## Procedure 
1. [fit.cc](#fit.cc)
2. [load_and_plot.cc](#load_and_plot.cc)
3. [generator.cc](#generator.cc)
3. [online_plot_tools.cc](#online_plot_tools.cc)
3. [rootlogon.C](#rootlogon.C)
3. [util_new.cc](#util_new.cc)

# fit.cc
The file `fit.cc` contains macros that are generally applicable to any data set.

# load_and_plot.cc
The file `load_and_plot.cc` contains macros for loading and plotting specific data sets.
They have been used by the author for various publications.
