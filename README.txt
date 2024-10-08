
This repository contains the data and analyses presented in 
"Teosinte populations exhibit weak local adaptation to their rhizosphere biota despite strong effects of biota source on teosinte fitness and traits"
by authors Anna M. O'Brien, Ruairidh J.H. Sawers, Jaime Gasca-Pineda, Ivan Baxter, Luis E. Eguiarte, Jeffrey Ross-Ibarra, and Sharon Y. Strauss.

Author Anna M. O'Brien wrote all analyses.

Script Teo_Collecting_Map_Figure.R is makes Figure S1.
It takes file "population locations.csv" and data downloaded from worldclim/bioclim as inputs. 
It is discouraged to re-share these data, but they are available freely on the distributer sites.

Script Teo-Biota_Ints_Main_Analysis.R performs the majority of analyses.
Please note sessionInfo.txt, which is a readout of the sessionInfo() function performed in R on the computer that ran these analyses.
This script takes all other .csv files in this repository as input, and produces all other .pdf files.

Input data files are described in R scripts at or near the line where they are read in. 

For permutations, pre-run outputs are provided, it is recommended to use these, rather than re-run if using a laptop or other time/RAM-restricted machine.
