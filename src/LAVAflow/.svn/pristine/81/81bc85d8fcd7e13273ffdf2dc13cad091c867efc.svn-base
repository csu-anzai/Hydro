
All of these analysis tools work by defining where the models reside on your local machine, giving each an index, defining where the analysis output will go, and what checkpoint file numbers are available for each model. These values are hard-coded. Then when running each tool, you can either hard-code the model index and checkpoint number and run the code, or you can specify the index of the model and the checkpoint number as command-line arguments. 

--------------------------------------
Setting up code on your local machine:
--------------------------------------

In CommonSNBinary.h:

This contains the model indices (the entries in the enum structure), the initial checkpoint numbers (not padded with zeros) and final checkpoint numbers. The inFolderLocation variable contains the folder containing these model folders. On my machine, for example, I have 'w7_ms38_l200-20070509' in the '/data1/Data/SNB_2007_Data/' folder. 

Change the inFolderLocation and outFolderLocation, and if you want to add another model, you'd have to add it the enum structure (anywhere before the NUMMODELS value, but ideally immediately before it) and add the initial and final checkpoint file numbers.

In CommonSNBinary.cxx:

This file defines the subfolder names (e.g. 'w7_ms38_l200-20070509' mentioned above). They correspond to the index in the enum structure. It should be obvious how to add another model here. The other functions here are just supporting functions used by some of the analysis tools. 

For any analysis tool here, you must also make subfolders in the model folder. This is where some output will be placed for each test. 

Example with all folder names:

.../w7_ms38_l200-20070509/columnDensities/
.../w7_ms38_l200-20070509/contamination/
.../w7_ms38_l200-20070509/eintTemp/
.../w7_ms38_l200-20070509/eintTime/
.../w7_ms38_l200-20070509/radialDistribution/
.../w7_ms38_l200-20070509/strippedMass/
.../w7_ms38_l200-20070509/velDist/

---------------------------------------
Setting up and running an example tool:
---------------------------------------

In StrippedMassSNBinary.cxx:

If you want to be able to just run the code without command-line arguments, hard-code the currentModel and checkNum values. Otherwise, run the code with command-line arguments like, e.g. "./StrippedMassSNBinary 2 40". This performs the stripped mass analysis on the 40th checkpoint file of model 'MS38_L200' (the third entry in the enum MODEL in CommonSNBinary.h). Note that the initial and final model numbers are constants, so if you're hard-coding, you can just use 'MS38_L200_F' as the final checkpoint file for MS38_L200.

The output for this tool is the simulation time and the total stripped mass in units of solar masses. It will also be output to a file in the model's folder under a 'strippedMass/' subfolder.

---------------------
Other tools included:
---------------------

The setup and execution of these tools is identical to that of the example above. Results are in the tools' respective subfolders and in command-line output. Note that some have hard-coded parameters that may need to be changed.

* ColDensSNBinary:

Calculates the column densities of the models over some user-defined (and hard-coded) range of angles and steps between those extremes. Most complex tool built for these models, so it may require some inspection of the code. 

You can change the ray origin location by hard-coding the defined ORIGINX, ORIGINY, and ORIGINZ values. Note that the ORIGINCHOICE value must then be USERDEF. Otherwise, the point of the highest density along the Y axis (X = 0) will be used. Can track either the companion mass, ejecta mass, or both.

* ContaminationSNBinary:

Finds how much the companion is contaminated by elements present only in the supernova. Oxygen, silicon, iron and nickel are tracked, but others could also be hard-coded in.

* EintTempSNBinary:

Builds a histogram of companion, ejecta, and unaffected ejecta (the bottom half of the domain) masses according to their temperatures. Minimum and maximum temperatures, as well as the number of bins, are hard-coded. Bins can have linear or log spacing (also hard-coded). 

* EintTimeSNBinary:

Calculates the total internal energy for both the (bound and unbound) companion material and the supernova ejecta. (The naming for this comes from the internal energy at a given time, i.e. checkpoint file).

* RadialDistributionSNBinary:

Uses much of the same code from the ColDensSNBinary tool to simply output the density at a given radius for each model. The resolution of the histogram output is dependent on the resolution of the model. 

* VelDistSNBinary:

Outputs velocity distribution for companion, supernova ejecta, and various isotopes. Number of bins, minimum and maximum velocities, and bin spacing are hard-coded but are easily modified. 