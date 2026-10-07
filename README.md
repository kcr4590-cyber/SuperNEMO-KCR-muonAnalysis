# SuperNEMO-KCR-muonAnalysis
Analysis code used to produce an angular distribution of the muon background at SuperNEMO. Made by Katherine Curtis-Rose.

# Introduction
In this repository you will find three folders:
- **__finalscripts__**, which contains the analysis scripts I used. Within this folder are two more folders, **__caloscripts__** and **__trackscripts__**, which contain the scripts used for the calorimeter clustering method and track reconstruction method respectively. These methods are explained in detail in my thesis and each script is explained below.
- **__bashscripts__**, which contains the bash scripts I used to execute my code on the CCIN2P3 server. A brief explanation of these scripts is given below.
**NOTE:** in some of these scripts files are opened/saved to my directory on the CCIN2P3 server. Please change the directories for your own use **before** running the code.
- **__csvfiles__**, which contains two csv files, **__OMGeometry.csv__** and **__timecalibration.csv__**.
- - **NOTE:** Move these files into the directory from which you are running the code. **__OMGeometry.csv__** can also be created by running the **__OMLookup.cpp__** script mentioned below. 

# caloscripts
These scripts were used to produce the angular distribution of muons by clustering calorimeter hits that are spatially and temporally close together and attempting to piece together the path of the muon. See **Section 3.5** of my dissertation for a more detailed explanation.
- **__OMLookup.cpp__**: RUN FIRST IF YOU DO NOT HAVE **__OMGeometry.csv__**. This script maps coordinates taken from the __gamma_om_x(y,z)__ branches in the ROOT tree in each data file to the corresponding OM and stores the mapped coordinates in a .csv file (**__OMGeometry.csv__**) to be used later. 
- **__caloanalysisPipeline.cpp__**: this script runs on each data file within the data directory used. When run, it applies a selection cut to the original data, calibrates the OMs in time (by reading the lookup table in **__timecalibration.csv__**), and clusters together OM hits that are close together in space and time.
- **__caloAnalysis3.cpp__**: 
