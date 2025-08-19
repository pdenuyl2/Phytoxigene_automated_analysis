# Phytoxigene Automated Analysis

#### A.	SCOPE AND APPLICATION

1.	Walkthrough of steps to process qPCR results generated from the NOAA-GLERL HABs monitoring program using R code.  No previous knowledge of R is required. The current version of code is limited to use with the Phytoxigene Total Cyanobacteria (Catalogue No: 205-0050) and Toxin (*mcyE* and *sxtA* targets) (Catalogue No: 205-0051) assays.  

#### B.	SUMMARY OF METHOD
1.	R code and the accompanying directory/files are used to take the “raw” results from  Phytoxigene Total Cyanobacteria and Toxin Gene assays, perform QC, then provide final concentrations of each respective target. Both qPCR assays are processed following the same steps, with an ability to analysis samples run for Total Cyanobacteria, *mcyE*, and *sxtA* toxin genes in unison.  Previous standards can be included into the analysis to allow for easy processing of projects that span multiple runs. 

#### C.	DEFINITIONS 
&nbsp;&nbsp;**IAC:**	Internal amplification control  

&nbsp;&nbsp;**Template volume:**	Volume (µl) of gDNA added to each qPCR reaction (Phytoxigene protocol is 5 µl)  

&nbsp;&nbsp;**Elution volume:**	Volume (µl) of eluant used to suspend the final gDNA extract  

&nbsp;&nbsp;**Dilution factor:**	gDNA is diluted to prevent amplification above the standard curve or reduce inhibition. Dilutions of 1:1 are written as a dilution factor of 2, 1:5 is a dilution factor of 5, &nbsp;&nbsp;1:10 is a dilution factor of 10, etc.  

&nbsp;&nbsp;**CPR:**	Copies per reaction  

#### D.	 Dependencies
&nbsp;&nbsp;1.)	QuantStudio6 Real Time PCR System  
&nbsp;&nbsp;2.)	Phytoxigene Total Cyanobacteria Assay (Catalogue No: 205-0050) or Phytoxigene Toxin Assay (Catalogue No: 205-0051)  
&nbsp;&nbsp;3.)	qPCR run performed using template: 					WLEWeekly_Template_TCToxin.edt (avaialable on GLERL PC, Lab 708)  
&nbsp;&nbsp;4.)	R/RStudio - https://www.rstudio.com  

#### E.	PREPARATION
Export Run .csv

1.	Locate exported data (.xls/.xlsx) from QuantStudio6 run of interest					Data should have been automatically exported if the appropriate qPCR run template was used.

2.	Click on ‘results’ tab in the .xls/xlsx run file

3.	Save results tab as a .csv into the unzipped working_dir directory 					File can be save under any name, but must have “.csv” file extension
run_file.csv included in directory as a placeholder

Setup std_curve.csv (optional)
4.	Copy and paste either assay’s standard curve data from previous run data into std_curve.csv		

5.	Rename std_cuve.csv file to include unique information related to the run of interest	i.e. data run on August 6, 2022, rename to std_curve_06Aug22.csv

6.	Alternately, copy a previously used std_curve.csv file into the working_dir.  

Do not add standard curve data from the run being processed into std_curve.csv, that data will be sourced from the exported run file.  

#### F.	Run Code
Define variables in code
1.	Open analysis_code.R in Rstudio.

2.	Edit ‘run_file.csv’ (Line 19) to match the results file name from step 4.

3.	Edit ‘std_curve.csv’ (Line 23) to match the standard curve file name from steps 6-7.

4.	Assign template volume (Line 32) (5µL is standard for Phytoxigene protocols) 			  If template volume is not consistent across all samples, replace ‘5’ with ‘NA’. Template volumes can be added on a sample-by-sample basis in later steps of the protocol.

5.	Assign elution volume (Line 36) (100µL is generally standard) 			  		  If elution volume is not consistent across all samples, replace ‘100’ with ‘NA’.  Elution volumes can be added on a sample-by-sample basis in later steps of the protocol. 

6.	Designate names of gDNA extraction controls (Line 41) 			                           
Example:  If two gDNA extraction controls were included in run, named: extCNTRL1 and extCNTRL2 change c('NA') to c(‘extCNTRL1’, ‘extCNTRL2’)

7.	Run code by pressing the Ctrl+Enter (Windows), Command+Enter (Mac), or use the Run toolbar button. Examine output files. After running the code, there will be multiple output files generated in the “output” folder of the working_dir.

8.	Open the file with the prefix “Run_results_summary” located in the “output” folder.

Rows 1-5 : metadata pertaining to the QuantStudio6 run

Rows 7-8 : QC results for Total Cyanobacteria assay standard curve, NTC, gDNA extraction controls, and environmental samples. Periods “.” indicate that the measure passed, while “FAIL” identifies that there as at least one sample or set of samples that did not meet the specifications identified in the test.  

Rows 10-11 : QC results for *mcyE* Toxin assay standard curve, NTC, gDNA extraction controls, and environmental samples. Periods “.” indicate that the measure passed, while “FAIL” identifies that there as at least one sample or set of samples that did not meet the specifications identified in the test.  

Rows 13-14 : QC results for *sxtA* Toxin assay standard curve, NTC, gDNA extraction controls, and environmental samples. Periods “.” indicate that the measure passed, while “FAIL” identifies that there as at least one sample or set of samples that did not meet the specifications identified in the test.  

**Test 1:**	Checks every standard curve reaction to see if it is within 1 cycle threshold (CT) of the respective average master curve point
-Identify outliers in standard curve
 
**Test 2:**	Performs modified Thompson Tau Test

**Test 3:**	Checks every standard curve reaction to see if less than + 3 standard deviations from respective average master curve point	

**Test 4:**	Checks every standard curve reaction to see if greater than - 3 standard deviations from respective average master curve point	

**Test 5:**	Confirms CT of NTC is above 36   	
-Identify amplification (contamination) in NTC

**Test 6:**	Confirms amplification is not greater than the limit of detection (100 copies per reaction). Determined by standard curve point NA015.   	
-Identify substantial amplification (contamination) in gDNA extraction control

**Test 7:**	Confirms IAC from each sample does not deviate more than 1.5 cycle threshold from NTC IAC  
-A sample IAC over 1.5 cycle threshold from the NTC IAC suggests inhibition  
-A sample IAC under 1.5 cycle threshold from the NTC IAC does not indicate any issue  

**Test 8:**	Checks every sample run in duplicate to confirm measures don’t differentiate beyond 0.5 cycle threshold (CT)  
-Identify disagreement between replicates

**Test 9:**	Checks every reaction to confirm it has a higher CT than the top point (NA026) calculated (average) for the master standard curve 	
-Confirms sample does not exceed upper bounds of master curve

**Test 10:**	[only included in manual template] Checks that the reaction is equal to or above the limit of detection for the assay (>=45 copies/reaction)  

**Test 11:**	[only included in manual template] Checks that the reaction is equal to or above the limit of quantification for the assay (>=100 copies/reaction)  

Rows 16-18 : Total Cyanobacteria assay standard curve parameters generated from the standard curve reactions included in the processed run as well as previous runs (std_curve.csv; optional).  

Rows 20-22 : *mcyE* Toxin assay standard curve parameters generated from the standard curve reactions included in the processed run as well as previous runs (std_curve.csv; optional).  

Rows 24-26 : *sxtA* Toxin assay standard curve parameters generated from the standard curve reactions included in the processed run as well as previous runs (std_curve.csv; optional). 

##### Standard curve tests  
R-squared	>0.985	Wells B15, B19  
Efficiency	90< %Efficiency <110
	Wells E15, B19

Rows 28+ : Manually enter dilution factor (column E: “dilution_factor”) for each sample.  Manually enter the volume (mL) of sample water processed through the filter used for gDNA extraction (column H: “volume_filtered_mL”) for each sample. 

If template volume (column F) or elution volume (column G) was designated as ‘NA’ in the user defined variables of the code, manually enter values for each sample.

##### Final values
Once all fields (column C-H) have been filled out a final value will be generated in column J. The standard deviation between replicates (if run) is generated in column K.  Both are provided in units: copies per mL.  

##### Undetermined samples  
Samples that were designated with “Undetermined” CT values during the run will be included in the Run_results and Run_results_summary outputs; however, if a sample is run in replicate and some of the replicates are “Undetermined”, while others have CT values, the Run_results_summary will only calculate CPR_mean for those samples with CT values.  In this example, the Run_results_summary will have a row for the CPR_mean calculation and an additional row indicating that the sample was also “Undetermined”.  In this case, it may be useful to refer to the Run_results file for better clarity.  

Note:
If all tests complete without failures, the file with the prefix “Run_results” located in the “output” folder provides accurate individual reaction results.  Do not use until all problems (“FAIL”) indicated in the run summary file are addressed.  

Troubleshooting
	Not all runs will complete free from failures.  The point of this code is to provide a means to trace and address these problems when they occur. Each test result is detailed in the respective “total”, "mcyE”, and "sxtA"/tests folders.  When there is a failure in the summary file, examine the data further to determine whether individual samples can be omitted, or a new run needs to be completed.  

##### When failures occur: 

**Tests 1-4**  Locate: output/[total|mcyE|sxtA]/tests/test1_2_3_4_results_[total|mcyE|sxtA].csv  
-Identify samples designated as “FAIL” in columns M-P  
-Consider removing failed samples from standard curve  

**Test 5**  Locate: output/[total|mcyE|sxtA]/tests/test5_result_[total|mcyE|sxtA].csv  
-Identify samples designated as “FAIL” in column H  
-Consider rerunning entire plate, as no template control may show signs of contamination  

**Test 6**  
(Optional)	Locate: output/[total|mcyE|sxtA]/tests/test6_result_[total|mcyE|sxtA].csv  
-Identify samples designated as “FAIL” starting at column H   
-Consider taking action to address contamination during the genomic DNA extraction process  

**Test 7**
(Total assay only)	Locate: output/total/tests/test7_result_total.csv
-Identify samples designated as “FAIL” starting at column H 
-Dilute failed samples (both assays) and rerun  

**Test 8**	
Locate: output/[total|mcyE|sxtA]/tests/test8_9_result_[total|mcyE|sxtA].csv  
-Identify samples designated as “FAIL” in column L  
-Investigate reasons why replicates are not in agreement  

**Test 9**	
Locate: output/[total|mcyE|sxtA]/tests/test8_9_result_[total|mcyE|sxtA].csv  
-Identify samples designated as “FAIL” in column M  
-Dilute samples and rerun  

**Omitting samples:**
After identifying samples or reactions you want to omit, return to ‘run_file.csv’ or equivalent in the primary working_dir folder. 

Starting at row 49, column C of the results file provides the option to remove a sample from analysis by changing “FALSE” to “TRUE”.  

Save edited run file, rerun code, and reexamine output
	  	                      
#### G.	REFERENCES
1.	Phytoxigene CyanoDTec manual - version 9, August 2019 
	https://static1.squarespace.com/static/531043b0e4b013842a3999f0/t/5d788d085bd75417004e0916/1568181527263/CyanoDTec+Procedure+Ver9.pdf
	
 
