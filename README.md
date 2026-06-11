# Phytoxigene Automated Analysis

A streamlined pipeline to process qPCR results generated from the NOAA-GLERL HABs monitoring program using R.

---

## A. Scope and Application

This guide provides a walkthrough of the steps required to process qPCR results using R code. 
* **Prerequisites:** No previous knowledge of R is required.
* **Compatibility:** The current version of the code is limited to use with the following assays:
  * **Phytoxigene Total Cyanobacteria** (Catalogue No: 205-0050)
  * **Phytoxigene Toxin** (*mcyE* and *sxtA* targets) (Catalogue No: 205-0051)

---

## B. Summary of Method

The provided R code and accompanying directory structure process "raw" results from Phytoxigene Total Cyanobacteria and Toxin Gene assays, perform Quality Control (QC), and calculate final target concentrations. 

Both qPCR assays follow identical processing steps and can analyze samples run for Total Cyanobacteria, *mcyE*, and *sxtA* toxin genes simultaneously. Historical standards can also be integrated into the analysis to allow for seamless data processing across projects spanning multiple instrument runs.

---

## C. Definitions

| Term | Definition |
| :--- | :--- |
| **IAC** | Internal Amplification Control |
| **Template Volume** | Volume (µL) of gDNA added to each qPCR reaction (Phytoxigene protocol standard is 5 µL) |
| **Elution Volume** | Volume (µL) of eluant used to suspend the final gDNA extract |
| **Dilution Factor** | Multiplier used to calculate original concentration. Undiluted samples have a factor of 1. Dilutions mixed at a 1:1 ratio (1 part sample + 1 part diluent) have a factor of 2. For higher dilutions, use the total volume denominator (e.g., 1:5 dilution = factor of 5). |
| **CPR** | Copies Per Reaction |

---

## D. Dependencies

1. **Hardware:** QuantStudio 6 Real-Time PCR System
2. **Assays:** Phytoxigene Total Cyanobacteria Assay (Cat No: 205-0050) or Phytoxigene Toxin Assay (Cat No: 205-0051)
3. **Template File:** qPCR runs must be performed using the template: `WLEWeekly_Template_TCToxin.edt` *(available on GLERL PC, Lab 708)*
4. **Software:** R / RStudio (Download at [rstudio.com](https://www.rstudio.com))

---

## E. Preparation

### Export Run `.csv`
1. Locate the exported data (`.xls`/`.xlsx`) from the target QuantStudio 6 run. 
2. Click on the **Results** tab in the excel file.
3. Save **only** this tab as a `.csv` file into your unzipped `working_dir` directory. 
   * *Note: The file can be saved under any name but must use the `.csv` extension. A template named `run_file.csv` is included in the directory as a placeholder.*

### Setup `std_curve.csv` (Optional)
4. Copy and paste standard curve data from a previous run into `std_curve.csv`.
5. Rename `std_curve.csv` to include unique identifying information related to the run (e.g., for data run on August 6, 2022, rename to `std_curve_06Aug22.csv`).
6. Alternatively, copy a previously generated `std_curve.csv` file into the `working_dir`.

> ⚠️ **CRITICAL:** Do not add standard curve data from the *current* run being processed into `std_curve.csv`. That data will be automatically sourced from the exported run file.

---

## F. Running the Code

### 1. Define Variables in the Script
1. Open `analysis_code.R` in RStudio.
2. **Line 19:** Edit `'run_file.csv'` to match your results file.
3. **Line 23 (Optional):** Edit `'std_curve.csv'` located in the working directory to include all standards from previous runs.
4. **Line 32 (Optional):** Assign the **Template Volume** (5 µL is standard). If volumes vary across samples, replace `5` with `NA`. They can be added on a sample-by-sample basis later.
5. **Line 36:** Assign the **Elution Volume** (100 µL is standard). If volumes vary, replace `100` with `NA`.
6. **Line 41 (Optional):** Designate the names of gDNA extraction controls.
   * *Example:* If two controls named `extCNTRL1` and `extCNTRL2` were used, change `c('NA')` to `c('extCNTRL1', 'extCNTRL2')`.

### 2. Execute the Code
Run the code by pressing `Ctrl + Enter` (Windows), `Cmd + Enter` (Mac), or using the **Run** toolbar button.

### 3. Examine Output Files
After completion, look in the `output` folder of your `working_dir`. Open the file prefixed with **`Run_results_summary`**:

* **Rows 1–5:** Metadata pertaining to the QuantStudio 6 run.
* **Rows 7–8:** QC results for **Total Cyanobacteria** assay.
* **Rows 10–11:** QC results for **_mcyE_ Toxin** assay.
* **Rows 13–14:** QC results for **_sxtA_ Toxin** assay.
  * *Note: A period (`.`) indicates a passing measure. `FAIL` identifies at least one sample/set that did not meet specifications.*

#### Quality Control Test Matrix
| Test ID | Objective / Description | Target / Threshold |
| :--- | :--- | :--- |
| **Test 1** | Identifies standard curve outliers | Checks if reactions are within 1 CT of the average master curve point |
| **Test 2** | Performs modified Thompson Tau Test | Statistical outlier detection |
| **Test 3** | Checks standard curve upper bound | Checks if reaction is $< +3$ SD from the average master curve point |
| **Test 4** | Checks standard curve lower bound | Checks if reaction is $> -3$ SD from the average master curve point |
| **Test 5** | Identifies NTC amplification (contamination) | Confirms CT of No Template Control (NTC) is $> 36$ |
| **Test 6** | Identifies extraction control contamination | Confirms extraction control amplification is $\le$ Limit of Detection (100 CPR, point NA015) |
| **Test 7** | Identifies sample inhibition | Confirms sample IAC does not deviate $> 1.5$ CT from NTC IAC |
| **Test 8** | Identifies disagreement between replicates | Checks if duplicates differentiate by $\le 0.5$ CT |
| **Test 9** | Confirms sample doesn't exceed upper bounds | Confirms sample has a higher CT than the top calculated point (NA026) |
| **Test 10** | *Not available yet* Limit of Detection | Checks if reaction is $\ge 45$ copies/reaction |
| **Test 11** | *Not available yet* Limit of Quantification | Checks if reaction is $\ge 100$ copies/reaction |

* **Rows 16–26:** Standard curve parameters for all three assays (Total, *mcyE*, *sxtA*) generated from current and historical runs.
  * **Acceptance Criteria:** $R^2 > 0.985$ (Wells B15, B19) | Efficiency: $90\% < \text{Efficiency} < 110\%$ (Wells E15, B19)
    
* **Rows 28+:** Manually enter meta-parameters for each sample:
  * Column E (`dilution_factor`)
  * Column H (`volume_filtered_mL`)
  * Column F & G (If template or elution volumes were set as `NA` in step 1)

### Final Values & Special Cases
* **Final Concentration:** Generated in **Column J** (Standard Deviation in **Column K**) in units of **copies per mL** once all fields are complete.
* **Undetermined Samples:** If a sample has partial amplification among replicates (some "Undetermined", some with CT values), `Run_results_summary` calculates `CPR_mean` using only the amplified wells and adds a row noting the undetermined wells. Check the `Run_results` file for full clarity.

> 📝 **NOTE:** If all tests pass, the file prefixed `Run_results` contains the final accurate individual reaction data. **Do not use this file until all `FAIL` marks in the summary file are resolved.**

---

## Troubleshooting

When failures occur, use the specific test files located in `output/[total|mcyE|sxtA]/tests/` to trace and remedy the errors.

* **Tests 1–4:** Open `test1_2_3_4_results_[assay].csv`
  * Identify samples marked `FAIL` in columns M–P. Consider removing failed points from the standard curve.
* **Test 5:** Open `test5_result_[assay].csv`
  * Identify samples marked `FAIL` in column H. Consider rerunning the plate due to NTC contamination.
* **Test 6:** Open `test6_result_[assay].csv`
  * Identify samples marked `FAIL` in column H. Address contamination protocols during genomic DNA extraction.
* **Test 7:** *(Total assay only)* Open `test7_result_total.csv`
  * Identify samples marked `FAIL` in column H. Dilute failed samples and rerun both assays to avoid inhibition.
* **Test 8:** Open `test8_9_result_[assay].csv`
  * Identify samples marked `FAIL` in column L. Investigate pipetting consistency between replicates.
* **Test 9:** Open `test8_9_result_[assay].csv`
  * Identify samples marked `FAIL` in column M. Dilute samples and rerun to bring within standard range.

### Omitting Samples from Analysis
If you decide to omit an outlier or failed reaction:
1. Open your primary input file (`run_file.csv` or equivalent).
2. Starting at **Row 49**, change the value in **Column C** from `FALSE` to `TRUE` for the target row.
3. Save the file, re-run the R script, and re-examine the newly generated summary output.

---

## G. References

* Phytoxigene CyanoDTec Manual - Version 9, August 2019. [Download PDF Link](https://static1.squarespace.com/static/531043b0e4b013842a3999f0/t/5d788d085bd75417004e0916/1568181527263/CyanoDTec+Procedure+Ver9.pdf)
