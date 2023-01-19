#Start with a fresh R environment 
rm(list=ls());if(is.null(dev.list()["RStudioGD"])){} else {dev.off(dev.list()["RStudioGD"])};cat("\014")

#Install dependencies (first use, remove "#")
#install.packages("tidyverse", type="binary")
#install.packages("ggplot2")
#install.packages("ggpmisc", type="binary")
#install.packages("dplyr")

#Load dependencies
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(dplyr)

######################START USER DEFINED VARIABLES######################

#Set file (.csv) - sourced from "results" tab
run_file <- 'results.csv'

###
#The following variables are optional and should only be used if the value is consistent across all samples being processed
#otherwise, leave blank or write 'NA'
###

#Set template volume (uL)
#designate NA if multiple template volumes were used for samples run
template_volume <- 5

#Set final volume from gDNA extraction elution
#designate NA if multiple elution volumes were used for samples run
elution_volume <- 100

#Designate gDNA Extraction controls
#example c("controlA","controlB")
#do not include NTC, designate "NA" if no extraction control was run
extraction_controls <- c("NA")

######################END USER DEFINED VARIABLES######################
#Set working directory to the location of this R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load run result file 
results <- read.csv(run_file)
#omit samples if designated in run result file
results <- results[results$X != "TRUE", ]  

#Isolate row 48 for use as a header
header_source <- c(as.data.frame(results)[48,])
header <- sub(" ", "_", header_source)
header <- sub("-", "_", header)

#Assign row 48 as header
colnames(results) <- header

#Remove the first set of rows that provide instrument and run details
run_file <- results[-c(1:48), ]

#Separate standards, NTC, and samples (unknown) - Total Cyanobacteria 16S rRNA gene Target
run_file_total_stds <- filter(run_file, Target_Name == "Total Cyanobacteria (16S rRNA gene)" & Task == "STANDARD")
run_file_total_ntc <- filter(run_file, Target_Name == "Total Cyanobacteria (16S rRNA gene)" & Task == "NTC" 
                                     | Target_Name == "Total Cyanobacteria (16S rRNA gene)" & Sample_Name == "NTC")
run_file_total_excntrl <- filter(run_file, Target_Name == "Total Cyanobacteria (16S rRNA gene)" & Sample_Name == extraction_controls)
run_file_total_unknown <- filter(run_file, Target_Name == "Total Cyanobacteria (16S rRNA gene)" & Task == "UNKNOWN" & !CT == "Undetermined")
run_file_total_unknown_undetermined <- filter(run_file, Target_Name == "Total Cyanobacteria (16S rRNA gene)" & Task == "UNKNOWN" & CT == "Undetermined")

#Separate standards, NTC, and samples (unknown) - Internal Amplification Control (IAC) Target 
run_file_iac_ntc <- filter(run_file, Target_Name == "Internal Amplification Control (IAC)" & Sample_Name == "NTC")
run_file_iac_stds_unknown <- filter(run_file, Target_Name == "Internal Amplification Control (IAC)" & Task == "UNKNOWN" & Sample_Name != "NTC")

#Separate standards, NTC, and samples (unknown) - mcyE Toxin gene Target
run_file_toxin_stds <- filter(run_file, Target_Name == "Microcystin/Nodularin (mcyE/ndaF)" & Task == "STANDARD")
run_file_toxin_ntc <- filter(run_file, Target_Name == "Microcystin/Nodularin (mcyE/ndaF)" & Task == "NTC" 
                                     | Target_Name == "Microcystin/Nodularin (mcyE/ndaF)" & Sample_Name == "NTC")
run_file_toxin_excntrl <- filter(run_file, Target_Name == "Microcystin/Nodularin (mcyE/ndaF)" & Sample_Name == extraction_controls)
run_file_toxin_unknown <- filter(run_file, Target_Name == "Microcystin/Nodularin (mcyE/ndaF)" & Task == "UNKNOWN" & !CT == "Undetermined")
run_file_toxin_unknown_undetermined <- filter(run_file, Target_Name == "Microcystin/Nodularin (mcyE/ndaF)" & Task == "UNKNOWN" & CT == "Undetermined")

#Load standard curve file 
std_curve <- read.csv(file = 'std_curve.csv', na.strings=c("NA",""))
std_curve[ , -which(names(std_curve) %in% c("Run_File"))]
colnames(std_curve) <- header

#Ensure consistent formatting between std_curve files
std_curve$Well<-as.character(std_curve$Well)
std_curve$Well_Position<-as.character(std_curve$Well_Position)
std_curve$Omit<-as.logical(std_curve$Omit)
std_curve$Sample_Name<-as.character(std_curve$Sample_Name)
std_curve$Target_Name<-as.character(std_curve$Target_Name)
std_curve$Task<-as.character(std_curve$Task)
std_curve$Reporter<-as.character(std_curve$Reporter)
std_curve$Quencher<-as.character(std_curve$Quencher)
std_curve$CT<-as.numeric(std_curve$CT)
std_curve$Ct_Mean<-as.numeric(std_curve$Ct_Mean)
std_curve$Ct_SD<-as.numeric(std_curve$Ct_SD)
std_curve$Quantity<-as.character(std_curve$Quantity)
std_curve$Quantity_Mean<-as.numeric(std_curve$Quantity_Mean)
std_curve$Quantity_SD<-as.numeric(std_curve$Quantity_SD)
std_curve$`Automatic_Ct Threshold`<-as.logical(std_curve$`Automatic_Ct Threshold`)
std_curve$Ct_Threshold<-as.numeric(std_curve$Ct_Threshold)
std_curve$Automatic_Baseline<-as.logical(std_curve$Automatic_Baseline)
std_curve$Baseline_Start<-as.numeric(std_curve$Baseline_Start)
std_curve$Baseline_End<-as.numeric(std_curve$Baseline_End)
std_curve$Comments<-as.character(std_curve$Comments)
std_curve$Y_Intercept<-as.numeric(std_curve$Y_Intercept)
std_curve$`R(superscript_2)`<-as.numeric(std_curve$`R(superscript_2)`)
std_curve$Slope<-as.numeric(std_curve$Slope)
std_curve$Efficiency<-as.numeric(std_curve$Efficiency)
std_curve$Amp_Score<-as.numeric(std_curve$Amp_Score)
std_curve$Cq_Conf<-as.numeric(std_curve$Cq_Conf)

run_file_total_stds$Well<-as.character(run_file_total_stds$Well)
run_file_total_stds$Well_Position<-as.character(run_file_total_stds$Well_Position)
run_file_total_stds$Omit<-as.logical(run_file_total_stds$Omit)
run_file_total_stds$Sample_Name<-as.character(run_file_total_stds$Sample_Name)
run_file_total_stds$Target_Name<-as.character(run_file_total_stds$Target_Name)
run_file_total_stds$Task<-as.character(run_file_total_stds$Task)
run_file_total_stds$Reporter<-as.character(run_file_total_stds$Reporter)
run_file_total_stds$Quencher<-as.character(run_file_total_stds$Quencher)
run_file_total_stds$CT<-as.numeric(run_file_total_stds$CT)
run_file_total_stds$Ct_Mean<-as.numeric(run_file_total_stds$Ct_Mean)
run_file_total_stds$Ct_SD<-as.numeric(run_file_total_stds$Ct_SD)
run_file_total_stds$Quantity<-as.character(run_file_total_stds$Quantity)
run_file_total_stds$Quantity_Mean<-as.numeric(run_file_total_stds$Quantity_Mean)
run_file_total_stds$Quantity_SD<-as.numeric(run_file_total_stds$Quantity_SD)
run_file_total_stds$`Automatic_Ct Threshold`<-as.logical(run_file_total_stds$`Automatic_Ct Threshold`)
run_file_total_stds$Ct_Threshold<-as.numeric(run_file_total_stds$Ct_Threshold)
run_file_total_stds$Automatic_Baseline<-as.logical(run_file_total_stds$Automatic_Baseline)
run_file_total_stds$Baseline_Start<-as.numeric(run_file_total_stds$Baseline_Start)
run_file_total_stds$Baseline_End<-as.numeric(run_file_total_stds$Baseline_End)
run_file_total_stds$Comments<-as.character(run_file_total_stds$Comments)
run_file_total_stds$Y_Intercept<-as.numeric(run_file_total_stds$Y_Intercept)
run_file_total_stds$`R(superscript_2)`<-as.numeric(run_file_total_stds$`R(superscript_2)`)
run_file_total_stds$Slope<-as.numeric(run_file_total_stds$Slope)
run_file_total_stds$Efficiency<-as.numeric(run_file_total_stds$Efficiency)
run_file_total_stds$Amp_Score<-as.numeric(run_file_total_stds$Amp_Score)
run_file_total_stds$Cq_Conf<-as.numeric(run_file_total_stds$Cq_Conf)
######################################################################
#Processing Stream - Total Cyanobacteria 16S rRNA gene Target - Start#
######################################################################
#Get all samples together to create a master total std curve 
std_curve_total <- filter(std_curve, Target_Name == "Total Cyanobacteria (16S rRNA gene)" & Task == "STANDARD")

ifelse(is.na(std_curve_total[1,9]), 
       std_curve_total_master <- run_file_total_stds,
       std_curve_total_master <- bind_rows(run_file_total_stds, std_curve_total))

#Add "Log Starting Quantity Machine" for standard curve generation
#Generate new corresponding values and column
std_curve_total_master$Log_starting_quantity_machine <- ifelse(grepl("NA026", std_curve_total_master$Sample_Name), 6, 
                                                               ifelse(grepl("NA015", std_curve_total_master$Sample_Name), 5,
                                                                      ifelse(grepl("NA014", std_curve_total_master$Sample_Name), 4,
                                                                             ifelse(grepl("NA013", std_curve_total_master$Sample_Name), 3,
                                                                                    ifelse(grepl("NA012", std_curve_total_master$Sample_Name), 2,"")))))

#Make sure all appropriate values are numeric (total std curve) 
std_curve_total_master$CT<-as.numeric(std_curve_total_master$CT)
std_curve_total_master$Ct_Mean<-as.numeric(std_curve_total_master$Ct_Mean)
std_curve_total_master$Ct_SD<-as.numeric(std_curve_total_master$Ct_SD)
std_curve_total_master$Y_Intercept<-as.numeric(std_curve_total_master$Y_Intercept)
std_curve_total_master$"R(superscript_2)"<-as.numeric(std_curve_total_master$"R(superscript_2)")
std_curve_total_master$Slope<-as.numeric(std_curve_total_master$Slope)
std_curve_total_master$Efficiency<-as.numeric(std_curve_total_master$Efficiency)
std_curve_total_master$Log_starting_quantity_machine <- as.numeric(std_curve_total_master$Log_starting_quantity_machine)

#Calculate and plot standard curve (total std curve)
std_curve_plot_total <- ggplot(std_curve_total_master,
                               aes(x=Log_starting_quantity_machine, y=CT)) +
  geom_point(size = 5) +
  geom_smooth(formula = y ~ x, method='lm', se=FALSE, color="black") +
  theme_bw() + 
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~~~")),
               formula = y~x, size = 6,
               parse = TRUE) +
  labs(x = "Log Starting Quantity (Copies Per Reaction)", y = "Ct") +
  scale_y_continuous(limits=c(0, 38), breaks = c(0,5,10,15,20,25,30,35))

ggsave(std_curve_plot_total, height = 5, width = 5, path = "output/Total", filename = "standard_curve_total.pdf", device = "pdf")

#Calculate and test r-squared of total std curve
std_curve_total_r2 <- tryCatch(summary(lm(CT ~ Log_starting_quantity_machine, data = std_curve_total_master))$r.squared, error=function(e){"NA"})
std_curve_total_r2 <- tryCatch(as.numeric(std_curve_total_r2), warning=function(e){std_curve_total_r2 <- "NA"})
std_curve_total_r2_test <- tryCatch(if(std_curve_total_r2>0.985) {
  "."
} else {
  "FAIL"
} , error=function(e){"NA"})
std_curve_total_r2_test <- ifelse(std_curve_total_r2 == "NA", std_curve_total_r2_test <- "NA", std_curve_total_r2_test)

#Calculate slope and intercept of total std curve
std_curve_total_slope <- tryCatch(summary(lm(CT ~ Log_starting_quantity_machine, data = std_curve_total_master))$coefficients[2], error=function(e){"NA"})
std_curve_total_intercept <- tryCatch(summary(lm(CT ~ Log_starting_quantity_machine, data = std_curve_total_master))$coefficients[1], error=function(e){"NA"})

#Calculate and test efficiency of total std curve
std_curve_total_efficiency <- tryCatch(((10^(-1/std_curve_total_slope))-1)*100, error=function(e){"NA"})
std_curve_total_efficiency_test <- tryCatch(if(std_curve_total_efficiency>90 & std_curve_total_efficiency<110) {
  "."
} else {
  "FAIL"
} , error=function(e){"NA"})
std_curve_total_efficiency_test <- ifelse(std_curve_total_efficiency == "NA", std_curve_total_efficiency_test <- "NA", std_curve_total_efficiency_test)

#Standard curve outlier check (tests 1-4)
#Create n_tau table
n_tau = matrix(c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,15, 16, 17,
                 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                 1.1511, 1.425, 1.5712, 1.6563, 1.711, 1.7491, 1.777, 1.7984, 1.8153, 1.829, 1.8403, 1.8498,
                 1.8579, 1.8649, 1.871, 1.8764, 1.8811, 1.8853, 1.8891, 1.8926, 1.8957, 1.8985, 1.9011, 
                 1.9035, 1.9057, 1.9078, 1.9096, 1.9114), ncol = 2, byrow=FALSE)
colnames(n_tau) <- c("n", "tau")

#Assign Tau to each sample based on (n) - Modified Thompson Tau - Test 2
std_curve_points <- c("NA026", "NA015", "NA014", "NA013", "NA012")
std_counts <- c()

for (i in std_curve_points) {
  std_counts[i] <- print(sum(std_curve_total_master$Sample_Name == i))
} 
std_counts <- as.matrix(std_counts)

std_curve_total_master$n <- ifelse(grepl("NA026", std_curve_total_master$Sample_Name), std_counts[1], 
                                   ifelse(grepl("NA015", std_curve_total_master$Sample_Name), std_counts[2],
                                          ifelse(grepl("NA014", std_curve_total_master$Sample_Name), std_counts[3],
                                                 ifelse(grepl("NA013", std_curve_total_master$Sample_Name), std_counts[4],
                                                        ifelse(grepl("NA012", std_curve_total_master$Sample_Name), std_counts[5],"")))))

std_curve_total_master$tau <- ifelse(grepl("3", std_curve_total_master$n), n_tau[1,2], 
                                     ifelse(grepl("4", std_curve_total_master$n), n_tau[2,2],
                                            ifelse(grepl("5", std_curve_total_master$n), n_tau[3,2],
                                                   ifelse(grepl("6", std_curve_total_master$n), n_tau[4,2],
                                                          ifelse(grepl("7", std_curve_total_master$n), n_tau[5,2],
                                                                 ifelse(grepl("8", std_curve_total_master$n), n_tau[6,2],
                                                                        ifelse(grepl("9", std_curve_total_master$n), n_tau[7,2],
                                                                               ifelse(grepl("10", std_curve_total_master$n), n_tau[8,2],
                                                                                      ifelse(grepl("11", std_curve_total_master$n), n_tau[9,2],
                                                                                             ifelse(grepl("12", std_curve_total_master$n), n_tau[10,2],
                                                                                                    ifelse(grepl("13", std_curve_total_master$n), n_tau[11,2],
                                                                                                           ifelse(grepl("14", std_curve_total_master$n), n_tau[12,2],
                                                                                                                  ifelse(grepl("15", std_curve_total_master$n), n_tau[13,2],
                                                                                                                         ifelse(grepl("16", std_curve_total_master$n), n_tau[14,2],
                                                                                                                                ifelse(grepl("17", std_curve_total_master$n), n_tau[15,2],
                                                                                                                                       ifelse(grepl("18", std_curve_total_master$n), n_tau[16,2],
                                                                                                                                              ifelse(grepl("19", std_curve_total_master$n), n_tau[17,2],
                                                                                                                                                     ifelse(grepl("20", std_curve_total_master$n), n_tau[18,2],
                                                                                                                                                            ifelse(grepl("21", std_curve_total_master$n), n_tau[19,2],
                                                                                                                                                                   ifelse(grepl("22", std_curve_total_master$n), n_tau[20,2],
                                                                                                                                                                          ifelse(grepl("23", std_curve_total_master$n), n_tau[21,2],
                                                                                                                                                                                 ifelse(grepl("24", std_curve_total_master$n), n_tau[22,2],
                                                                                                                                                                                        ifelse(grepl("25", std_curve_total_master$n), n_tau[23,2],
                                                                                                                                                                                               ifelse(grepl("26", std_curve_total_master$n), n_tau[24,2],
                                                                                                                                                                                                      ifelse(grepl("27", std_curve_total_master$n), n_tau[25,2],
                                                                                                                                                                                                             ifelse(grepl("28", std_curve_total_master$n), n_tau[26,2],
                                                                                                                                                                                                                    ifelse(grepl("29", std_curve_total_master$n), n_tau[27,2],
                                                                                                                                                                                                                          ifelse(grepl("30", std_curve_total_master$n), n_tau[28,2],""))))))))))))))))))))))))))))
###############################################################
# Total Cyanobacteria #
###############################################################
# TESTS #
#########
#Std. curve outlier checks (test 1-4) - Total Cyanobacteria
#Calculate avg. ct all runs
NA026_all_total_CT <- filter(std_curve_total_master, Sample_Name == "NA026") 
NA015_all_total_CT <- filter(std_curve_total_master, Sample_Name == "NA015") 
NA014_all_total_CT <- filter(std_curve_total_master, Sample_Name == "NA014") 
NA013_all_total_CT <- filter(std_curve_total_master, Sample_Name == "NA013") 
NA012_all_total_CT <- filter(std_curve_total_master, Sample_Name == "NA012") 

std_curve_total_master$CT_avg_allruns <- ifelse(grepl("NA026", std_curve_total_master$Sample_Name), mean(as.numeric(NA026_all_total_CT$CT)), 
                                                ifelse(grepl("NA015", std_curve_total_master$Sample_Name), mean(as.numeric(NA015_all_total_CT$CT)), 
                                                       ifelse(grepl("NA014", std_curve_total_master$Sample_Name), mean(as.numeric(NA014_all_total_CT$CT)),
                                                              ifelse(grepl("NA013", std_curve_total_master$Sample_Name), mean(as.numeric(NA013_all_total_CT$CT)),
                                                                     ifelse(grepl("NA012", std_curve_total_master$Sample_Name), mean(as.numeric(NA012_all_total_CT$CT)),"")))))

std_curve_total_master$Ct_SD_allruns <- ifelse(grepl("NA026", std_curve_total_master$Sample_Name), mean(as.numeric(NA026_all_total_CT$Ct_SD)), 
                                               ifelse(grepl("NA015", std_curve_total_master$Sample_Name), mean(as.numeric(NA015_all_total_CT$Ct_SD)), 
                                                      ifelse(grepl("NA014", std_curve_total_master$Sample_Name), mean(as.numeric(NA014_all_total_CT$Ct_SD)),
                                                             ifelse(grepl("NA013", std_curve_total_master$Sample_Name), mean(as.numeric(NA013_all_total_CT$Ct_SD)),
                                                                    ifelse(grepl("NA012", std_curve_total_master$Sample_Name), mean(as.numeric(NA012_all_total_CT$Ct_SD)),"")))))
#Test1 
#abs(CT_avg_allruns-CT)<1.00 
std_curve_total_master$CT_avg_allruns <- as.numeric(std_curve_total_master$CT_avg_allruns)
std_curve_total_master$test1 <- abs(std_curve_total_master$CT_avg_allruns-std_curve_total_master$CT)
tryCatch(for (i in 1:length(std_curve_total_master$test1)) {
  std_curve_total_master$test1[i] <- ifelse(std_curve_total_master$test1[i]<1.00, ".","Fail")
}, error=function(e){"NA"})

#Test2
#(CT_avg_allruns-CT)-(tau*Ct_SD_allruns)<0 
std_curve_total_master$Ct_SD_allruns <- as.numeric(std_curve_total_master$Ct_SD_allruns)
std_curve_total_master$test2 <- tryCatch(
  ((std_curve_total_master$CT_avg_allruns-std_curve_total_master$CT)-
     (std_curve_total_master$tau*std_curve_total_master$Ct_SD_allruns)),
  error = function(e)
    "-1000"
)
tryCatch(for (i in 1:length(std_curve_total_master$test2)) {
  std_curve_total_master$test2[i] <- ifelse(std_curve_total_master$test2[i]<0, ".","FAIL")
}, error=function(e){"NA"})

#Test3
#(CT_avg_allruns+(Ct_SD_allruns*3))>CT
std_curve_total_master$test3 <- tryCatch((std_curve_total_master$CT_avg_allruns+
                                            (std_curve_total_master$Ct_SD_allruns*3)), error=function(e){"NA"})
tryCatch(for (i in 1:length(std_curve_total_master$test3)) {
  std_curve_total_master$test3[i] <- ifelse(std_curve_total_master$test3[i]>std_curve_total_master$CT[i], ".","Fail")
}, error=function(e){"NA"})

#Test4
#(CT_avg_allruns-(Ct_SD_allruns*3))<CT
tryCatch(std_curve_total_master$test4 <- (std_curve_total_master$CT_avg_allruns-
                                            (std_curve_total_master$Ct_SD_allruns*3)), error=function(e){"NA"})

tryCatch(for (i in 1:length(std_curve_total_master$test4)) {
  std_curve_total_master$test4[i] <- ifelse(std_curve_total_master$test4[i]<std_curve_total_master$CT[i], ".","Fail")
}, error=function(e){"NA"})

test1_4_results_total <- std_curve_total_master %>% select(Well_Position, Sample_Name, Target_Name, Task, CT, Ct_Mean, Ct_SD, Quantity, n, tau, CT_avg_allruns, Ct_SD_allruns, test1, test2, test3, test4) 

write.csv(
  test1_4_results_total, file='output/total/tests/test1_2_3_4_results_total.csv', row.names=FALSE)

#NTC Amplification Check (Test 5)
tryCatch(for (i in 1:length(run_file_total_ntc$CT)) {
  run_file_total_ntc$test5[i] <- ifelse(run_file_total_ntc$CT[i]>36, ".","FAIL")
}, error=function(e){"NA"})

tryCatch(test5_result <- run_file_total_ntc %>% select(Well_Position, Sample_Name, Target_Name, Task, CT, Ct_Mean, Ct_SD, test5), error=function(e){"NA"})

tryCatch(write.csv(
  test5_result, file='output/total/tests/test5_result_total.csv', row.names=FALSE), error=function(e){"NA"})

#Extraction Blank Amplification Check (Test 6)
tryCatch(for (i in 1:length(run_file_total_excntrl$CT)) {
  run_file_total_excntrl$test6[i] <- ifelse(run_file_total_excntrl$CT[i]>as.numeric(NA012_all_total_CT$Ct_Mean[1]), ".","FAIL")
}, error=function(e){"NA"})

tryCatch(test6_result <- run_file_total_excntrl %>% select(Well_Position, Sample_Name, Target_Name, Task, CT, Ct_Mean, Ct_SD, test6), error=function(e){"NA"})

tryCatch(write.csv(
  test6_result, file='output/total/tests/test6_result_total.csv', row.names=FALSE), error=function(e){"NA"})

#IAC Check (Test 7)
#Separate standards, NTC, and samples (unknown) - Internal Amplification Control (IAC) Target
#Change "Undetermined" to a Ct of 42 to allow for absent Ct's to be identified
run_file_iac_stds_unknown_test7 <-  run_file_iac_stds_unknown
run_file_iac_stds_unknown_test7$CT[run_file_iac_stds_unknown_test7$CT == "Undetermined"] <- "42"

tryCatch(for(i in 1:length(run_file_iac_ntc$CT)) {                                                     # Head of for-loop                                            # Create new column
  run_file_iac_stds_unknown[ , ncol(run_file_iac_stds_unknown) + 1] <- ifelse(abs(as.numeric(run_file_iac_ntc$CT[i])-as.numeric(run_file_iac_stds_unknown_test7$CT))>1.5,"FAIL",".")  # Append new column
  
  colnames(run_file_iac_stds_unknown)[ncol(run_file_iac_stds_unknown)] <- paste0("test7_NTC", i)    # Rename column name
}, warning=function(e){""})
  
test7_result <- subset(run_file_iac_stds_unknown, select=-c(Well, Omit, Reporter, Quencher, Quantity, Quantity_Mean,
                                                            Quantity_SD, Ct_Threshold, Automatic_Baseline, Baseline_Start, Baseline_End,
                                                            Comments, Y_Intercept, Slope, Efficiency, Amp_Score))
test7_result <- test7_result[-c(0,8:9)]

write.csv(
  test7_result, file='output/total/tests/test7_result_total.csv', row.names=FALSE)
#######################
# SAMPLE CALCULATION #
#######################
# Total Cyanobacteria #

#Calculate copies per reaction (CPR) for each sample - no correction
tryCatch(run_file_total_unknown$CPR <- 10^((as.numeric(run_file_total_unknown$CT)-as.numeric(std_curve_total_intercept))/as.numeric(std_curve_total_slope)), warning=function(e){"NA"})

#Calculate new CPR mean without omitted samples
new_means <- run_file_total_unknown %>% group_by(Sample_Name, Target_Name)  %>% summarise(value=mean(CPR))
colnames(new_means)[colnames(new_means) == 'value'] <- 'CPR_mean'

#Recombine data frames
run_file_total_unknown <- merge(run_file_total_unknown, new_means, by=c("Sample_Name","Target_Name"))

#Calculate sample Stdev for each CPR calculated
tryCatch(sd_CPR <- run_file_total_unknown %>% group_by(Sample_Name) %>% summarise(value=sd(CPR)) , error=function(e){"NA"})
colnames(sd_CPR)[colnames(sd_CPR) == 'value'] <- 'Stdev_CPR'

run_file_total_unknown <- merge(run_file_total_unknown, sd_CPR, by="Sample_Name")
  
run_file_total_unknown$RSD <- abs(as.numeric(run_file_total_unknown$Stdev_CPR)/as.numeric(run_file_total_unknown$CPR))

#Difference between replicates <0.5 CT (Test 8)
tryCatch(for (i in 1:length(run_file_total_unknown$CT)) {
  run_file_total_unknown$test8[i] <- ifelse(abs(as.numeric(run_file_total_unknown$Ct_Mean[i])-as.numeric(run_file_total_unknown$CT[i]))<0.25, ".","FAIL")
}, error=function(e){"NA"})

#Samples within limit of highest std curve (Test 9)
tryCatch(for (i in 1:length(run_file_total_unknown$CT)) {
  run_file_total_unknown$test9[i] <- ifelse(as.numeric(run_file_total_unknown$CT[i])>as.numeric(NA026_all_total_CT$Ct_Mean[1]), ".","FAIL")
}, error=function(e){"NA"})
tryCatch(test8_9_result <- run_file_total_unknown %>% select(Well_Position, Sample_Name, Target_Name, Task, CT, Ct_Mean, Ct_SD, CPR, CPR_mean, Stdev_CPR, RSD, test8, test9), error=function(e){"NA"})

tryCatch(write.csv(
  test8_9_result, file='output/total/Tests/test8_9_result_total.csv', row.names=FALSE), error=function(e){"NA"})

######################################################################
#Processing Stream - Total Cyanobacteria 16S rRNA gene Target - End#
######################################################################

######################################################################
#Processing Stream - mcyE Toxin gene gene Target - Start#
######################################################################
#Ensure consistent formatting between std_curve files
run_file_toxin_stds$Well<-as.character(run_file_toxin_stds$Well)
run_file_toxin_stds$Well_Position<-as.character(run_file_toxin_stds$Well_Position)
run_file_toxin_stds$Omit<-as.logical(run_file_toxin_stds$Omit)
run_file_toxin_stds$Sample_Name<-as.character(run_file_toxin_stds$Sample_Name)
run_file_toxin_stds$Target_Name<-as.character(run_file_toxin_stds$Target_Name)
run_file_toxin_stds$Task<-as.character(run_file_toxin_stds$Task)
run_file_toxin_stds$Reporter<-as.character(run_file_toxin_stds$Reporter)
run_file_toxin_stds$Quencher<-as.character(run_file_toxin_stds$Quencher)
run_file_toxin_stds$CT<-as.numeric(run_file_toxin_stds$CT)
run_file_toxin_stds$Ct_Mean<-as.numeric(run_file_toxin_stds$Ct_Mean)
run_file_toxin_stds$Ct_SD<-as.numeric(run_file_toxin_stds$Ct_SD)
run_file_toxin_stds$Quantity<-as.character(run_file_toxin_stds$Quantity)
run_file_toxin_stds$Quantity_Mean<-as.numeric(run_file_toxin_stds$Quantity_Mean)
run_file_toxin_stds$Quantity_SD<-as.numeric(run_file_toxin_stds$Quantity_SD)
run_file_toxin_stds$`Automatic_Ct Threshold`<-as.logical(run_file_toxin_stds$`Automatic_Ct Threshold`)
run_file_toxin_stds$Ct_Threshold<-as.numeric(run_file_toxin_stds$Ct_Threshold)
run_file_toxin_stds$Automatic_Baseline<-as.logical(run_file_toxin_stds$Automatic_Baseline)
run_file_toxin_stds$Baseline_Start<-as.numeric(run_file_toxin_stds$Baseline_Start)
run_file_toxin_stds$Baseline_End<-as.numeric(run_file_toxin_stds$Baseline_End)
run_file_toxin_stds$Comments<-as.character(run_file_toxin_stds$Comments)
run_file_toxin_stds$Y_Intercept<-as.numeric(run_file_toxin_stds$Y_Intercept)
run_file_toxin_stds$`R(superscript_2)`<-as.numeric(run_file_toxin_stds$`R(superscript_2)`)
run_file_toxin_stds$Slope<-as.numeric(run_file_toxin_stds$Slope)
run_file_toxin_stds$Efficiency<-as.numeric(run_file_toxin_stds$Efficiency)
run_file_toxin_stds$Amp_Score<-as.numeric(run_file_toxin_stds$Amp_Score)
run_file_toxin_stds$Cq_Conf<-as.numeric(run_file_toxin_stds$Cq_Conf)

#Get all samples together to create a master toxin std curve ,
std_curve_toxin <- filter(std_curve, Target_Name == "Microcystin/Nodularin (mcyE/ndaF)" & Task == "STANDARD")

ifelse(is.na(std_curve_toxin[1,9]), 
       std_curve_toxin_master <- run_file_toxin_stds,
       std_curve_toxin_master <- bind_rows(run_file_toxin_stds, std_curve_toxin))

#Add "Log Starting Quantity Machine" for standard curve generation
#Generate new corresponding values and column
std_curve_toxin_master$Log_starting_quantity_machine <- ifelse(grepl("NA026", std_curve_toxin_master$Sample_Name), 6, 
                                                               ifelse(grepl("NA015", std_curve_toxin_master$Sample_Name), 5,
                                                                      ifelse(grepl("NA014", std_curve_toxin_master$Sample_Name), 4,
                                                                             ifelse(grepl("NA013", std_curve_toxin_master$Sample_Name), 3,
                                                                                    ifelse(grepl("NA012", std_curve_toxin_master$Sample_Name), 2,"")))))

#Make sure all appropriate values are numeric (toxin std curve) 
std_curve_toxin_master$CT<-as.numeric(std_curve_toxin_master$CT)
std_curve_toxin_master$Ct_Mean<-as.numeric(std_curve_toxin_master$Ct_Mean)
std_curve_toxin_master$Ct_SD<-as.numeric(std_curve_toxin_master$Ct_SD)
std_curve_toxin_master$Y_Intercept<-as.numeric(std_curve_toxin_master$Y_Intercept)
std_curve_toxin_master$"R(superscript_2)"<-as.numeric(std_curve_toxin_master$"R(superscript_2)")
std_curve_toxin_master$Slope<-as.numeric(std_curve_toxin_master$Slope)
std_curve_toxin_master$Efficiency<-as.numeric(std_curve_toxin_master$Efficiency)
std_curve_toxin_master$Log_starting_quantity_machine <- as.numeric(std_curve_toxin_master$Log_starting_quantity_machine)

#Calculate and plot standard curve (toxin std curve)
std_curve_plot_toxin <- ggplot(std_curve_toxin_master,
                               aes(x=Log_starting_quantity_machine, y=CT)) +
  geom_point(size = 5) +
  geom_smooth(formula = y ~ x, method='lm', se=FALSE, color="black") +
  theme_bw() + 
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~~~")),
               formula = y~x, size = 6,
               parse = TRUE) +
  labs(x = "Log Starting Quantity (Copies Per Reaction)", y = "Ct") +
  scale_y_continuous(limits=c(0, 38), breaks = c(0,5,10,15,20,25,30,35))

ggsave(std_curve_plot_toxin, height = 5, width = 5, path = "output/toxin", filename = "standard_curve_toxin.pdf", device = "pdf")

#Calculate and test r-squared of toxin std curve
std_curve_toxin_r2 <- tryCatch(summary(lm(CT ~ Log_starting_quantity_machine, data = std_curve_toxin_master))$r.squared, error=function(e){"NA"})
std_curve_toxin_r2 <- tryCatch(as.numeric(std_curve_toxin_r2), warning=function(e){std_curve_toxin_r2 <- "NA"})
std_curve_toxin_r2_test <- tryCatch(if(std_curve_toxin_r2>0.985) {
  "."
} else {
  "FAIL"
} , error=function(e){"NA"})
std_curve_toxin_r2_test <- ifelse(std_curve_toxin_r2 == "NA", std_curve_toxin_r2_test <- "NA", std_curve_toxin_r2_test)

#Calculate slope and intercept of toxin std curve
std_curve_toxin_slope <- tryCatch(summary(lm(CT ~ Log_starting_quantity_machine, data = std_curve_toxin_master))$coefficients[2], error=function(e){"NA"})
std_curve_toxin_intercept <- tryCatch(summary(lm(CT ~ Log_starting_quantity_machine, data = std_curve_toxin_master))$coefficients[1], error=function(e){"NA"})

#Calculate and test efficiency of toxin std curve
std_curve_toxin_efficiency <- tryCatch(((10^(-1/std_curve_toxin_slope))-1)*100, error=function(e){"NA"})
std_curve_toxin_efficiency_test <- tryCatch(if(std_curve_toxin_efficiency>90 & std_curve_toxin_efficiency<110) {
  "."
} else {
  "FAIL"
} , error=function(e){"NA"})
std_curve_toxin_efficiency_test <- ifelse(std_curve_toxin_efficiency == "NA", std_curve_toxin_efficiency_test <- "NA", std_curve_toxin_efficiency_test)

#Standard curve outlier check (tests 1-4)
#Create n_tau table
n_tau = matrix(c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,15, 16, 17,
                 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
                 1.1511, 1.425, 1.5712, 1.6563, 1.711, 1.7491, 1.777, 1.7984, 1.8153, 1.829, 1.8403, 1.8498,
                 1.8579, 1.8649, 1.871, 1.8764, 1.8811, 1.8853, 1.8891, 1.8926, 1.8957, 1.8985, 1.9011, 
                 1.9035, 1.9057, 1.9078, 1.9096, 1.9114), ncol = 2, byrow=FALSE)
colnames(n_tau) <- c("n", "tau")

#Assign Tau to each sample based on (n) - Modified Thompson Tau - Test 2
std_curve_points <- c("NA026", "NA015", "NA014", "NA013", "NA012")
std_counts <- c()

for (i in std_curve_points) {
  std_counts[i] <- print(sum(std_curve_toxin_master$Sample_Name == i))
} 
std_counts <- as.matrix(std_counts)

std_curve_toxin_master$n <- ifelse(grepl("NA026", std_curve_toxin_master$Sample_Name), std_counts[1], 
                                   ifelse(grepl("NA015", std_curve_toxin_master$Sample_Name), std_counts[2],
                                          ifelse(grepl("NA014", std_curve_toxin_master$Sample_Name), std_counts[3],
                                                 ifelse(grepl("NA013", std_curve_toxin_master$Sample_Name), std_counts[4],
                                                        ifelse(grepl("NA012", std_curve_toxin_master$Sample_Name), std_counts[5],"")))))

std_curve_toxin_master$tau <- ifelse(grepl("3", std_curve_toxin_master$n), n_tau[1,2], 
                                     ifelse(grepl("4", std_curve_toxin_master$n), n_tau[2,2],
                                            ifelse(grepl("5", std_curve_toxin_master$n), n_tau[3,2],
                                                   ifelse(grepl("6", std_curve_toxin_master$n), n_tau[4,2],
                                                          ifelse(grepl("7", std_curve_toxin_master$n), n_tau[5,2],
                                                                 ifelse(grepl("8", std_curve_toxin_master$n), n_tau[6,2],
                                                                        ifelse(grepl("9", std_curve_toxin_master$n), n_tau[7,2],
                                                                               ifelse(grepl("10", std_curve_toxin_master$n), n_tau[8,2],
                                                                                      ifelse(grepl("11", std_curve_toxin_master$n), n_tau[9,2],
                                                                                             ifelse(grepl("12", std_curve_toxin_master$n), n_tau[10,2],
                                                                                                    ifelse(grepl("13", std_curve_toxin_master$n), n_tau[11,2],
                                                                                                           ifelse(grepl("14", std_curve_toxin_master$n), n_tau[12,2],
                                                                                                                  ifelse(grepl("15", std_curve_toxin_master$n), n_tau[13,2],
                                                                                                                         ifelse(grepl("16", std_curve_toxin_master$n), n_tau[14,2],
                                                                                                                                ifelse(grepl("17", std_curve_toxin_master$n), n_tau[15,2],
                                                                                                                                       ifelse(grepl("18", std_curve_toxin_master$n), n_tau[16,2],
                                                                                                                                              ifelse(grepl("19", std_curve_toxin_master$n), n_tau[17,2],
                                                                                                                                                     ifelse(grepl("20", std_curve_toxin_master$n), n_tau[18,2],
                                                                                                                                                            ifelse(grepl("21", std_curve_toxin_master$n), n_tau[19,2],
                                                                                                                                                                   ifelse(grepl("22", std_curve_toxin_master$n), n_tau[20,2],
                                                                                                                                                                          ifelse(grepl("23", std_curve_toxin_master$n), n_tau[21,2],
                                                                                                                                                                                 ifelse(grepl("24", std_curve_toxin_master$n), n_tau[22,2],
                                                                                                                                                                                        ifelse(grepl("25", std_curve_toxin_master$n), n_tau[23,2],
                                                                                                                                                                                               ifelse(grepl("26", std_curve_toxin_master$n), n_tau[24,2],
                                                                                                                                                                                                      ifelse(grepl("27", std_curve_toxin_master$n), n_tau[25,2],
                                                                                                                                                                                                             ifelse(grepl("28", std_curve_toxin_master$n), n_tau[26,2],
                                                                                                                                                                                                                    ifelse(grepl("29", std_curve_toxin_master$n), n_tau[27,2],
                                                                                                                                                                                                                           ifelse(grepl("30", std_curve_toxin_master$n), n_tau[28,2],""))))))))))))))))))))))))))))



###############################################################
# mcyE Toxin gene #
###############################################################
# TESTS #
#########
#Std. curve outlier checks (test 1,2, 3, 4) - toxin Cyanobacteria
#Calculate avg. ct all runs
NA026_all_toxin_CT <- filter(std_curve_toxin_master, Sample_Name == "NA026") 
NA015_all_toxin_CT <- filter(std_curve_toxin_master, Sample_Name == "NA015") 
NA014_all_toxin_CT <- filter(std_curve_toxin_master, Sample_Name == "NA014") 
NA013_all_toxin_CT <- filter(std_curve_toxin_master, Sample_Name == "NA013") 
NA012_all_toxin_CT <- filter(std_curve_toxin_master, Sample_Name == "NA012") 

std_curve_toxin_master$CT_avg_allruns <- ifelse(grepl("NA026", std_curve_toxin_master$Sample_Name), mean(as.numeric(NA026_all_toxin_CT$CT)), 
                                                ifelse(grepl("NA015", std_curve_toxin_master$Sample_Name), mean(as.numeric(NA015_all_toxin_CT$CT)), 
                                                       ifelse(grepl("NA014", std_curve_toxin_master$Sample_Name), mean(as.numeric(NA014_all_toxin_CT$CT)),
                                                              ifelse(grepl("NA013", std_curve_toxin_master$Sample_Name), mean(as.numeric(NA013_all_toxin_CT$CT)),
                                                                     ifelse(grepl("NA012", std_curve_toxin_master$Sample_Name), mean(as.numeric(NA012_all_toxin_CT$CT)),"")))))

std_curve_toxin_master$Ct_SD_allruns <- ifelse(grepl("NA026", std_curve_toxin_master$Sample_Name), mean(as.numeric(NA026_all_toxin_CT$Ct_SD)), 
                                               ifelse(grepl("NA015", std_curve_toxin_master$Sample_Name), mean(as.numeric(NA015_all_toxin_CT$Ct_SD)), 
                                                      ifelse(grepl("NA014", std_curve_toxin_master$Sample_Name), mean(as.numeric(NA014_all_toxin_CT$Ct_SD)),
                                                             ifelse(grepl("NA013", std_curve_toxin_master$Sample_Name), mean(as.numeric(NA013_all_toxin_CT$Ct_SD)),
                                                                    ifelse(grepl("NA012", std_curve_toxin_master$Sample_Name), mean(as.numeric(NA012_all_toxin_CT$Ct_SD)),"")))))

#Test1
#abs(CT_avg_allruns-CT)<1.00 
std_curve_toxin_master$Ct_Mean <- as.numeric(std_curve_toxin_master$Ct_Mean)
std_curve_toxin_master$test1 <- abs(std_curve_toxin_master$Ct_Mean-std_curve_toxin_master$CT)
tryCatch(for (i in 1:length(std_curve_toxin_master$test1)) {
  std_curve_toxin_master$test1[i] <- ifelse(std_curve_toxin_master$test1[i]<1.00, ".","Fail")
}, error=function(e){"NA"})

#Test2
#(CT_avg_allruns-CT)-(tau*Ct_SD_allruns)<0 
std_curve_toxin_master$Ct_SD_allruns <- as.numeric(std_curve_toxin_master$Ct_SD_allruns)
std_curve_toxin_master$test2 <- tryCatch(
  ((std_curve_toxin_master$CT_avg_allruns-std_curve_toxin_master$CT)-
     (std_curve_toxin_master$tau*std_curve_toxin_master$Ct_SD_allruns)),
  error = function(e)
    "-1000"
)
tryCatch(for (i in 1:length(std_curve_toxin_master$test2)) {
  std_curve_toxin_master$test2[i] <- ifelse(std_curve_toxin_master$test2[i]<0, ".","FAIL")
}, error=function(e){"NA"})

#Test3
#(CT_avg_allruns+(Ct_SD_allruns*3))>CT
std_curve_toxin_master$test3 <- tryCatch((std_curve_toxin_master$Ct_Mean+
                                            (std_curve_toxin_master$Ct_SD*3)), error=function(e){"NA"})
tryCatch(for (i in 1:length(std_curve_toxin_master$test3)) {
  std_curve_toxin_master$test3[i] <- ifelse(std_curve_toxin_master$test3[i]>std_curve_toxin_master$CT[i], ".","Fail")
}, error=function(e){"NA"})

#Test4
#(CT_avg_allruns-(Ct_SD_allruns*3))<CT
tryCatch(std_curve_toxin_master$test4 <- (std_curve_toxin_master$Ct_Mean-
                                            (std_curve_toxin_master$Ct_SD*3)), error=function(e){"NA"})

tryCatch(for (i in 1:length(std_curve_toxin_master$test4)) {
  std_curve_toxin_master$test4[i] <- ifelse(std_curve_toxin_master$test4[i]<std_curve_toxin_master$CT[i], ".","Fail")
}, error=function(e){"NA"})

test1_4_results_toxin <- std_curve_toxin_master %>% select(Well_Position, Sample_Name, Target_Name, Task, CT, Ct_Mean, Ct_SD, Quantity, n, tau, CT_avg_allruns, Ct_SD_allruns, test1, test2, test3, test4) 
 
write.csv(
  test1_4_results_toxin, file='output/toxin/tests/test1_2_3_4_results_toxin.csv', row.names=FALSE)

#NTC amplification check (test 5)
tryCatch(for (i in 1:length(run_file_toxin_ntc$CT)) {
  run_file_toxin_ntc$test5[i] <- ifelse(run_file_toxin_ntc$CT[i]>36, ".","FAIL")
}, error=function(e){"NA"})

tryCatch(test5_result_toxin <- run_file_toxin_ntc %>% select(Well_Position, Sample_Name, Target_Name, Task, CT, Ct_Mean, Ct_SD, test5), error=function(e){"NA"})

tryCatch(write.csv(
  test5_result_toxin, file='output/toxin/tests/test5_result_toxin.csv', row.names=FALSE), error=function(e){"NA"})

#Extraction Blank Amplification Check (Test 6)
tryCatch(for (i in 1:length(run_file_toxin_excntrl$CT)) {
  run_file_toxin_excntrl$test6[i] <- ifelse(run_file_toxin_excntrl$CT[i]>as.numeric(NA012_all_toxin_CT$Ct_Mean[1]), ".","FAIL")
}, error=function(e){"NA"})

tryCatch(test6_result_toxin <- run_file_toxin_excntrl %>% select(Well_Position, Sample_Name, Target_Name, Task, CT, Ct_Mean, Ct_SD, test6), error=function(e){"NA"})

tryCatch(write.csv(
  test6_result_toxin, file='output/toxin/tests/test6_result_toxin.csv', row.names=FALSE), error=function(e){"NA"})

#IAC Check (Test 7) - NA for Toxin Assay

#######################
# SAMPLE CALCULATION #
#######################
# mcyE Toxin Gene #

#Calculate copies per reaction (CPR) for each sample - no correction
tryCatch(run_file_toxin_unknown$CPR <- 10^((as.numeric(run_file_toxin_unknown$CT)-as.numeric(std_curve_toxin_intercept))/as.numeric(std_curve_toxin_slope)), warning=function(e){"NA"})

#Calculate new CPR mean without omitted samples
tryCatch(new_means <- run_file_toxin_unknown %>% group_by(Sample_Name, Target_Name)  %>% summarise(value=mean(CPR)), error=function(e){"NA"})
colnames(new_means)[colnames(new_means) == 'value'] <- 'CPR_mean'

#Recombine data frames
run_file_toxin_unknown <- merge(run_file_toxin_unknown, new_means, by=c("Sample_Name","Target_Name"))

#Calculate sample mean for each CPR calculated
tryCatch(sd_CPR <- run_file_toxin_unknown %>% group_by(Sample_Name) %>% summarise(value=sd(CPR)), error=function(e){"NA"})
colnames(sd_CPR)[colnames(sd_CPR) == 'value'] <- 'Stdev_CPR'

run_file_toxin_unknown <- merge(run_file_toxin_unknown, sd_CPR, by="Sample_Name")

run_file_toxin_unknown$RSD <- abs(as.numeric(run_file_toxin_unknown$Stdev_CPR)/as.numeric(run_file_toxin_unknown$CPR))

#Difference between replicates <0.5 CT (Test 8)
tryCatch(for (i in 1:length(run_file_toxin_unknown$CT)) {
  run_file_toxin_unknown$test8[i] <- ifelse(abs(as.numeric(run_file_toxin_unknown$Ct_Mean[i])-as.numeric(run_file_toxin_unknown$CT[i]))<0.25, ".","FAIL")
}, error=function(e){"NA"})

#Samples within limit of highest std curve (Test 9)
tryCatch(for (i in 1:length(run_file_toxin_unknown$CT)) {
  run_file_toxin_unknown$test9[i] <- ifelse(as.numeric(run_file_toxin_unknown$CT[i])>as.numeric(NA026_all_toxin_CT$Ct_Mean[1]), ".","FAIL")
}, error=function(e){"NA"})
tryCatch(test8_9_result_toxin <- run_file_toxin_unknown %>% select(Well_Position, Sample_Name, Target_Name, Task, CT, Ct_Mean, Ct_SD, CPR, CPR_mean, Stdev_CPR, RSD, test8, test9), error=function(e){"NA"})

tryCatch(write.csv(
  test8_9_result_toxin, file='output/toxin/Tests/test8_9_result_toxin.csv', row.names=FALSE), error=function(e){"NA"})


######################################################################
#Processing Stream - mcyE Toxin gene gene Target - End#
######################################################################

###############################################################
# Print Summary - Total Cyanobacteria and mcyE Toxin Gene #
###############################################################
std_curve_summary_total_headers <- c("R-squared", "Slope", "Y-intercept", "Efficiency")
std_curve_summary_toxin_headers <- c("R-squared", "Slope", "Y-intercept", "Efficiency")
std_curve_summary_total_values <- c(std_curve_total_r2, std_curve_total_slope, std_curve_total_intercept, std_curve_total_efficiency)
std_curve_summary_toxin_values <- c(std_curve_toxin_r2, std_curve_toxin_slope, std_curve_toxin_intercept, std_curve_toxin_efficiency)
std_curve_total_summary <- rbind(std_curve_summary_total_headers, std_curve_summary_total_values)
std_curve_toxin_summary <- rbind(std_curve_summary_toxin_headers, std_curve_summary_toxin_values)

#Results summary: tests 1-4
std_curve_outlier_tests_headers <- c("test1", "test2", "test3", "test4")

std_curve_outlier_tests_results_total <- c(ifelse(any(std_curve_total_master$test1=="FAIL") == "TRUE", "FAIL", "."),
                                           ifelse(any(std_curve_total_master$test2=="FAIL") == "TRUE", "FAIL", "."),
                                           ifelse(any(std_curve_total_master$test3=="FAIL") == "TRUE", "FAIL", "."),
                                           ifelse(any(std_curve_total_master$test4=="FAIL") == "TRUE", "FAIL", "."))

std_curve_outlier_tests_results_toxin <- c(ifelse(any(std_curve_toxin_master$test1=="FAIL") == "TRUE", "FAIL", "."),
                                           ifelse(any(std_curve_toxin_master$test2=="FAIL") == "TRUE", "FAIL", "."),
                                           ifelse(any(std_curve_toxin_master$test3=="FAIL") == "TRUE", "FAIL", "."),
                                           ifelse(any(std_curve_toxin_master$test4=="FAIL") == "TRUE", "FAIL", "."))

std_curve_total_master_test_summary <- rbind(std_curve_outlier_tests_headers, std_curve_outlier_tests_results_total)
std_curve_toxin_master_test_summary <- rbind(std_curve_outlier_tests_headers, std_curve_outlier_tests_results_toxin)

#Results summary: test 5
run_file_ntc_test_header <- c("test5")

tryCatch(run_file_total_ntc_result <- run_file_total_ntc %>% select(Sample_Name, test5), error=function(e){"NA"})
tryCatch(run_file_toxin_ntc_result <- run_file_toxin_ntc %>% select(Sample_Name, test5), error=function(e){"NA"})

tryCatch(run_file_total_ntc_summary <- ifelse(any(run_file_total_ntc_result$test5=="FAIL") == "TRUE", "FAIL", "."), error=function(e){"NA"})
tryCatch(run_file_toxin_ntc_summary <- ifelse(any(run_file_toxin_ntc_result$test5=="FAIL") == "TRUE", "FAIL", "."), error=function(e){"NA"})

tryCatch(run_file_total_ntc_test_summary <- rbind(run_file_ntc_test_header, run_file_total_ntc_summary), error=function(e){"NA"})
tryCatch(run_file_toxin_ntc_test_summary <- rbind(run_file_ntc_test_header, run_file_toxin_ntc_summary), error=function(e){"NA"})

#Results summary: test 6
run_file_total_excntrl_test_header <- c("test6")
tryCatch(run_file_total_excntrl_summary <- ifelse(any(run_file_total_excntrl$test6=="FAIL") == "TRUE", "FAIL", "."), error=function(e){"NA"})
tryCatch(run_file_toxin_excntrl_summary <- ifelse(any(run_file_toxin_excntrl$test6=="FAIL") == "TRUE", "FAIL", "."), error=function(e){"NA"})

run_file_toxin_excntrl_test_header <- c("test6")
tryCatch(run_file_total_excntrl_summary <- ifelse(any(run_file_total_excntrl$test6=="FAIL") == "TRUE", "FAIL", "."), error=function(e){"NA"})
tryCatch(run_file_toxin_excntrl_summary <- ifelse(any(run_file_toxin_excntrl$test6=="FAIL") == "TRUE", "FAIL", "."), error=function(e){"NA"})

tryCatch(run_file_total_excntrl_test_summary <- rbind(run_file_total_excntrl_test_header, run_file_total_excntrl_summary), error=function(e){"NA"})
tryCatch(run_file_toxin_excntrl_test_summary <- rbind(run_file_toxin_excntrl_test_header, run_file_toxin_excntrl_summary), error=function(e){"NA"})

#Results summary: test 7
run_file_iac_stds_unknown_test_header <- c("test7")

run_file_iac_stds_unknown_result <- t(tail(t(run_file_iac_stds_unknown), n=length(run_file_iac_ntc$CT)))
run_file_iac_stds_unknown_result <- as.data.frame(run_file_iac_stds_unknown_result)
tryCatch(run_file_iac_stds_unknown_summary <- ifelse(any(run_file_iac_stds_unknown_result=="FAIL") == "TRUE", "FAIL", "."), error=function(e){"NA"})
tryCatch(run_file_iac_stds_unknown_test_summary <- rbind(run_file_iac_stds_unknown_test_header, run_file_iac_stds_unknown_summary), error=function(e){"NA"})

run_file_iac_stds_unknown_result_toxin <- "NA"    

run_file_iac_stds_unknown_test_summary_toxin <- rbind(run_file_iac_stds_unknown_test_header, run_file_iac_stds_unknown_result_toxin)

#Results summary: test 8 and 9
run_file_unknown_test_header <- c("test8", "test9")

#Change "NA" to "." for samples that are UNDETERMINED
run_file_total_unknown$test8 <- gsub('NA', '.', run_file_total_unknown$test8)
run_file_total_unknown$test9 <- gsub('NA', '.', run_file_total_unknown$test9)

run_file_total_unknown_result <- c(ifelse(any(run_file_total_unknown$test8=="FAIL") == "TRUE", "FAIL", "."),
                                   ifelse(any(run_file_total_unknown$test9=="FAIL") == "TRUE", "FAIL", "."))

run_file_total_unknown_test_summary <- rbind(run_file_unknown_test_header, run_file_total_unknown_result)

#Change "NA" to "." for samples that are UNDETERMINED
run_file_toxin_unknown$test8 <- gsub('NA', '.', run_file_toxin_unknown$test8)
run_file_toxin_unknown$test9 <- gsub('NA', '.', run_file_toxin_unknown$test9)
run_file_toxin_unknown_result <- c(ifelse(any(run_file_toxin_unknown$test8=="FAIL") == "TRUE", "FAIL", "."),
                                   ifelse(any(run_file_toxin_unknown$test9=="FAIL") == "TRUE", "FAIL", "."))

run_file_toxin_unknown_test_summary <- rbind(run_file_unknown_test_header, run_file_toxin_unknown_result)

#Format row names
row.names(std_curve_total_master_test_summary) <- c("Test", "result")
tryCatch(row.names(run_file_total_ntc_test_summary) <- c("Test", "result"), error=function(e){"NA"})
tryCatch(row.names(run_file_iac_stds_unknown_test_summary) <- c("Test", "result"), error=function(e){"NA"})
row.names(run_file_total_unknown_test_summary) <- c("Test", "result")

row.names(std_curve_toxin_master_test_summary) <- c("Test", "result")
tryCatch(row.names(run_file_toxin_ntc_test_summary) <- c("Test", "result"), error=function(e){"NA"})
row.names(run_file_total_unknown_test_summary) <- c("Test", "result")

#############################
#Summary of all test results#
#############################
tryCatch(all_test_summary_total <-cbind(std_curve_total_master_test_summary, run_file_total_ntc_test_summary, run_file_total_excntrl_test_summary, run_file_iac_stds_unknown_test_summary, run_file_total_unknown_test_summary), error=function(e){"NA"})
tryCatch(all_test_summary_toxin <-cbind(std_curve_toxin_master_test_summary, run_file_toxin_ntc_test_summary, run_file_toxin_excntrl_test_summary, run_file_iac_stds_unknown_test_summary_toxin, run_file_toxin_unknown_test_summary), error=function(e){"NA"})

#Create empty vector 
summary_output <- data.frame(matrix(ncol = 100, nrow = 100))
summary_output[is.na(summary_output)] <- ""

#Print test 1-9 summary - Total Cyanobacteria
tryCatch(summary_output[7:8,2:10] <- all_test_summary_total, error=function(e){"NA"})
summary_output[7,1]<- "Total Cyanobacteria Assay Test Results"

#Add filler test 1-9 summary if Total Cyanobacteria Assay not used
replacement_headers <- c("test1", "test2", "test3", "test4", "test5", "test6", "test7", "test8", "test9")
ifelse(summary_output[7,2]=="", summary_output[7,2:10] <- replacement_headers, print("NA"))
replacement_results <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
ifelse(summary_output[8,2]=="", summary_output[8,2:10] <- replacement_results, print("NA"))

#Print test 1-9 summary - mcyE Toxin Gene
tryCatch(summary_output[10:11,2:10] <- all_test_summary_toxin, error=function(e){"NA"})
summary_output[10,1] <- "mcyE Toxin Gene Assay Test Results"

#Add filler test 1-9 summary if mcyE Toxin Gene not used
replacement_headers <- c("test1", "test2", "test3", "test4", "test5", "test6", "test7", "test8", "test9")
ifelse(summary_output[10,2]=="", summary_output[10,2:10] <- replacement_headers, print("NA"))
replacement_results <- c("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
ifelse(summary_output[11,2]=="", summary_output[11,2:10] <- replacement_results, print("NA"))

#Summary of run/machine data
#Calibration expired?
calibration_state <- rbind(results[1,2], results[5,2], results[7,2], results[11,2])
calibration_state_result <- ifelse(any(calibration_state=="No") == "TRUE", ".", "FAIL")

summary_output[1,1]<- "Instrument Calibration"
summary_output[1,2] <- calibration_state_result

#User name
user_summary <- results[46,2]
summary_output[2,1] <- "Technician"
summary_output[2,2] <- user_summary

#Experiment name
name_summary <- results[36,2]
summary_output[3,1] <- "Experiment Name"
summary_output[3,2] <- name_summary

#Experiment file location
location_summary <- results[35,2]
summary_output[4,1] <- "Run File Location"
summary_output[4,2] <- location_summary

#Date created
date_summary <- results[32,2]
summary_output[5,1] <- "Run Start Date"
summary_output[5,2] <- date_summary

#Standard curve summary - Total Cyanobacteria
summary_output[13:14,2:5] <- std_curve_total_summary
summary_output[13,1] <- "Total Cyanobacteria Assay Std. Curve Metrics"
#Standard curve Pass/FAIL - Total Cyanobacteria
summary_output[15,2] <- std_curve_total_r2_test
summary_output[15,5] <- std_curve_total_efficiency_test

#Standard curve summary - mcyE Toxin Gene
summary_output[17:18,2:5] <- std_curve_toxin_summary
summary_output[17,1] <- "mcyE Toxin Gene Assay Std. Curve Metrics"
#Standard curve Pass/FAIL - mcyE Toxin Gene
summary_output[19,2] <- std_curve_toxin_r2_test
summary_output[19,5] <- std_curve_toxin_efficiency_test

#Results: sample concentrations and stats 
run_file_total_unknown_sample_results <- tryCatch(run_file_total_unknown %>% select(Well_Position, Sample_Name, Target_Name, CPR, CPR_mean, Stdev_CPR), error=function(e){"NA"})
tryCatch(if(run_file_total_unknown_sample_results[1]=="NA") {
  run_file_total_unknown_sample_results <- run_file_total_unknown
}, error=function(e){"NA"})
  
run_file_toxin_unknown_sample_results <- tryCatch(run_file_toxin_unknown %>% select(Well_Position, Sample_Name, Target_Name, CPR, CPR_mean, Stdev_CPR), error=function(e){"NA"})
tryCatch(if(run_file_toxin_unknown_sample_results[1]=="NA") {
  run_file_toxin_unknown_sample_results <- run_file_toxin_unknown
  }, error=function(e){"NA"})
  
run_file_unknown_sample_results_head <- c("Well_Position", "Sample_Name", "Target_Name", "CPR", "CPR_mean")
run_file_unknown_sample_results_head_trim <- c("Sample_Name", "Target_Name", "CPR_mean", "Stdev_CPR")

run_file_unknown_sample_results <- rbind(run_file_total_unknown_sample_results, run_file_toxin_unknown_sample_results)
run_file_unknown_sample_results %>% select(sort(names(.)))

run_file_unknown_sample_results_trim <- run_file_unknown_sample_results[!duplicated(run_file_unknown_sample_results$Sample_Name), ]
run_file_unknown_sample_results_trim <- run_file_unknown_sample_results_trim %>% select(Sample_Name, Target_Name, CPR_mean, Stdev_CPR)
run_file_unknown_sample_results_trim[,3:4] <- round(run_file_unknown_sample_results_trim[,3:4])

#Add Undetermined samples to summary (trim) file 
run_file_total_unknown_undetermined_sub_trim <- tryCatch(run_file_total_unknown_undetermined %>% select(Sample_Name, Target_Name, CT, Ct_SD), error=function(e){"NA"})
names(run_file_total_unknown_undetermined_sub_trim) <- c("Sample_Name", "Target_Name", "CPR_mean", "Stdev_CPR")
run_file_toxin_unknown_undetermined_sub_trim <- tryCatch(run_file_toxin_unknown_undetermined %>% select(Sample_Name, Target_Name, CT, Ct_SD), error=function(e){"NA"})
names(run_file_toxin_unknown_undetermined_sub_trim) <- c("Sample_Name", "Target_Name", "CPR_mean", "Stdev_CPR")

run_file_unknown_sample_results_trim <- rbind(run_file_unknown_sample_results_trim, run_file_total_unknown_undetermined_sub_trim, run_file_toxin_unknown_undetermined_sub_trim)

insert_length <- dim(run_file_unknown_sample_results_trim)[1]+21
insert_width <- dim(run_file_unknown_sample_results_trim)[2]

summary_output[22:insert_length,1:insert_width] <- run_file_unknown_sample_results_trim
summary_output[21, 1:insert_width] <- run_file_unknown_sample_results_head_trim

run_file_unknown_sample_results[,4:6] <- round(run_file_unknown_sample_results[,4:6])

#Add Undetermined samples to results file 
tryCatch(run_file_total_unknown_undetermined$blankVar <- "", error=function(e){"NA"})
run_file_total_unknown_undetermined_sub <- tryCatch(run_file_total_unknown_undetermined %>% select(Well_Position, Sample_Name, Target_Name, CT, blankVar, Ct_SD), error=function(e){"NA"})
tryCatch(names(run_file_total_unknown_undetermined_sub) <- c("Well_Position", "Sample_Name", "Target_Name", "CPR", "CPR_mean", "Stdev_CPR"), error=function(e){"NA"})
ifelse(run_file_total_unknown_undetermined_sub=="NA", run_file_total_unknown_undetermined_sub <- run_file_total_unknown_undetermined, print(""))

tryCatch(run_file_toxin_unknown_undetermined$blankVar <- "", error=function(e){"NA"})
run_file_toxin_unknown_undetermined_sub <- tryCatch(run_file_toxin_unknown_undetermined %>% select(Well_Position, Sample_Name, Target_Name, CT, blankVar, Ct_SD), error=function(e){"NA"})
tryCatch(names(run_file_toxin_unknown_undetermined_sub) <- c("Well_Position", "Sample_Name", "Target_Name", "CPR", "CPR_mean", "Stdev_CPR"), error=function(e){"NA"})
ifelse(run_file_toxin_unknown_undetermined_sub=="NA", run_file_toxin_unknown_undetermined_sub <- run_file_toxin_unknown_undetermined, print(""))

run_file_unknown_sample_results <- rbind(run_file_unknown_sample_results, run_file_total_unknown_undetermined_sub, run_file_toxin_unknown_undetermined_sub)
run_file_unknown_sample_results %>% select(sort(names(.)))

write.table(
  run_file_unknown_sample_results, file='output/Run_results.csv', row.names=FALSE, col.names=TRUE, sep=",")

#Include column for template volume (uL)
summary_output[21,6] <- "template_volume_uL"
tryCatch(summary_output[22:insert_length,6] <- template_volume, error=function(e){"NA"})

#Include column for elution volume (uL)
summary_output[21,7] <- "elution_volume_uL"
tryCatch(summary_output[22:insert_length,7] <- elution_volume, error=function(e){"NA"})

#Include column for manual additions (dilution factor, volume filtered)
summary_output[21,5] <- "dilution_factor"
summary_output[21,8] <- "volume_filtered_mL"

#Include column for final calculation in excel
summary_output[21,10] <- "FINAL Sample Value"
summary_output[21,11] <- "FINAL Sample Stdev"
summary_output[21,12] <- "Units"
summary_output[22:insert_length,12] <- "copies per mL"

formula_length<-nrow(run_file_unknown_sample_results_trim)
     
#Add formula(s) for automatic calculation of sample value and stdev in excel
for (i in 22:(21+formula_length)) {
summary_output[i,10] <- paste0("=((C",i,"*E",i,")/F",i,")*(G",i,"/H",i,")")
}

#Add formula(s) for automatic calculation of sample value and stdev in excel
for (i in 22:(21+formula_length)) {
  summary_output[i,11] <- paste0("=((D",i,"*E",i,")/F",i,")*(G",i,"/H",i,")")
}

write.table(
  summary_output, file='output/Run_results_summary.csv', row.names=FALSE, col.names=FALSE, sep=",")

#Version edits
#1.2 - added the ability of the code to identify NTC if the samples name is "NTC", regardless of assigned task
#    - corrected Test 7 to include IAC failures labeled "Undetermined" 

#1.3 -  "Undetermined" ct values were preventing final results from being generated.  Undetermined ct values are 
#now removed at the beginning of this process and then added back to the results prior to export.  
#    - "Well_Position" is added to Run_results.csv output

###END###









