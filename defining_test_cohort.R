#install.packages("dplyr", dependencies=TRUE)
#install.packages("ggplot2", dependencies=TRUE)
#install.packages("rmarkdown", dependencies=TRUE)
#install.packages("shiny", dependencies=TRUE)
#install.packages("readxl", dependencies=TRUE)
#install.packages("tidyverse", dependencies=TRUE)
#install.packages("magrittr", dependencies=TRUE)
#install.packages("ggpubr", dependencies=TRUE)
#install.packages("cowplot", dependencies=TRUE)
#install.packages("ggsci", dependencies=TRUE)
#install.packages("lubridate", dependencies=TRUE)
#install.packages("zoo", dependencies=TRUE)
#install.packages("formattable", dependencies=TRUE)
#install.packages("caret", dependencies=TRUE)
#install.packages("plotROC", dependencies=TRUE)
#install.packages("eeptools", dependencies=TRUE)
#install.packages("EpiStats", dependencies=TRUE)
#install.packages("pacman", dependencies=TRUE)
#install.packages("tableone", dependencies=TRUE)
#install.packages("xlsx")
pacman::p_load(
  rio,         
  here,         
  tidyverse,     
  stringr,      
  purrr,       
  gtsummary,    
  broom,    
  lmtest,   
  parameters, 
  see       
)

library(readxl)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(cowplot)
library(ggsci)
library(lubridate)
library(zoo)
library(formattable)
library(caret)
library(plotROC)
library(dplyr)
library(ggplot2)
library(shiny)
library(readr)
library(eeptools)
library(survival)
library(survminer)
library(EpiStats)
library(tableone)
library(dplyr)

#
setwd("/Volumes/Extreme SSD/Survival final")

## Explore CD4 count data
library(readr)
DR_Labs_HIV_update <- read_delim("DR_Labs_HIV_update.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
View(DR_Labs_HIV_update)
table(DR_Labs_HIV_update$common_concept)
HIV_Patients <- distinct(DR_Labs_HIV_update, study_id, .keep_all = FALSE)
DR_Labs_HIV_update$lab_assumed_taken_date <- as.Date(DR_Labs_HIV_update$lab_assumed_taken_date, "%Y-%m-%d")
Active_Patients <- subset(DR_Labs_HIV_update, lab_assumed_taken_date >= "2017-01-01" & lab_assumed_taken_date <= "2021-03-31")
Active_Patients <- distinct(Active_Patients, study_id, .keep_all = FALSE)
CD4_Counts <- subset(DR_Labs_HIV_update, lab_assumed_taken_date >= "2017-01-01" & lab_assumed_taken_date <= "2021-03-31")
CD4_Counts <- subset(CD4_Counts, common_concept == "Absolute CD4")
CD4_Counts <- select(CD4_Counts, study_id, lab_result_numerical)
CD4_Counts_Dist <- distinct(CD4_Counts, study_id, .keep_all = TRUE)

library(readr)
DR_Labs_HIV_CD4 <- read_delim("/Volumes/Extreme SSD/Survival final/DR_Labs_HIV_CD4.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
CD4_Data <- subset(DR_Labs_HIV_CD4, lab_assumed_taken_date >= "2017-01-01" & lab_assumed_taken_date <= "2021-03-31")
CD4_Data_Dist <- distinct(CD4_Data, study_id, .keep_all = TRUE)
CD4_Low <- subset(CD4_Data, lab_result_numerical < 200)
CD4_Low_Dist <- distinct(CD4_Low, study_id, .keep_all = TRUE)
summary(CD4_Data$lab_result_numerical)

##
###### Summary: There are 97 624 unique  patients in the overarching cohort, 72 102 of whom had 
#HIV related investigations within the study period. 48602 unique patients had a CD4 count done within the study period
#Of whom 13582 patients had a CD4 of less than 200
#Of all CD4 counts done in the study period, the median was 350,0 ; mean 400.6 ; 1Q 181 ; 3Q 559

###


#import CD4 data
library(readr)
DR_Labs_HIV_CD4 <- read_delim("DR_Labs_HIV_CD4.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
View(DR_Labs_HIV_CD4)

DR_Labs_HIV_CD4$lab_assumed_taken_date <- as.Date(DR_Labs_HIV_CD4$lab_assumed_taken_date, "%Y-%m-%d")
DR_Labs_HIV_CD4 <- subset(DR_Labs_HIV_CD4, lab_assumed_taken_date >= "2017-01-01")
AHD_CD4 <- subset(DR_Labs_HIV_CD4, lab_result_numerical < 200)
AHD_CD4 <- AHD_CD4[order(AHD_CD4$study_id, AHD_CD4$lab_assumed_taken_date),]
AHD_CD4 <- distinct(AHD_CD4, study_id, .keep_all = TRUE)
#check no enumerations after censor date
AHD_CD4$Censored_Date <- NA
AHD_CD4$Censored_Date <- "2021-03-31"
AHD_CD4$Censored_Date <- as.Date(AHD_CD4$Censored_Date, "%Y-%m-%d")
AHD_CD4$Time <- AHD_CD4$Censored_Date - AHD_CD4$lab_assumed_taken_date
#NO enumerations after censor date!

#enumerate cohort
Cohort <- distinct(AHD_CD4, study_id, .keep_all = FALSE)
#Add enumeration date
Cohort$Enumeration_date <- 
  AHD_CD4$lab_assumed_taken_date[match(Cohort$study_id, AHD_CD4$study_id)]
#Add Mortality data
library(readxl)
Deaths_Table_2021 <- read_excel("Deaths_Table_2021.xlsx")
View(Deaths_Table_2021)
Deaths_Table_2021$DATE_OF_DEATH <- as.Date(Deaths_Table_2021$DATE_OF_DEATH, "%Y-%m-%d")
Cohort$Censored_Date <- Deaths_Table_2021$DATE_OF_DEATH[match(Cohort$study_id, Deaths_Table_2021$STUDY_ID)]
#Add Censoring as binary
Cohort$Censored = 1
Cohort$Censored[is.na(Cohort$Censored_Date)] = 0
#make 31 March 2021 censoring date
Cohort$Censored_Date[is.na(Cohort$Censored_Date)] = "2021-03-31"
Cohort$Censored_Date <- as.Date(Cohort$Censored_Date, "%Y-%m-%d")
#Add linked column
Cohort$Linked <- NA
#Add time variable
Cohort$Enumeration_date <- as.Date(Cohort$Enumeration_date, "%Y-%m-%d")
Cohort$Time <- Cohort$Censored_Date - Cohort$Enumeration_date
#5 people have negative time periods
#Drop those with negative time periods
Cohort <- subset(Cohort, Time >= 0)

#######
#Add covariates 

##### CD4 counts data
Cohort$Enumeration_CD4 <- 
  AHD_CD4$lab_result_numerical[match(Cohort$study_id, AHD_CD4$study_id)]
Cohort$Enumeration_Facility <- 
  AHD_CD4$fac_name[match(Cohort$study_id, AHD_CD4$study_id)]
#add latest CD4
AHD_CD4 <- subset(DR_Labs_HIV_CD4, lab_result_numerical < 200)
AHD_CD4 <- AHD_CD4[rev(order(AHD_CD4$study_id, AHD_CD4$lab_assumed_taken_date)),]
AHD_CD4 <- distinct(AHD_CD4, study_id, .keep_all = TRUE)
Cohort$Latest_CD4 <- 
  AHD_CD4$lab_result_numerical[match(Cohort$study_id, AHD_CD4$study_id)]
## Cat CD4
#Make CD4 categories <50, 50-100, >100
Cohort$Enum_CD4_Cat <-cut(Cohort$Enumeration_CD4, c(-Inf,50,100, Inf))
Cohort$Enumeration_CD4_Cat <- NA
Cohort$Enumeration_CD4_Cat[Cohort$Enum_CD4_Cat=="(-Inf,50]"] <- "0-50"
Cohort$Enumeration_CD4_Cat[Cohort$Enum_CD4_Cat=="(50,100]"] <- "051-100" 
Cohort$Enumeration_CD4_Cat[Cohort$Enum_CD4_Cat=="(100, Inf]"] <- "101-199"
Cohort$E_CD4_Cat[Cohort$Enum_CD4_Cat=="(-Inf,50]"] <- 2
Cohort$E_CD4_Cat[Cohort$Enum_CD4_Cat=="(50,100]"] <- 1
Cohort$E_CD4_Cat[Cohort$Enum_CD4_Cat=="(100, Inf]"] <- 0

#Emueration CD4 done at hopsital or primary care
table(Cohort$Enumeration_Facility)
Hospital_Enumerations <- subset(Cohort, Enumeration_Facility == "Tygerberg Hospital" | Enumeration_Facility == "Khayelitsha Hospital"
                                | Enumeration_Facility == "Brooklyn Chest Hospital" | Enumeration_Facility == "New Somerset Hospital" 
                                | Enumeration_Facility == "Groote Schuur Hospital" | Enumeration_Facility == 	"Eerste River Hospital" 
                                | Enumeration_Facility == "Vredenburg Hospital" | Enumeration_Facility == "Victoria Hospital"
                                | Enumeration_Facility == "Lentegeur Hospital" | Enumeration_Facility == "Helderberg Hospital"
                                | Enumeration_Facility == "Mitchells Plain Hospital" | Enumeration_Facility == "False Bay Hospital"
                                | Enumeration_Facility == "Karl Bremer Hospital" | Enumeration_Facility == "Wesfleur Hospital"
                                | Enumeration_Facility == "Vredendal Hospital" | Enumeration_Facility == "Mowbray Maternity Hospital"
                                | Enumeration_Facility == "Swartland Hospital" | Enumeration_Facility == "Paarl Hospital"
                                | Enumeration_Facility == "Valkenberg Hospital" | Enumeration_Facility == "Red Cross War Memorial Children's Hospital"
                                | Enumeration_Facility == "Worcester Hospital" | Enumeration_Facility == "Hermanus Hospital"
                                | Enumeration_Facility == "SAMHS 2 Military Hospital" | Enumeration_Facility == "UCT Private Academic Hospital"
                                | Enumeration_Facility == "Stellenbosch Hospital" | Enumeration_Facility == "Knysna Hospital"
                                | Enumeration_Facility == "George Hospital" | Enumeration_Facility == "Mossel Bay Hospital"
                                | Enumeration_Facility == "DP Marais TB Hospital" | Enumeration_Facility == "Caledon Hospital"
                                | Enumeration_Facility == "Clanwilliam Hospital" | Enumeration_Facility == "Brewelskloof TB Hospital"
                                | Enumeration_Facility == "Ceres Hospital" | Enumeration_Facility == "Stikland Hospital"
                                | Enumeration_Facility == "Robertson Hospital" | Enumeration_Facility == "Riversdale Hospital"
                                | Enumeration_Facility == "Radie Kotze Hospital")
Hospital_Enumerations$Hospital_Enumerations <- 1
Cohort$Hospital_Enumerations <- Hospital_Enumerations$Hospital_Enumerations[match(Cohort$study_id, 
                                                                                  Hospital_Enumerations$study_id)]
Cohort$Hospital_Enumerations[is.na(Cohort$Hospital_Enumerations)] = 0
table(Cohort$Hospital_Enumerations)


############


#Viral load data
DR_Labs_HIV_VL <- read_delim("DR_Labs_HIV_VL.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE) 
DR_Labs_HIV_VL$Enumeration_Date <- Cohort$Enumeration_date[match(DR_Labs_HIV_VL$study_id, Cohort$study_id)]
DR_Labs_HIV_VL$Enumeration_Date <- as.Date(DR_Labs_HIV_VL$Enumeration_Date, "%Y-%m-%d")
DR_Labs_HIV_VL$lab_assumed_taken_date <- as.Date(DR_Labs_HIV_VL$lab_assumed_taken_date, "%Y-%m-%d")
DR_Labs_HIV_VL$Time <- DR_Labs_HIV_VL$Enumeration_Date - DR_Labs_HIV_VL$lab_assumed_taken_date
VL_Incident <- subset(DR_Labs_HIV_VL, Time < -60)
VL_Prevalent <- subset(DR_Labs_HIV_VL, Time  < 366 & Time > -61)
VL_Incident_limited <- select(VL_Incident, study_id, Time, lab_result_numerical, lab_result_assigned_numerical,
                              lab_result_text)
VL_Incident_limited$lab_result_numerical <- 
  replace(VL_Incident_limited$lab_result_numerical, VL_Incident_limited$lab_result_text == ">10000000", 750000)
VL_Incident_limited$lab_result_numerical <- 
  replace(VL_Incident_limited$lab_result_numerical, VL_Incident_limited$lab_result_numerical == "NULL", 0)
VL_Prevalent_limited <- select(VL_Prevalent, study_id, Time, lab_result_numerical, lab_result_assigned_numerical,
                               lab_result_text)
VL_Prevalent_limited$lab_result_numerical <- 
  replace(VL_Prevalent_limited$lab_result_numerical, VL_Prevalent_limited$lab_result_text == ">10000000", 750000)
VL_Prevalent_limited$lab_result_numerical <- 
  replace(VL_Prevalent_limited$lab_result_numerical, VL_Prevalent_limited$lab_result_numerical == "NULL", 0)

VL_Incident_limited <- VL_Incident_limited[order(VL_Incident_limited$study_id, VL_Incident_limited$Time),]
VL_Prevalent_limited <- VL_Prevalent_limited[order(VL_Prevalent_limited$study_id, VL_Prevalent_limited$Time),]
table(VL_Incident_limited$study_id)
table(VL_Prevalent_limited$study_id)
VL_I <- distinct(VL_Incident_limited, study_id, .keep_all = TRUE)
VL_P <- distinct(VL_Prevalent_limited, study_id, .keep_all = TRUE)
Cohort$VL_Enum <- VL_P$lab_result_numerical[match(Cohort$study_id, VL_P$study_id)]
Cohort$VL_Latest <- VL_I$lab_result_numerical[match(Cohort$study_id, VL_I$study_id)]

#Make enumeration VL cat
summary(Cohort$VL_Enum)
#Is character -> make numeric
Cohort$VL_Enum <- as.numeric(Cohort$VL_Enum)
summary(Cohort$VL_Enum)
Cohort$VL_E_cat <-cut(Cohort$VL_Enum, c(-Inf,99,1000, Inf))
summary(Cohort$VL_E_cat)
# First create the new field
Cohort$VL_E_Cat <- NA
# Then recode the old field into the new one for the specified rows
Cohort$VL_E_Cat[Cohort$VL_E_cat=="(-Inf,99]"] <- 0
Cohort$VL_E_Cat[Cohort$VL_E_cat=="(99,1e+03]"] <- 1
Cohort$VL_E_Cat[Cohort$VL_E_cat=="(1e+03, Inf]"] <- 2
table(Cohort$VL_E_Cat)

Cohort$VL_E_Not_Done[is.na(Cohort$VL_E_Cat)] = 1
Cohort$VL_E_Not_Done[is.na(Cohort$VL_E_Not_Done)] = 0


#Make Incident VL cat
summary(Cohort$VL_Latest)
#Is character -> make numeric
Cohort$VL_Latest <- as.numeric(Cohort$VL_Latest)
summary(Cohort$VL_Latest)
Cohort$VL_L_cat <-cut(Cohort$VL_Latest, c(-Inf,99,1000, Inf))
summary(Cohort$VL_L_cat)
# First create the new field
Cohort$VL_L_Cat <- NA
# Then recode the old field into the new one for the specified rows
Cohort$VL_L_Cat[Cohort$VL_L_cat=="(-Inf,99]"] <- 0
Cohort$VL_L_Cat[Cohort$VL_L_cat=="(99,1e+03]"] <- 1
Cohort$VL_L_Cat[Cohort$VL_L_cat=="(1e+03, Inf]"] <- 2
Cohort$VL_L_Cat[is.na(Cohort$VL_L_Cat)] = 3
table(Cohort$VL_L_Cat)

###### END of VL

######################

#Crypto Data
#import Crypto Data
Cryptococcus <- read_excel("Cryptococcus.xlsx")
table(Cryptococcus$common_method)
#Crypto positive (any test)
Crypto_Positive <- subset(Cryptococcus, common_result =="Positive")
table(Crypto_Positive$lab_test_description)
Crypto_Positive$lab_assumed_taken_date <- as.Date(Crypto_Positive$lab_assumed_taken_date, "%Y-%m-%d")
Crypto_Positive$Enumeration_date <- Cohort$Enumeration_date[match(Crypto_Positive$study_id, Cohort$study_id)]
Crypto_Positive$Enumeration_date <- as.Date(Crypto_Positive$Enumeration_date, "%Y-%m-%d")
Crypto_Positive$Time <- Crypto_Positive$Enumeration_date - Crypto_Positive$lab_assumed_taken_date
Crypto_Pos_Incident <- subset(Crypto_Positive, Time < -60)
Crypto_Pos_Prevalent <- subset(Crypto_Positive, Time  < 61 & Time > -61)
Crypto_Pos_Previous <- subset(Crypto_Positive, Time > 60)

#Incident cases
Crypto_Pos_Incident <- distinct(Crypto_Pos_Incident, study_id, .keep_all = TRUE)
#214 patients had microbiological evidence of Crpyto more than 60 days after enumeration
Crypto_Pos_Incident$Incident_Crypto <- 1
Cohort$Incident_Crypto <- Crypto_Pos_Incident$Incident_Crypto[match(Cohort$study_id, 
                                                                    Crypto_Pos_Incident$study_id)]
Cohort$Incident_Crypto[is.na(Cohort$Incident_Crypto)] = 0
table(Cohort$Incident_Crypto)
#Blood or CSF
Crypto_Pos_Incident <- subset(Crypto_Positive, Time < -60)
Crypto_I_Blood_Pos <- subset(Crypto_Pos_Incident, lab_specimen_type != "CSF")
Crypto_I_Blood_Pos <- subset(Crypto_I_Blood_Pos, lab_specimen_type != "BCSF")
Crypto_I_Blood_Pos <- subset(Crypto_I_Blood_Pos, lab_specimen_type != "FLUID CSF")
Crypto_I_Blood_Pos <- distinct(Crypto_I_Blood_Pos, study_id, .keep_all = TRUE)
Crypto_I_Blood_Pos$Blood_Pos_I <- 1
Crypto_Pos_Incident$Blood_Pos_I <- Crypto_I_Blood_Pos$Blood_Pos_I[match(Crypto_Pos_Incident$study_id,
                                                                        Crypto_I_Blood_Pos$study_id)]
Crypto_Pos_Incident$Blood_Pos_I[is.na(Crypto_Pos_Incident$Blood_Pos_I)] = 0
Cohort$Blood_Pos_I <- Crypto_Pos_Incident$Blood_Pos_I[match(Cohort$study_id,
                                                            Crypto_Pos_Incident$study_id)]
table(Cohort$Blood_Pos_I)
#197 people had serum evidence of incident crypto while 17 had only CSF evidence. 13363 NAs
#reflex testing
Crypto_Pos_Reflex <- subset(Crypto_Positive, common_method =="LFA Reflex Screening")
Crypto_Pos_R_I <- subset(Crypto_Pos_Reflex, Time < -60)
Crypto_Pos_R_I <- distinct(Crypto_Pos_R_I, study_id, .keep_all = TRUE)
Crypto_Pos_R_I$Reflex_Pos_I <- 1
Crypto_Pos_Incident$Reflex_Pos_I <- Crypto_Pos_R_I$Reflex_Pos_I[match(Crypto_Pos_Incident$study_id,
                                                                      Crypto_Pos_R_I$study_id)]
Crypto_Pos_Incident$Reflex_Pos_I[is.na(Crypto_Pos_Incident$Reflex_Pos_I)] = 0
Cohort$Reflex_Pos_I <- Crypto_Pos_Incident$Reflex_Pos_I[match(Cohort$study_id,
                                                              Crypto_Pos_Incident$study_id)]
table(Cohort$Reflex_Pos_I)
#168 had an incident reflex LFA pos


#Prevalent cases
Crypto_Pos_Prevalent <- distinct(Crypto_Pos_Prevalent, study_id, .keep_all = TRUE)
#229 patients had microbiological evidence of Crypto within 60 days before or after enumeration
Crypto_Pos_Prevalent$Prevalent_Crypto <- 1
Cohort$Prevalent_Crypto <- Crypto_Pos_Prevalent$Prevalent_Crypto[match(Cohort$study_id, 
                                                                       Crypto_Pos_Prevalent$study_id)]
Cohort$Prevalent_Crypto[is.na(Cohort$Prevalent_Crypto)] = 0
table(Cohort$Prevalent_Crypto)
#Blood or CSF
Crypto_Pos_Prevalent <- subset(Crypto_Positive, Time  < 61 & Time > -61)
Crypto_P_Blood_Pos <- subset(Crypto_Pos_Prevalent, lab_specimen_type != "CSF")
Crypto_P_Blood_Pos <- subset(Crypto_P_Blood_Pos, lab_specimen_type != "BCSF")
Crypto_P_Blood_Pos <- subset(Crypto_P_Blood_Pos, lab_specimen_type != "FLUID CSF")
Crypto_P_Blood_Pos <- distinct(Crypto_P_Blood_Pos, study_id, .keep_all = TRUE)
Crypto_P_Blood_Pos$Blood_Pos_P <- 1
Crypto_Pos_Prevalent$Blood_Pos_P <- Crypto_P_Blood_Pos$Blood_Pos_P[match(Crypto_Pos_Prevalent$study_id,
                                                                         Crypto_P_Blood_Pos$study_id)]
Crypto_Pos_Prevalent$Blood_Pos_P[is.na(Crypto_Pos_Prevalent$Blood_Pos_P)] = 0
Cohort$Blood_Pos_P <- Crypto_Pos_Prevalent$Blood_Pos_P[match(Cohort$study_id,
                                                             Crypto_Pos_Prevalent$study_id)]
table(Cohort$Blood_Pos_P)
#217 people had serum evidence of Prevalent crypto while 12 had only CSF evidence. 13348 NAs
#reflex testing
Crypto_Pos_R_P <- subset(Crypto_Pos_Reflex, Time  < 61 & Time > -61)
Crypto_Pos_R_P <- distinct(Crypto_Pos_R_P, study_id, .keep_all = TRUE)
Crypto_Pos_R_P$Reflex_Pos_P <- 1
Crypto_Pos_Prevalent$Reflex_Pos_P <- Crypto_Pos_R_P$Reflex_Pos_P[match(Crypto_Pos_Prevalent$study_id,
                                                                       Crypto_Pos_R_P$study_id)]
Crypto_Pos_Prevalent$Reflex_Pos_P[is.na(Crypto_Pos_Prevalent$Reflex_Pos_P)] = 0
Cohort$Reflex_Pos_P <- Crypto_Pos_Prevalent$Reflex_Pos_P[match(Cohort$study_id,
                                                               Crypto_Pos_Prevalent$study_id)]
table(Cohort$Reflex_Pos_P)
#204 had an prevalent reflex LFA pos

#Previous cases
Crypto_Pos_Previous <- distinct(Crypto_Pos_Previous, study_id, .keep_all = TRUE)
#122 patients had microbioligcal evidence of Crypto more than 60 days before enumeration
Crypto_Pos_Previous$Previous_Crypto <- 1
Cohort$Previous_Crypto <- Crypto_Pos_Previous$Previous_Crypto[match(Cohort$study_id, 
                                                                    Crypto_Pos_Previous$study_id)]
Cohort$Previous_Crypto[is.na(Cohort$Previous_Crypto)] = 0
table(Cohort$Previous_Crypto)
#Blood or CSF
Crypto_Pos_Previous <- subset(Crypto_Positive, Time > 60)
Crypto_Prev_Blood_Pos <- subset(Crypto_Pos_Previous, lab_specimen_type != "CSF")
Crypto_Prev_Blood_Pos <- subset(Crypto_Prev_Blood_Pos, lab_specimen_type != "BCSF")
Crypto_Prev_Blood_Pos <- subset(Crypto_Prev_Blood_Pos, lab_specimen_type != "FLUID CSF")
Crypto_Prev_Blood_Pos <- distinct(Crypto_Prev_Blood_Pos, study_id, .keep_all = TRUE)
Crypto_Prev_Blood_Pos$Blood_Pos_Prev <- 1
Crypto_Pos_Previous$Blood_Pos_Prev <- Crypto_Prev_Blood_Pos$Blood_Pos_Prev[match(Crypto_Pos_Previous$study_id,
                                                                                 Crypto_Prev_Blood_Pos$study_id)]
Crypto_Pos_Previous$Blood_Pos_Prev[is.na(Crypto_Pos_Previous$Blood_Pos_Prev)] = 0
Cohort$Blood_Pos_Prev <- Crypto_Pos_Previous$Blood_Pos_Prev[match(Cohort$study_id,
                                                                  Crypto_Pos_Previous$study_id)]
table(Cohort$Blood_Pos_Prev)
#81 people had serum evidence of Previous crypto while 41  had only CSF evidence. 13455 NAs
#reflex testing
Crypto_Pos_R_Prev <- subset(Crypto_Pos_Reflex, Time > 60)
Crypto_Pos_R_Prev <- distinct(Crypto_Pos_R_Prev, study_id, .keep_all = TRUE)
Crypto_Pos_R_Prev$Reflex_Pos_Prev <- 1
Crypto_Pos_Previous$Reflex_Pos_Prev <- Crypto_Pos_R_Prev$Reflex_Pos_Prev[match(Crypto_Pos_Previous$study_id,
                                                                               Crypto_Pos_R_Prev$study_id)]
Crypto_Pos_Previous$Reflex_Pos_Prev[is.na(Crypto_Pos_Previous$Reflex_Pos_Prev)] = 0
Cohort$Reflex_Pos_Prev <- Crypto_Pos_Previous$Reflex_Pos_Prev[match(Cohort$study_id,
                                                                    Crypto_Pos_Previous$study_id)]
table(Cohort$Reflex_Pos_Prev)
#29 had an Previous reflex LFA pos

#Reflex LFA overall
Crypto_Reflex <- subset(Cryptococcus, common_method =="LFA Reflex Screening")
#11648 observations
table(Crypto_Reflex$common_result)
#11141 Negative and 506 positive
Crypto_Reflex_Patients <- distinct(Crypto_Reflex, study_id, .keep_all = TRUE)
#7489 unique patients
Crypto_Reflex_Patients$Reflex_LFA_Done <- 1
Cohort$Reflex_LFA_Done <- Crypto_Reflex_Patients$Reflex_LFA_Done[match(Cohort$study_id, 
                                                                       Crypto_Reflex_Patients$study_id)]
Cohort$Reflex_LFA_Done[is.na(Cohort$Reflex_LFA_Done)] = 0
table(Cohort$Reflex_LFA_Done)
#6959 patients with a reflex LFA done

#Reflex LFA positive
Crypto_Reflex_Positive <- subset(Crypto_Reflex, common_result =="Positive")
Crypto_Reflex_Positive_Patients <- distinct(Crypto_Reflex_Positive, study_id, .keep_all = TRUE)
#349 unique patients with pos reflex LFA
Crypto_Reflex_Positive_Patients$Reflex_LFA_Positive <- 1
Cohort$Reflex_LFA_Positive <- Crypto_Reflex_Positive_Patients$Reflex_LFA_Positive[match(Cohort$study_id, 
                                                                                        Crypto_Reflex_Positive_Patients$study_id)]
table(Cohort$Reflex_LFA_Positive)
#335 patients with a reflex LFA Positive

##### END of crypto

## ART Data!
library(dplyr)

# 1. Load ART data
ART_Data <- read.csv("/Volumes/Extreme SSD/Survival final/ART_Data.csv")

# 2. Parse dates and bring in enumeration
ART_Data <- ART_Data %>%
  mutate(
    pharm_issue_date = as.Date(pharm_issue_date, "%Y-%m-%d"),
    Enumeration_date  = as.Date(Cohort$Enumeration_date[match(study_id, Cohort$study_id)])
  )

# 3. Compute first-ever ART date per patient
first_ART <- ART_Data %>%
  group_by(study_id) %>%
  summarise(
    first_ART_date = min(pharm_issue_date, na.rm = TRUE),
    .groups = "drop"
  )

# 4. Flag any ART dispense within 180 days prior to enumeration
recent_ART <- ART_Data %>%
  filter(
    pharm_issue_date >= Enumeration_date - 90,
    pharm_issue_date <  Enumeration_date
  ) %>%
  distinct(study_id) %>%
  mutate(recent_ART = TRUE)

# 5. Merge into Cohort and classify into Naive / On ART / LTFU
Cohort <- Cohort %>%
  left_join(first_ART,   by = "study_id") %>%
  left_join(recent_ART,  by = "study_id") %>%
  mutate(
    recent_ART = replace_na(recent_ART, FALSE),
    ART_Exp_C = case_when(
      # first ART < 180 days before enumeration → truly naive
      !is.na(first_ART_date) & first_ART_date >= Enumeration_date - 90       ~ "Naive",
      # started >180 days ago and still refilling → on ART
      !is.na(first_ART_date) & first_ART_date <  Enumeration_date - 90 & recent_ART ~ "On ART",
      # started >180 days ago but no refill in last 180 days → LTFU
      !is.na(first_ART_date) & first_ART_date <  Enumeration_date - 90 & !recent_ART ~ "LTFU",
      # never dispensed → naive
      is.na(first_ART_date)                                                     ~ "Naive"
    )
  )

# quick sanity check
table(Cohort$ART_Exp_C)

# 6. Define groups for downstream analyses
On_ART <- subset(Cohort, ART_Exp_C == "On ART")
Naive  <- subset(Cohort, ART_Exp_C == "Naive")
LTFU   <- subset(Cohort, ART_Exp_C == "LTFU")

# 7. Rebuild ART_AHD_Data for time-on-ART summaries
ART_AHD_Data <- ART_Data %>%
  filter(study_id %in% Cohort$study_id) %>%
  mutate(
    Time = as.numeric(Enumeration_date - pharm_issue_date)
  ) %>%
  arrange(study_id, pharm_issue_date)

# 8. Time on ART among those On ART
filtered_on_ART_ART_AHD_Data <- ART_AHD_Data %>%
  filter(study_id %in% On_ART$study_id)

oldest_ART_AHD_Data <- filtered_on_ART_ART_AHD_Data %>%
  group_by(study_id) %>%
  slice_max(Time, with_ties = FALSE) %>%
  ungroup()

combined_data <- oldest_ART_AHD_Data %>%
  select(study_id, pharm_issue_date) %>%
  left_join(On_ART %>% select(study_id, Enumeration_date, Censored), by = "study_id")

summary_data <- combined_data %>%
  mutate(Time_Difference = as.numeric(Enumeration_date - pharm_issue_date)) %>%
  group_by(Censored) %>%
  summarise(
    Average_Time = mean(Time_Difference, na.rm = TRUE),
    Median_Time  = median(Time_Difference, na.rm = TRUE),
    Min_Time     = min(Time_Difference, na.rm = TRUE),
    Max_Time     = max(Time_Difference, na.rm = TRUE),
    Q1_Time      = quantile(Time_Difference, 0.25, na.rm = TRUE),
    Q3_Time      = quantile(Time_Difference, 0.75, na.rm = TRUE),
    .groups      = "drop"
  )

print(summary_data)

# 9. Of those Naive, how many were subsequently initiated immediately (Time < 1)
Subs_ART <- ART_AHD_Data %>%
  filter(Time < 1) %>%
  distinct(study_id) %>%
  mutate(Subs_ART = 1)

Naive <- Naive %>%
  left_join(Subs_ART, by = "study_id") %>%
  mutate(Subs_ART = replace_na(Subs_ART, 0))

Naive_Subs <- subset(Naive, Subs_ART == 1)

Subs_ART_N <- ART_AHD_Data %>%
  filter(Time < 1, study_id %in% Naive_Subs$study_id) %>%
  group_by(study_id) %>%
  slice_max(Time, with_ties = FALSE) %>%
  ungroup() %>%
  select(study_id, pharm_issue_date)

Naive <- Naive %>%
  left_join(Subs_ART_N, by = "study_id") %>%
  rename(Date_First_Initiated = pharm_issue_date) %>%
  mutate(
    Time_To_First_ART = as.numeric(Date_First_Initiated - Enumeration_date)
  )

Cohort$Subs_Initiated      <- Naive$Subs_ART[match(Cohort$study_id, Naive$study_id)]
Cohort$Date_First_Initiated <- Naive$Date_First_Initiated[match(Cohort$study_id, Naive$study_id)]
Cohort$Time_To_First_ART    <- Naive$Time_To_First_ART[match(Cohort$study_id, Naive$study_id)]

# 10. Of those LTFU, how many re-engaged (Time < 1)
Subs_ART_L <- ART_AHD_Data %>%
  filter(Time < 1, study_id %in% LTFU$study_id) %>%
  distinct(study_id) %>%
  mutate(Subs_ART = 1)

LTFU <- LTFU %>%
  left_join(Subs_ART_L, by = "study_id") %>%
  mutate(Subs_ART = replace_na(Subs_ART, 0))

LTFU_Subs <- subset(LTFU, Subs_ART == 1)

Subs_ART_L_dates <- ART_AHD_Data %>%
  filter(Time < 1, study_id %in% LTFU_Subs$study_id) %>%
  group_by(study_id) %>%
  slice_max(Time, with_ties = FALSE) %>%
  ungroup() %>%
  select(study_id, pharm_issue_date)

LTFU <- LTFU %>%
  left_join(Subs_ART_L_dates, by = "study_id") %>%
  rename(Date_Reinitiated = pharm_issue_date) %>%
  mutate(
    Time_Reinitiation = as.numeric(Date_Reinitiated - Enumeration_date)
  )

Cohort$Reinitiated       <- LTFU$Subs_ART[match(Cohort$study_id, LTFU$study_id)]
Cohort$Date_Reinitiated  <- LTFU$Date_Reinitiated[match(Cohort$study_id, LTFU$study_id)]
Cohort$Time_Reinitiation <- LTFU$Time_Reinitiation[match(Cohort$study_id, LTFU$study_id)]
table(Cohort$Reinitiated)

#Add bactrim
##Bactrim data
#Co_tri_Data <- subset(DR_Pharmarcy, pharm_description =="sulfamethoxazole and trimethoprim")
#write_csv(Co_tri_Data,"/Volumes/Extreme SSD/Survival final/Co_tri_Data.csv")
C#o_tri_patients <- distinct(Co_tri_Data, study_id, .keep_all = TRUE)
##47220 unique patients have recieved co-trimoxazole
#Co_tri_patients$Recieved_Co_trimoxazole <- 1
#Cohort$Recieved_Co_trimoxazole <- Co_tri_patients$Recieved_Co_trimoxazole[match(Cohort$study_id, 
#                                                                                Co_tri_patients$study_id)]
#Cohort$Recieved_Co_trimoxazole[is.na(Cohort$Recieved_Co_trimoxazole)] = 0
#table(Cohort$Recieved_Co_trimoxazole)
#9447 patients have recieved co-trimoxazole before

#Add Fluconazole
#table(DR_Pharmarcy$pharm_description)
#flu_data <- subset(DR_Pharmarcy, pharm_description =="fluconazole")
#write_csv(flu_data,"/Volumes/Extreme SSD/Survival final/Cohort_fluconazole_data.csv")


##

### Demographics

#####Demographics
PAT_Site_B_youth <- read_excel("PAT Site B youth.xlsx")
PAT_Site_C_Youth <- read_excel("PAT Site C Youth.xlsx")
PAT_Town_2 <- read_excel("PAT Town 2.xlsx")
PAT_Zakhele <- read_excel("PAT Zakhele.xlsx")
PAT_Khuyasa_CDC <- read_excel("PAT_Khuyasa CDC.xlsx")
PAT_Khuyasa_males <- read_excel("PAT_Khuyasa males.xlsx")
PAT_Luvuyo <- read_excel("PAT_Luvuyo.xlsx")
PAT_Mathew_Goniwe <- read_excel("PAT_Mathew Goniwe.xlsx")
PAT_Mayenzeke <- read_excel("PAT_Mayenzeke.xlsx")
PAT_Michael_M <- read_excel("PAT_Michael M.xlsx")
PAT_Nolungile <- read_excel("PAT_Nolungile.xlsx")
PAT_Site_B_Males <- read_excel("PAT_Site B Males.xlsx")
PAT_Site_B_MOU <- read_excel("PAT_Site B MOU.xlsx")
PAT_Site_B_PHCIS <- read_excel("PAT_Site B PHCIS.xlsx")
PAT_Site_B_PHDC <- read_excel("PAT_Site B PHDC.xlsx")

Demographics <- rbind(PAT_Site_B_youth, PAT_Site_C_Youth) 
Demographics <- rbind(Demographics, PAT_Town_2)
Demographics <- rbind(Demographics, PAT_Zakhele)
Demographics <- rbind(Demographics, PAT_Khuyasa_CDC)
Demographics <- rbind(Demographics, PAT_Khuyasa_males)
Demographics <- rbind(Demographics, PAT_Luvuyo)
Demographics <- rbind(Demographics, PAT_Mathew_Goniwe)
Demographics <- rbind(Demographics, PAT_Mayenzeke)
PAT_Michael_M <- select(PAT_Michael_M, -c(LAST_CONTACT_DMY))
Demographics <- rbind(Demographics, PAT_Michael_M)
Demographics <- rbind(Demographics, PAT_Site_B_Males)
PAT_Nolungile_D <- select(PAT_Nolungile, PATIENT, COHORT, FACILITY, BIRTH_DMY, GENDER)
PAT_Site_B_MOU_D <- select(PAT_Site_B_MOU, PATIENT, COHORT, FACILITY, BIRTH_DMY, GENDER)
PAT_Site_B_PHCIS_D <- select(PAT_Site_B_PHCIS, PATIENT, COHORT, FACILITY, BIRTH_DMY, GENDER)
PAT_Site_B_PHDC_D <- select(PAT_Site_B_PHDC, PATIENT, COHORT, FACILITY, BIRTH_DMY, GENDER)
D <- rbind(PAT_Nolungile_D, PAT_Site_B_MOU_D, PAT_Site_B_PHCIS_D, PAT_Site_B_PHDC_D)
Demographics <- select(Demographics, PATIENT, COHORT, FACILITY, BIRTH_DMY, GENDER)
Demographics <- rbind(Demographics, D)
Demographics$BIRTH_DMY <- as.Date(Demographics$BIRTH_DMY, "%Y-%m-%d")
Cohort$DOB <- Demographics$BIRTH_DMY[match(Cohort$study_id, Demographics$PATIENT)]
Cohort$DOB <- as.Date(Cohort$DOB, "%Y-%m-%d")
Cohort$Enumeration_date <- as.Date(Cohort$Enumeration_date, "%Y-%m-%d")
Cohort$Sex <- Demographics$GENDER[match(Cohort$study_id, Demographics$PATIENT)]
Cohort$DOB[is.na(Cohort$DOB)] = "1900-01-01"
Cohort <- Cohort[Cohort$DOB < "2017-01-01",]
Cohort$Age <- NA
Cohort$Age <- age_calc(Cohort$DOB, enddate = Cohort$Enumeration_date, units = "years", precise = TRUE)
Cohort_Final <- Cohort[Cohort$Age < 100 & Cohort$Age >= 18,]
Cohort_Final$AgeCat <-cut(Cohort_Final$Age, c(0,40,100))
summary(Cohort_Final$AgeCat)
Cohort_Final$Age_Cat <- NA
Cohort_Final$Age_Cat[Cohort_Final$AgeCat=="(0,40]"] <- 0
Cohort_Final$Age_Cat[Cohort_Final$AgeCat=="(40,100]"] <- 1
table(Cohort_Final$Age_Cat)
Cohort_Final$Age_cat[Cohort_Final$AgeCat=="(0,40]"] <- "=< 40"
Cohort_Final$Age_cat[Cohort_Final$AgeCat=="(40,100]"] <- "> 40"
table(Cohort_Final$Age_cat)

Cohort_Final$Sex <- as.numeric(Cohort_Final$Sex)
Cohort_Final$Sex_B[Cohort_Final$Sex==1] <- 1
Cohort_Final$Sex_B[Cohort_Final$Sex==2] <- 0
Cohort_Final <- Cohort_Final[Cohort_Final$Sex_B == 0 | Cohort_Final$Sex_B== 1,]
Cohort_Final$Sex_B[is.na(Cohort_Final$Sex_B)] = 99
Cohort_Final <- subset(Cohort_Final, Sex_B < 2)

#lose 214 people due to age issues
#Lose 5 people due to unkown sexes



##


library(readr)
write_csv(Cohort,"/Volumes/Extreme SSD/Survival final/Cohort_test.csv")
Cohort <- read_csv("Cohort_test.csv")
Cohort$Age <- as.numeric(Cohort$Age)

write_csv(Cohort_Final,"/Volumes/Extreme SSD/Survival final/Cohort_Final_test.csv")
Cohort_Final <- read_csv("Cohort_Final_test.csv")
View(Cohort_Final)

NPR_Linkage <- read_excel("NPR_Link.xlsx")

Cohort_Final$Linked <- NPR_Linkage$LINKED[match(Cohort_Final$study_id, NPR_Linkage$PATIENT)]
Cohort_Final$NPR_Vital_Status <- NPR_Linkage$VITAL_STATUS[match(Cohort_Final$study_id, NPR_Linkage$PATIENT)]
Cohort_Final$NPR_Date <- NPR_Linkage$NPR_DEATH_D[match(Cohort_Final$study_id, NPR_Linkage$PATIENT)]

Cohort_Check <- select(Cohort_Final, study_id, Censored_Date, Censored, Linked, NPR_Vital_Status, NPR_Date)
Check_1 <- subset(Cohort_Check, Censored == 1 & Linked == 0)
Check_2 <- subset(Cohort_Check, Censored == 1 & Linked == 1)
Check_3 <- subset(Cohort_Check, Censored == 0 & Linked == 1)

Episodes <- read_excel("PHDC_Episode_Dates.xlsx")
Episodes$episode_start_date <- as.Date(Episodes$episode_start_date, "%Y-%m-%d")
table(Episodes$episode_short_description)

#Diabetes
Diabetes <- subset(Episodes, episode_short_description == "Diabetes Mellitus")
Diabetes <- Diabetes[order(Diabetes$study_id, Diabetes$episode_start_date),]
Diabetes$Enumeration_Date <- Cohort_Final$Enumeration_date[match(Diabetes$study_id, Cohort_Final$study_id)]
Diabetes$Enumeration_Date <- as.Date(Diabetes$Enumeration_Date, "%Y-%m-%d")
Diabetes$Time <- Diabetes$Enumeration_Date - Diabetes$episode_start_date
Current_Diabetes <- subset(Diabetes, Time >= -60)
Current_Diabetes$Current_Diabetes <- 1
Cohort_Final$Current_Diabetes <- Current_Diabetes$Current_Diabetes[match(Cohort_Final$study_id, Current_Diabetes$study_id)]
Cohort_Final$Current_Diabetes[is.na(Cohort_Final$Current_Diabetes)] = 0
table(Cohort_Final$Current_Diabetes)

#Hypertension
Hypertension <- subset(Episodes, episode_short_description == "Hypertension")
Hypertension <- Hypertension[order(Hypertension$study_id, Hypertension$episode_start_date),]
Hypertension$Enumeration_Date <- Cohort_Final$Enumeration_date[match(Hypertension$study_id, Cohort_Final$study_id)]
Hypertension$Enumeration_Date <- as.Date(Hypertension$Enumeration_Date, "%Y-%m-%d")
Hypertension$Time <- Hypertension$Enumeration_Date - Hypertension$episode_start_date
Current_Hypertension <- subset(Hypertension, Time >= -60)
Current_Hypertension$Current_Hypertension <- 1
Cohort_Final$Current_Hypertension <- Current_Hypertension$Current_Hypertension[match(Cohort_Final$study_id, Current_Hypertension$study_id)]
Cohort_Final$Current_Hypertension[is.na(Cohort_Final$Current_Hypertension)] = 0
table(Cohort_Final$Current_Hypertension)

#CKD
CKD <- subset(Episodes, episode_short_description == "Chronic Kidney Disease")
CKD <- CKD[order(CKD$study_id, CKD$episode_start_date),]
CKD$Enumeration_Date <- Cohort_Final$Enumeration_date[match(CKD$study_id, Cohort_Final$study_id)]
CKD$Enumeration_Date <- as.Date(CKD$Enumeration_Date, "%Y-%m-%d")
CKD$Time <- CKD$Enumeration_Date - CKD$episode_start_date
Current_CKD <- subset(CKD, Time >= -60)
Current_CKD$Current_CKD <- 1
Cohort_Final$Current_CKD <- Current_CKD$Current_CKD[match(Cohort_Final$study_id, Current_CKD$study_id)]
Cohort_Final$Current_CKD[is.na(Cohort_Final$Current_CKD)] = 0
table(Cohort_Final$Current_CKD)

#ARF
AKI <- subset(Episodes, episode_short_description == "Acute Kidney Injury")
AKI <- AKI[order(AKI$study_id, AKI$episode_start_date),]
AKI$Enumeration_Date <- Cohort_Final$Enumeration_date[match(AKI$study_id, Cohort_Final$study_id)]
AKI$Enumeration_Date <- as.Date(AKI$Enumeration_Date, "%Y-%m-%d")
AKI$Time <- AKI$Enumeration_Date - AKI$episode_start_date
Current_AKI <- subset(AKI, Time >= -60 & Time <= 60)
Current_AKI$Current_AKI <- 1
Cohort_Final$Current_AKI <- Current_AKI$Current_AKI[match(Cohort_Final$study_id, Current_AKI$study_id)]
Cohort_Final$Current_AKI[is.na(Cohort_Final$Current_AKI)] = 0
table(Cohort_Final$Current_AKI)

Previous_AKI <- subset(AKI, Time > 60)
Previous_AKI$Previous_AKI <- 1
Cohort_Final$Previous_AKI <- Previous_AKI$Previous_AKI[match(Cohort_Final$study_id, Previous_AKI$study_id)]
Cohort_Final$Previous_AKI[is.na(Cohort_Final$Previous_AKI)] = 0
table(Cohort_Final$Previous_AKI)

Incident_AKI <- subset(AKI, Time < -60)
Incident_AKI$Incident_AKI <- 1
Cohort_Final$Incident_AKI <- Incident_AKI$Incident_AKI[match(Cohort_Final$study_id, Incident_AKI$study_id)]
Cohort_Final$Incident_AKI[is.na(Cohort_Final$Incident_AKI)] = 0
table(Cohort_Final$Incident_AKI)

#TB (PHDC)
TB_PHDC <- subset(Episodes, episode_short_description == "TB positive")
TB_PHDC <- TB_PHDC[order(TB_PHDC$study_id, TB_PHDC$episode_start_date),]
TB_PHDC$Enumeration_Date <- Cohort_Final$Enumeration_date[match(TB_PHDC$study_id, Cohort_Final$study_id)]
TB_PHDC$Enumeration_Date <- as.Date(TB_PHDC$Enumeration_Date, "%Y-%m-%d")
TB_PHDC$Time <- TB_PHDC$Enumeration_Date - TB_PHDC$episode_start_date
Current_TB_PHDC <- subset(TB_PHDC, Time >= -60 & Time <= 60)
Current_TB_PHDC$Current_TB_PHDC <- 1
Cohort_Final$Current_TB_PHDC <- Current_TB_PHDC$Current_TB_PHDC[match(Cohort_Final$study_id, Current_TB_PHDC$study_id)]
Cohort_Final$Current_TB_PHDC[is.na(Cohort_Final$Current_TB_PHDC)] = 0
table(Cohort_Final$Current_TB_PHDC)

Previous_TB_PHDC <- subset(TB_PHDC, Time > 60)
Previous_TB_PHDC$Previous_TB_PHDC <- 1
Cohort_Final$Previous_TB_PHDC <- Previous_TB_PHDC$Previous_TB_PHDC[match(Cohort_Final$study_id, Previous_TB_PHDC$study_id)]
Cohort_Final$Previous_TB_PHDC[is.na(Cohort_Final$Previous_TB_PHDC)] = 0
table(Cohort_Final$Previous_TB_PHDC)

Incident_TB_PHDC <- subset(TB_PHDC, Time < -60)
Incident_TB_PHDC$Incident_TB_PHDC <- 1
Cohort_Final$Incident_TB_PHDC <- Incident_TB_PHDC$Incident_TB_PHDC[match(Cohort_Final$study_id, Incident_TB_PHDC$study_id)]
Cohort_Final$Incident_TB_PHDC[is.na(Cohort_Final$Incident_TB_PHDC)] = 0
table(Cohort_Final$Incident_TB_PHDC)

write_csv(Cohort_Final,"/Volumes/Extreme SSD/Survival final/Cohort_Final_test.csv")
Cohort_Final <- read_csv("Cohort_Final_test.csv")
View(Cohort_Final)

