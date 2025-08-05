Co_tri_Data <- read.csv("/Volumes/Extreme SSD/Survival final/Co_tri_Data.csv")
flu_data <- read.csv("/Volumes/Extreme SSD/Survival final/Cohort_fluconazole_data.csv")
crypto_data <- read_excel("/Volumes/Extreme SSD/Survival final/Cryptococcus.xlsx")
cohort <- read.csv("/Volumes/Extreme SSD/Survival final/Cohort_Final_1.csv")


#####
#reflex LFA positive data:
colnames(flu_data)
colnames(cohort)
table(crypto_data$common_method)

# Step 1: Filter crypto_data for LFA Reflex Screening
lfa_reflex_data <- crypto_data[crypto_data$common_method == "LFA Reflex Screening", ]

# Step 2: Keep only people who are in the study cohort
lfa_reflex_in_cohort <- lfa_reflex_data[lfa_reflex_data$study_id %in% cohort$study_id, ]

# View the number of records and unique patients
nrow(lfa_reflex_in_cohort)          # number of test records
length(unique(lfa_reflex_in_cohort$study_id))  # number of unique people

##only include test around enumeration date:
# Step 1: Make sure dates are proper Date class
lfa_reflex_in_cohort$lab_assumed_taken_date <- as.Date(lfa_reflex_in_cohort$lab_assumed_taken_date)
cohort$Enumeration_date <- as.Date(cohort$Enumeration_date)

# Step 2: Merge Enumeration_date onto the LFA dataset
library(dplyr)

lfa_reflex_in_cohort <- lfa_reflex_in_cohort %>%
  left_join(cohort %>% select(study_id, Enumeration_date, Enumeration_CD4), by = "study_id")

# Step 3: Calculate difference between test date and enumeration date
lfa_reflex_in_cohort <- lfa_reflex_in_cohort %>%
  mutate(diff_days = as.numeric(lab_assumed_taken_date - Enumeration_date))

# Step 4: Filter for tests within -7 and +14 days
lfa_reflex_restricted <- lfa_reflex_in_cohort %>%
  filter(diff_days >= 0 & diff_days <= 14)

# Optional: how many patients?
length(unique(lfa_reflex_restricted$study_id))

# seemingly 100 people did not have reflex LFAs done?
# Step 1: Start from the cohort with CD4 < 100
cohort_under_100_cd4 <- subset(cohort, Enumeration_CD4 < 100)

# Step 2: Find study_ids of people who had LFA in window
study_ids_with_lfa <- unique(lfa_reflex_restricted$study_id)
length(study_ids_with_lfa)
# Step 3: Find people in cohort_under_100_cd4 who did NOT have an LFA
cohort_under_100_no_lfa <- cohort_under_100_cd4[!(cohort_under_100_cd4$study_id %in% study_ids_with_lfa), ]

# Step 4: Check how many
nrow(cohort_under_100_no_lfa)


#112 of 5903 people with enumeration CD4 counts <100 didn't have a reflex LFA done

## How many were positive and subsequently had a CSF crytpo invesitgation within 30 days?
table(crypto_data$lab_specimen_type)
table(crypto_data$common_result)

## Positives
# Keep only positive Reflex LFA results
lfa_reflex_restricted_positive <- lfa_reflex_restricted %>%
  filter(common_result == "Positive")

# Now, keep only positives among people with Enumeration CD4 < 100
lfa_reflex_positive_under_100 <- lfa_reflex_restricted_positive %>%
  filter(study_id %in% cohort_under_100_cd4$study_id)

# How many positives?
nrow(lfa_reflex_positive_under_100)
length(unique(lfa_reflex_positive_under_100$study_id))  # unique patients

## CSF
# First, define the CSF-related specimen types
csf_specimen_types <- c("CSF", "FLUID CSF", "BLOOD AND CSF", "BCSF", "BLC,C,CF", "C,CF")

# Subset for proper CSF specimens OR "UNKNOWN" + India Ink
csf_tests <- crypto_data %>%
  filter(
    lab_specimen_type %in% csf_specimen_types |
      (lab_specimen_type == "UNKNOWN" & common_method == "India Ink")
  )

# Confirm
table(csf_tests$lab_specimen_type)

# Make sure the test dates are correct
csf_tests$lab_assumed_taken_date <- as.Date(csf_tests$lab_assumed_taken_date)

# First make a lookup table of positive LFA study_ids and dates
lfa_positive_dates <- lfa_reflex_positive_under_100 %>%
  select(study_id, lab_assumed_taken_date) %>%
  rename(reflex_date = lab_assumed_taken_date)

# Join CSF tests to positive LFAs by study_id
csf_after_lfa <- csf_tests %>%
  inner_join(lfa_positive_dates, by = "study_id") %>%
  mutate(diff_days = as.numeric(lab_assumed_taken_date - reflex_date)) %>%
  filter(diff_days >= 0 & diff_days <= 30)
# How many patients had a CSF test within 30 days?
length(unique(csf_after_lfa$study_id))

## Fluconazole
colnames(flu_data)
# Make sure the pharmacy issue dates are Date class
flu_data$pharm_issue_date <- as.Date(flu_data$pharm_issue_date)
# Join fluconazole prescriptions to positive Reflex LFA patients
fluconazole_after_lfa <- flu_data %>%
  inner_join(lfa_positive_dates, by = "study_id") %>%
  mutate(diff_days = as.numeric(pharm_issue_date - reflex_date)) %>%
  filter(diff_days >= 0 & diff_days <= 90)

length(unique(fluconazole_after_lfa$study_id))

#135 people recieved fluconazole within 90 days after a positive reflex LFA.

##### Bactrim
colnames(Co_tri_Data)
# Convert pharm_issue_date to Date if needed
Co_tri_Data$pharm_issue_date <- as.Date(Co_tri_Data$pharm_issue_date)
# Create a lookup table of study_id and Enumeration_date
enumeration_dates <- cohort %>%
  select(study_id, Enumeration_date)
library(dplyr)

# Join prescriptions to Enumeration Dates
co_trimox_with_enum <- Co_tri_Data %>%
  inner_join(enumeration_dates, by = "study_id") %>%
  mutate(diff_days = as.numeric(pharm_issue_date - Enumeration_date))

# Keep only prescriptions from 0 to 90 days after Enumeration
co_trimox_within_90d <- co_trimox_with_enum %>%
  filter(diff_days >= 0 & diff_days <= 90)
# How many unique patients
length(unique(co_trimox_within_90d$study_id))

#5857 recieved Bactrim within 90 days of being enrolled into the cohort.

## Any bactrim:
# How many cohort patients ever received Bactrim (at any time)?
co_trimox_ever <- Co_tri_Data %>%
  filter(study_id %in% cohort$study_id)

length(unique(co_trimox_ever$study_id))

#9277 at any point during the study period

#### Different KM curves:
cohort <- read.csv("/Volumes/Extreme SSD/Survival final/Cohort_Final_1.csv")
summary(cohort)
table(cohort$ART_Exp_C)
Cohort <- read.csv("/Volumes/Extreme SSD/Survival final/Cohort_Final_test.csv")
summary(Cohort)
table(Cohort$ART_Exp_C)
library(survival)
library(survminer)
library(dplyr)

# Prepare cohort 1
cohort <- cohort %>%
  mutate(
    Time = ifelse(Time > 1500, 1500, Time),
    Censored = ifelse(Time == 1500, 0, Censored),
    ART_Exp_C = as.factor(ART_Exp_C)
  )

# Prepare cohort 2
Cohort <- Cohort %>%
  mutate(
    Time = ifelse(Time > 1500, 1500, Time),
    Censored = ifelse(Time == 1500, 0, Censored),
    ART_Exp_C = as.factor(ART_Exp_C)
  )

# Fit survival models separately
fit1 <- survfit(Surv(Time, Censored) ~ ART_Exp_C, data = cohort)
fit2 <- survfit(Surv(Time, Censored) ~ ART_Exp_C, data = Cohort)

# Plot for Original cohort
plot1 <- ggsurvplot(
  fit1,
  data = cohort,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  xlim = c(0, 1500),
  ylim = c(0.75, 1.0),
  break.time.by = 250,
  ggtheme = theme_bw(),
  title = "Original Cohort: Kaplan-Meier by ART Experience",
  legend.title = "ART Experience",
  risk.table.y.text = FALSE
)

# Plot for Reclassified cohort
plot2 <- ggsurvplot(
  fit2,
  data = Cohort,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  xlim = c(0, 1500),
  ylim = c(0.75, 1.0),
  break.time.by = 250,
  ggtheme = theme_bw(),
  title = "Reclassified Cohort (90d naive): Kaplan-Meier by ART Experience",
  legend.title = "ART Experience",
  risk.table.y.text = FALSE
)

# Display both plots
print(plot1)
print(plot2)
