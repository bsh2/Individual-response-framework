# ==============================================================================
# CONTACT
# ==============================================================================

# b.stephens-hemingway[at]rgu.ac.uk

# ==============================================================================
# IMPORT DATA
# ==============================================================================

mock_data <- read.csv('mock_data.csv')

# ==============================================================================
# CALCULATE TYPICAL ERROR FROM TEST-RETEST RELIABILITY TRIALS
# ==============================================================================
# Step 1: Calculate test-retest difference scores [Test(2)-Test(1)] for each
#         individual across both outcome measures (S7SF, MCARN)
# Step 2: Calculate standard deviation of the difference scores
# Step 3: Estimate TE by dividing the standard deviation of the difference
#         scores by the positive square root of 2.

calculate_TE <- function(RT1,RT2){
  diff <- RT2 - RT1
  stdev <- sd(diff)
  te <- stdev/sqrt(2)
  return(te)
}

S7SF_TE <- calculate_TE(mock_data$S7SF_RT1, mock_data$S7SF_RT2)
MCARN_TE <- calculate_TE(mock_data$MCARN_RT1, mock_data$MCARN_RT2)

# ==============================================================================
# COMPUTE TYPICAL ERROR FROM PUBLISHED RELIABILITY DATA (CoV)
# ==============================================================================
# See SF-S6 of the supplementary spreadsheet and paper for detailed explanation
# Take CV% as 4.94% for total work done (kJ) in the CCT110 exhaustion test
# TE Estimate = [CV * Mean(Baseline)]/100

CCT110_TE <- (mean(mock_data$CCT110_PRE) * 4.94) / 100

# ==============================================================================
# COMPUTE n-ADJUSTED TRUE-SCORE CONFIDENCE INTERVALS OF ANY WIDTH FOR BASELINE
# VALUES, USING THE T-DISTRIBUTION
# ==============================================================================
# Requires ---------------------------------------------------------------------
  ci_width_ts <- 90          # Adjusted True Score CI width in % (1-99)
# ------------------------------------------------------------------------------

TS_CI_Adj <- function(base_vals, TE){
    # Inputs *******************************************************************
    #   1. base_vals: Vector of individual pre-intervention testing values for
    #      measure being used. DF = Length(base_vals)-1
    #   2. TE: Typical error value already calculated for the measure
    # **************************************************************************
    tdf <- data.frame(participant = c(1:length(base_vals)), base_vals)
    t_c <- abs(qt((((100-ci_width_ts)/100)/2), length(tdf$base_vals)-1))
    tdf$LB <- tdf$base_vals - (t_c * TE)
    tdf$UB <- tdf$base_vals + (t_c * TE)
    return(tdf[,c(1,3,4)])
    # Returns dataframe containing upper and lower bound values
  }

S7SF_TSCI <- TS_CI_Adj(mock_data$S7SF_PRE, S7SF_TE)
MCARN_TSCI <- TS_CI_Adj(mock_data$MCARN_PRE, MCARN_TE)
CCT110_TSCI <- TS_CI_Adj(mock_data$CCT110_PRE, CCT110_TE)

# ==============================================================================
# COMPUTE n-ADJUSTED CHANGE SCORE CIs (OF ANY WIDTH) USING PRE-POST DIFFERENCE
# SCORES FOLLOWING THE INTERVENTION, USING THE T-DISTRIBUTION.
# ==============================================================================

# Requires ---------------------------------------------------------------------
  ci_width_cs <- 90          # Adjusted Change Score CI width in % (1-99)
  # Typical error value calculated in one of previous sections
  # Pre and post-intervention measured values
# ------------------------------------------------------------------------------

CS_CI_Adj <- function(pre_vals, post_vals, TE){
    # Inputs *******************************************************************
    #   1. pre/post_vals : Pre-and-post intervention measured values
    #   2. TE: Typical error value already calculated for the measure
    # **************************************************************************
    tdf <- data.frame(participant = c(1:length(pre_vals)), difference = 
              (post_vals - pre_vals))
    t_c <- abs(qt((((100-ci_width_cs)/100)/2), length(tdf$participant)-1))
    tdf$LB <- tdf$difference - (sqrt(2) * t_c * TE)
    tdf$UB <- tdf$difference + (sqrt(2) * t_c * TE)
    return(tdf)
    # Returns dataframe containing upper and lower bound values of the interval
    # of the change score.
  }

S7SF_CSCI <- CS_CI_Adj(mock_data$S7SF_PRE, mock_data$S7SF_POST, S7SF_TE)
MCARN_CSCI <- CS_CI_Adj(mock_data$MCARN_PRE, mock_data$MCARN_POST, MCARN_TE)
CCT110_CSCI <- CS_CI_Adj(mock_data$CCT110_PRE, mock_data$CCT110_POST, CCT110_TE)

# ==============================================================================
# ZERO-BASED, AND SWC-BASED THRESHOLDS FOR ASSESSING MEANINGFUL CHANGE
# ==============================================================================

# Requires ---------------------------------------------------------------------
  # Adjusted change score CI (upper and lower bound values) calculated in the 
  # previous section [for desired CI width]
  SWC_factor <- 0.2   
    # You could pick < 0.2 (Trivial), 0.2-0.5 (small), 0.5-0.8 (moderate) or 
    # >0.8 (large). See Hopkins (2004) for more detail.
  SWC_MCARN <- sd(mock_data$MCARN_PRE) * SWC_factor
  SWC_S7SF <- -sd(mock_data$S7SF_PRE) * SWC_factor
    # Note the - because we are looking for skinfold reduction post-intervention
    # So threshold will be LB < SWC and UB < SWC for S7SF. We use ifelse() to 
    # control for this.
  SWC_CCT110 <- sd(mock_data$CCT110_PRE) * SWC_factor
# ------------------------------------------------------------------------------

CSCI_Threshold <- function(LB, UB, SWC, negative){
  # negative: TRUE means that improvement in the negative direction is
  # considered good for the measure considered. FALSE implies improvement is
  # considered as positive change in the measure.
  ifelse(negative == TRUE, zero_satisfy <- ifelse(LB < 0 & UB < 0, TRUE, FALSE),
    zero_satisfy <- ifelse(LB > 0 & UB > 0, TRUE, FALSE))
  ifelse(SWC > 0, SWC_satisfy <- ifelse(LB > SWC & UB > SWC, TRUE, FALSE), 
    SWC_satisfy <- ifelse(LB < SWC & UB < SWC, TRUE, FALSE))
  tempdf <- data.frame(participant = c(1:length(LB)), LB, UB, zero_satisfy,
    SWC_satisfy)
  return(tempdf)
  # Function simply returns a true/false in the relevant 'satisfy' column, 
  # indicating whether both bounds of the CI lie to the correct side of the 
  # threshold considered.
}

S7SF_Threshold <- CSCI_Threshold(S7SF_CSCI$LB, S7SF_CSCI$UB, SWC_S7SF, T)
MCARN_Threshold <- CSCI_Threshold(MCARN_CSCI$LB, MCARN_CSCI$UB, SWC_MCARN, F)
CCT110_Threshold <- CSCI_Threshold(CCT110_CSCI$LB, CCT110_CSCI$UB, SWC_CCT110,F)

# ==============================================================================
# ESTIMATING PROPORTION OF RESPONSE
# ==============================================================================

# To estimate the proportion of 'response' and 'non-response' to a chronic 
# intervention, we first estimate the variability in change scores directly
# attributable to the intervention, and the variability attributed to 
# biological variation and measurement error. This partitioning of variance
# is achieved by comparing standard deviation of change scores between the
# intervention group and control group.

mock_data$S7SF_DIFF <- mock_data$S7SF_POST - mock_data$S7SF_PRE
mock_data$MCARN_DIFF <- mock_data$MCARN_POST - mock_data$MCARN_PRE
mock_data$CCT110_DIFF <- mock_data$CCT110_POST - mock_data$CCT110_PRE

# A bit pushed for time here, hence the repetitive code, 
# the better implementation would be to write a function to do this, 
# rather than my repetitive code! When I have time I will come back and rewrite.

response <- data.frame(
  Measure = c("S7SF", "MCARN", "CCT110"),
  Int_SD = c(sd(unlist(subset(mock_data, GROUP == 1, S7SF_DIFF))),
             sd(unlist(subset(mock_data, GROUP == 1, MCARN_DIFF))),
             sd(unlist(subset(mock_data, GROUP == 1, CCT110_DIFF)))),
  Int_Mean = c(mean(unlist(subset(mock_data, GROUP == 1, S7SF_DIFF))),
               mean(unlist(subset(mock_data, GROUP == 1, MCARN_DIFF))),
               mean(unlist(subset(mock_data, GROUP == 1, CCT110_DIFF)))),
  Ctr_SD = c(sd(unlist(subset(mock_data, GROUP == 2, S7SF_DIFF))),
             sd(unlist(subset(mock_data, GROUP == 2, MCARN_DIFF))),
             sd(unlist(subset(mock_data, GROUP == 2, CCT110_DIFF)))),
  SWC = c(SWC_S7SF, SWC_MCARN, SWC_CCT110)
)

# Calculating SD_IR using SD_IR = sqrt((SD_INT)^2-(SD_CTR)^2)
response$SD_IR <- sqrt((descriptive_stats$Int_SD)^2 - 
  (descriptive_stats$Ctr_SD)^2)

# The distribution of interest is modelled as a normal distribution, centered at
# the mean observed change score with standard deviation equal to SD_IR.

response$Proportion <- ifelse(response$SWC < 0, ((pnorm(response$SWC, 
  response$Int_Mean,response$SD_IR)) * 100) , 
    ((1-pnorm(response$SWC,response$Int_Mean, response$SD_IR)) * 100))

# ==============================================================================
# BOOTSTRAPPING
# ==============================================================================

# As an additional series of steps, you can obtain confidence limits for the 
# estimated proportion of response using bootstrapping, which is a resampling
# technique that uses the observed sample to model the population distribution.

library(boot)
boot_ci <- 90
nbootstraps <- 1000

propboot <- function(ctrl_data, int_data, swc, ci_width, bootstraps){
  # Inputs:
  #   ctrl_data: Vector of control group difference scores
  #   int_data: Vector of intervention group difference scores
  #   swc: single numerical value for SWC
  #   ci_width: single numerical value (1-99) [units %] for CI width
  #   boostraps: number of bootstraps to compute
  counter <- 0
  proportion_response <- as.numeric()
  for (i in 1:bootstraps){
    tv_ctr <- sample(ctrl_data, length(ctrl_data), replace = TRUE)
    tv_int <- sample(int_data, length(int_data), replace = TRUE)
    mean_ctr <- mean(tv_ctr)
    mean_int <- mean(tv_int)
    sd_ctr <- sd(tv_ctr)
    sd_int <- sd(tv_int)
    if (sd_int > sd_ctr){
      sd_ir <- sqrt( (sd_int)^2 - (sd_ctr)^2 )
        if (swc < 0){
          prop_resp <- as.numeric((pnorm(swc, mean_int, sd_ir) * 100))
        } 
        if (swc >= 0){
          prop_resp <- as.numeric(((1 - pnorm(swc, mean_int, sd_ir)) * 100))
        }
      proportion_response[length(proportion_response)+1] <- prop_resp
    } else {counter <- counter + 1}
  }
  print("Number of bootstraps completed:")
  print(bootstraps - counter)
  output <- quantile(proportion_response, 
   probs = c(((1-(ci_width/100))/2),((ci_width/100) + ((1 - (ci_width/100))/2))))
  return(output)
}

ctr_S7SF <- as.numeric(unlist(subset(mock_data, GROUP == 2, S7SF_DIFF)))
int_S7SF <- as.numeric(unlist(subset(mock_data, GROUP == 1, S7SF_DIFF)))
ctr_MCARN <- as.numeric(unlist(subset(mock_data, GROUP == 2, MCARN_DIFF)))
int_MCARN <- as.numeric(unlist(subset(mock_data, GROUP == 1, MCARN_DIFF)))
ctr_CCT110 <- as.numeric(unlist(subset(mock_data, GROUP == 2, CCT110_DIFF)))
int_CCT110 <- as.numeric(unlist(subset(mock_data, GROUP == 1, CCT110_DIFF)))

# Run the bootstraps

propboot(ctr_S7SF, int_S7SF, SWC_S7SF, boot_ci, nbootstraps)
propboot(ctr_MCARN, int_MCARN, SWC_MCARN, boot_ci, nbootstraps)
propboot(ctr_CCT110, int_CCT110, SWC_CCT110, boot_ci, nbootstraps)