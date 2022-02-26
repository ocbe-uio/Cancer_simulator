# Set wd
setwd("~/Documents/projects/cellmigRation/modeling/v01/")

# Load R code
source("modeling_backend_v001.R")


# Initialize or load Vessel Data
VesselsList <- initializeBloodVessels(BloVesN = 36, Rrows = 4, VesselsXY2CSVFile = "ves_xy.csv",
                                      LatDF2CSVFile = "lat_df.csv", seed = 1234)

# Alternatively
VesselsList <- loadBloodVessels(VesselsXY_CSVFile = "ves_xy.csv",
                                LatDF_CSVFile = "lat_df.csv", plot = TRUE)



# Define Simulation params
TimeInterval=1            # in hours
DT= 156                   # in hours
DTinR=108                 # in hours
G1=60                     # in hours
SimTime=1600              # in hours

minSpe <- 1.5             # um/h
maxSpe <- 64              # um/h
desired_meanSpe <- 15     # um/h
desired_sdSpe <- 11       # um/h

PER= 0.7                  # The persistence value (a value between 0 and 1)
PeriodForSenescence=2     # number of cell-cycles periods that lead to cell Senescence if the division did not happen during that period
PeriodForDeath=3          # number of cell-cycles periods that lead to cell death if the division did not happen during that period
LeavingRzone=0.5          # Threshold for leaving red zone    ( a value between 0 and 1)
IntravasationProb=0.01    # For a cell to be able to intravasate it should get a random number (between 0 and 1) that is equal to the "IntravasationProb" or less.


# Run Modeling Simulation
RunCellmigSimulation(TimeInterval = TimeInterval,
                     DT = DT,
                     DTinR = DTinR,
                     G1 = G1,
                     SimTime = SimTime,
                     minSpe = minSpe,
                     maxSpe = maxSpe,
                     desired_meanSpe = desired_meanSpe,
                     desired_sdSpe = desired_sdSpe,
                     VesselsList = VesselsList,
                     PER = PER,
                     PeriodForSenescence = PeriodForSenescence,
                     PeriodForDeath = PeriodForDeath,
                     LeavingRzone = LeavingRzone,
                     IntravasationProb =  IntravasationProb)


















