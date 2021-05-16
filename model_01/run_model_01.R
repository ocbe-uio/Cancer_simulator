# Load R code
source("backend_v01.R")

# Initialize or load Vessel Data
VesselsList.26 <- initializeBloodVessels(plot = TRUE, BloVesN = 26, Rrows = 4, seed = 111)
saveRDS(object = VesselsList.26, file = "VesselsList.26.rds", ascii = FALSE, compress = 'xz')

# Alternatively, load from RDS
VesselsList <- readRDS(file = "VesselsList.26.rds")


#Same for all runs:
TimeInterval=1            # in hours
DT= 156                   # in hours
DTinR=108                 # in hours
G1=60                     # in hours
SimTime=1500              # in hours  (90 days)
PeriodForSenescence=2     # number of cell-cycles periods that lead to cell Senescence if the division did not happen during that period
PeriodForDeath=3          # number of cell-cycles periods that lead to cell death if the division did not happen during that period
IntravasationProb=0.01    # For a cell to be able to intravasate it should get a random number (between 0 and 1) that is equal to the "IntravasationProb" or less.
PER= 0.45                        # The persistence value (a value between 0 and 1)


#Simulation #7:
minSpe <- 3             # um/h
maxSpe <- 130              # um/h
desired_meanSpe <- 35     # um/h
desired_sdSpe <- 25       # um/h
LeavingRzone=0.5          # Threshold for leaving red zone    ( a value between 0 and 1)
BloVesN=VesselsList$vesselParams$BloVesN


y <- tryCatch({
    RunCellmigSimulation(
        TimeInterval = TimeInterval,
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
        SetCellN = 100000,
        IntravasationProb =  IntravasationProb)
}, error = function(e) {e})


