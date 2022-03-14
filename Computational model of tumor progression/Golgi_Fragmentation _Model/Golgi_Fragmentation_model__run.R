# Set wd
# You may need to change this PATH depending on the local installation folder
setwd("~/Documents/projects/GitHub_salim_models/export/v03/Model_02/")

# Load R code (backend)
source("model_02__backend.R")

# Initialize or load Vessel Data
# Uncomment to build a new VesselList and save an RDS file
#VesselsList <- initializeBloodVessels(plot = TRUE, BloVesN = 36, Rrows = 4, seed = 111)
#saveRDS(object = VesselsList, file = "VesselsList.rds", ascii = FALSE, compress = 'xz')

# or simply load a file you generated before
VesselsList <- readRDS(file = "VesselsList.rds")

# Initialize a subfolder to save run results
SUBF <- 'test_results_2'
dir.create(SUBF)
setwd(SUBF)


# Define params
TimeInterval=1            # in hours
DT= 156                   # Doubling Time in hours
DTinR=108                 # Doubling Time in Red zones in hours
G1=60                     # in hours
SimTime=1900               # Simulation time in hours

PeriodForDeath=3          # number of cell-cycles periods that lead to cell death if the division did not happen during that period
PeriodForSenescence=2     # number of cell-cycles periods that lead to cell Senescence if the division did not happen during that period
LeavingRzone=0.5          # Threshold for leaving red zone    ( a value between 0 and 1)

IntravasationProb=0.01    # For a cell to be able to intravasate it should get a random number (between 0 and 1) that is equal to the "IntravasationProb" or less.
IntraVas_inc_By_GF =33    # the increment in the IntravasationProb due to Golgi fragmentation
minRange_IntraVas_inc_By_GF=0.5
maxRange_IntraVas_inc_By_GF=1.5

# Golgi freq
Golgi_Fra_Frequency= 80

# Export params
exportSP=FALSE
exportRdata=FALSE


# Run Model 2
y <- RunCellmigSimulation_M2(
    TimeInterval = TimeInterval,
    DT = DT,
    DTinR = DTinR,
    G1 = G1,
    SimTime = SimTime,
    PeriodForSenescence = PeriodForSenescence,
    PeriodForDeath = PeriodForDeath,
    LeavingRzone = LeavingRzone,
    IntravasationProb = IntravasationProb,
    IntraVas_inc_By_GF = IntraVas_inc_By_GF,
    minRange_IntraVas_inc_By_GF = minRange_IntraVas_inc_By_GF,
    maxRange_IntraVas_inc_By_GF = maxRange_IntraVas_inc_By_GF,
    Golgi_Fra_Frequency = Golgi_Fra_Frequency,
    VesselsList = VesselsList,
    exportSP = exportSP,
    exportRdata = exportRdata)


# Expected results
#> Number of Metastatic Cells =  359
#> Number of Cells after  1900 hours =  6329
#> Total Number of Cells generated after  1900 hours =  12657

