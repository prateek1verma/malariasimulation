####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
outFolder <- "mgdrive"

####################
# Simulation Parameters
####################
# days to run the simulation
tMax <- 365

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 1-node network where mosquitoes do not leave
moveMat <- matrix(data = 1, nrow = 1, ncol = 1)

# parameters of the population equilibrium
adultPopEquilibrium <- 4000
sitesNumber <- nrow(moveMat)

####################
# Basic Inheritance pattern
####################
# Mendelian cube with standard (default) parameters
cube <- cubeMendelian()

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
#  Releases start at time 25, occur every day, for 1 day.
#  There are 10 mosquitoes released every time.
releasesParameters <- list(releasesStart=50,
                           releasesNumber=10,
                           releasesInterval=7,
                           releaseProportion=200)

# generate release vector
releasesVector <- generateReleaseVector(driveCube=cube,
                                        releasesParameters=releasesParameters)

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- releasesVector


# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                      sexProbs=c(.5,.5),
                                      numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# set MGDrivE to run deterministic
setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we havee a network of 1 population.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, sampTime = 1, nPatch=sitesNumber,
                              beta=bioParameters$betaK, muAd=bioParameters$muAd,
                              popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                              tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                              AdPopEQ=adultPopEquilibrium, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                          driveCube=cube,
                          patchReleases=patchReleases,
                          migrationMale=moveMat,
                          migrationFemale=moveMat,
                          migrationBatch=batchMigration,
                          directory=outFolder,
                          verbose=FALSE)

# run simulation
MGDrivESim$oneRun(verbose = FALSE)

####################
# Post Analysis
####################
# split output by patch
#  Required for plotting later
splitOutput(readDir = outFolder, remFile = TRUE, verbose = FALSE)

# aggregate females by their mate choice
#  This reduces the female file to have the same columns as the male file
aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID,
                 remFile = TRUE, verbose = FALSE)

# plot output to see effect
plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)

