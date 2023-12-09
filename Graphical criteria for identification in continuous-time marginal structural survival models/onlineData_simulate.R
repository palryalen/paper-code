
# Install and load the devtools and miscFunctions packages
library(devtools)
# install_github("palryalen/miscFunctions")
library(miscFunctions)
library(data.table)

# Transition hazards are defined as smooth.spline objects. These are automatically loaded when running this script. Users who want to simulate realizations according to their own choice of transition hazards can do so by building on the following simple example:
# time.points = 1:10
# haz.values = c( rep(0.1,5), rep(0.3,5) )
# haz_DNA = smooth.spline(time.points, haz.values)

# Load simulation_transition_hazards.RData file
# save(simulation_transition_hazards,file="simulation_transition_hazards.RData")
load("simulation_transition_hazards.RData")

# Extract relevant hazard functions
haz_DNA = simulation_transition_hazards$haz_DNA
haz_RNA = simulation_transition_hazards$haz_RNA
haz_DNA_ci2p = simulation_transition_hazards$haz_DNA_ci2p
haz_RNA_ci2p = simulation_transition_hazards$haz_RNA_ci2p



# Set the end time for simulation
endTime <- 4

# Create a transition matrix with different possible transition states
transitionMatrix <- vector("list",length = 5)

# Define transition probabilities and neighbors for each state
transitionMatrix[[1]]$smoothspline[[1]] <- list(NULL)
transitionMatrix[[1]]$smoothspline[[2]] <- list(haz_DNA)
transitionMatrix[[1]]$smoothspline[[3]] <- list(NULL)
transitionMatrix[[1]]$smoothspline[[4]] <- list(NULL)
transitionMatrix[[1]]$smoothspline[[5]] <- list(NULL)
transitionMatrix[[1]]$neighbours <- c(2)
transitionMatrix[[2]]$smoothspline[[1]] <- list(NULL)
transitionMatrix[[2]]$smoothspline[[2]] <- list(NULL)
transitionMatrix[[2]]$smoothspline[[3]] <- list(NULL)
transitionMatrix[[2]]$smoothspline[[4]] <- list(NULL)
transitionMatrix[[2]]$smoothspline[[5]] <- list(haz_DNA_ci2p)
transitionMatrix[[2]]$neighbours <- c(5)
transitionMatrix[[3]]$smoothspline[[1]] <- list(NULL)
transitionMatrix[[3]]$smoothspline[[2]] <- list(NULL)
transitionMatrix[[3]]$smoothspline[[3]] <- list(NULL)
transitionMatrix[[3]]$smoothspline[[4]] <- list(haz_RNA)
transitionMatrix[[3]]$smoothspline[[5]] <- list(NULL)
transitionMatrix[[3]]$neighbours <- c(4)
transitionMatrix[[4]]$smoothspline[[1]] <- list(NULL)
transitionMatrix[[4]]$smoothspline[[2]] <- list(NULL)
transitionMatrix[[4]]$smoothspline[[3]] <- list(NULL)
transitionMatrix[[4]]$smoothspline[[4]] <- list(NULL)
transitionMatrix[[4]]$smoothspline[[5]] <- list(haz_RNA_ci2p)
transitionMatrix[[4]]$neighbours <- c(5)
transitionMatrix[[5]]$smoothspline[[1]] <- list(NULL)
transitionMatrix[[5]]$smoothspline[[2]] <- list(NULL)
transitionMatrix[[5]]$smoothspline[[3]] <- list(NULL)
transitionMatrix[[5]]$smoothspline[[4]] <- list(NULL)
transitionMatrix[[5]]$smoothspline[[5]] <- list(NULL)
transitionMatrix[[5]]$neighbours <- NULL

# Number of individuals in each group
n_DNA <- 858 
n_RNA <- 878

# Number of individuals under follow-up at time 0
n_DNA__ <- 4
n_RNA__ <- 11

# Simulate data for each group and time under follow-up
dfr1 <- data.table(generateFrame(n_DNA-n_DNA__, endTime, transitionMatrix, 1))
dfr1__ <- data.table(generateFrame(n_DNA__, endTime, transitionMatrix, 2))
dfr2 <- data.table(generateFrame(n_RNA-n_RNA__, endTime, transitionMatrix, 3))
dfr2__ <- data.table(generateFrame(n_RNA__, endTime, transitionMatrix, 3))

# Adjust IDs to avoid conflicts between groups
dfr1__[, id := id + max(dfr1$id)]
dfr2[, id := id + max(dfr1__$id)]
dfr2__[, id := id + max(dfr2$id)]

# Combine datasets
dfr1 <- rbind(dfr1, dfr1__)
dfr2 <- rbind(dfr2, dfr2__)
dfr <- rbind(dfr1, dfr2)

# Remove terminating state 5
dfr <- dfr[!(dfr$from.state == 5), ]

# Label states and technologies
dfr[, f.state := ""]
dfr[, t.state := ""]
dfr[, technology := ""]
dfr[from.state %in% c(1, 2), technology := "DNA"]
dfr[from.state %in% c(3, 4), technology := "RNA"]
dfr[from.state %in% c(2, 4), f.state := "follow-up"]
dfr[to.state %in% c(2, 4), t.state := "follow-up"]
dfr[to.state %in% c(5), t.state := "cin2+"]
dfr[, from.state := f.state]
dfr[, to.state := t.state]
dfr[from.state == to.state, to.state := ""]
dfr <- subset(dfr, select = !(names(dfr) %in% c("f.state", "t.state")))


# sim_data is the simulated data set. The analysis can be run in example_analysis.R
sim_data = dfr

# Save simulated data:
# save(sim_data, file = "sim_data.RData")


