# script to simulate overdominance without migration...
# Richmond change...
setwd("~/R")
# file name info...
scenarioname <- "scenario1"
dir.create(scenarioname)
# what is the total pop size
population_size <- 100
# how many offspring does each mating produce?
number_offspring <- 10
# how many generations are we simulating
number_generations <- 500
# how many loci are being modeled
number_loci <- 2
# does recombination happen in males? TRUE = yes, FALSE = no
male_recombination <- TRUE
# how many replicates to run
replicates <- 100
# start clock
st<-proc.time()[1]
# need to loop over a range of recombination rates...
for (recombination_rate in seq(0.1,0.1,0.1)) {
  message("recombination_rate = ", recombination_rate) # to track simulation progress
  # need to loop over a range of mortality rates due to homozygosity
  for (m in seq(0.15, 0.15)) {
    message("m = ", m) # to track simulation progress
    # need to loop through a series of replicate runs...
    # building starting population
    for (i in 1:replicates) {
      message("i = ", i) # to track simulation progress
      alleles_names <- list()
      alleles_names[[1]] <- c("A1", "A2") # Name the alleles from locus under selection 
      alleles_names[[2]] <- c("B1", "B2") # Name the alleles from neutral locus
      alleles_frequencies <- list()
      alleles_frequencies[[1]] <- c(0.2, 0.8) # Allele frequencies for locus A
      alleles_frequencies[[2]] <- c(0.5, 0.5) # Allele frequencies for locus B
      # Create a data.frame to track details...
      starting_population <- as.data.frame(matrix(ncol = number_loci * 2 + 6, nrow = population_size))
      starting_population[1:(population_size / 2), 1] <- "Male" # first column stores sex
      starting_population[(population_size / 2 + 1):population_size, 1] <- "Female" # first column stores sex
      starting_population[, 2] <- 1 # second column stores population
      starting_population[, 3] <- 0 # third column stores generation
      starting_population[, 4] <- 0 # fourth column stores male parent
      starting_population[, 5] <- 0 # fifth column stores female parent
      starting_population[, 6] <- rownames(starting_population) # sixth column stores generation individual number
      # generate initial genotypes for population 1
      # going to put the two loci on a single stretch of chromosome one after the other in the data table
      for (k in 1:number_loci) {
        starting_population[1:population_size, (k + 6)] <- sample(alleles_names[[k]], size = population_size, replace = TRUE, prob = alleles_frequencies[[k]])
        starting_population[1:population_size, (k + 8)] <- sample(alleles_names[[k]], size = population_size, replace = TRUE, prob = alleles_frequencies[[k]])
      }
      # We have our population genomes, now we have to simulate the population over time...
      pop1 <- starting_population
      
      # heterozygosity variables
      het_tot_A = list()
      het_tot_B = list()
      generation_tot = list()
      alle_freq_tot_A1 = list()
      
      for (g in 1:number_generations) {
        message("g = ", g)
        # males and females from pop 1
        male_parents <- which(starting_population$V1 == "Male")
        female_parents <- which(starting_population$V1 == "Female")
        
        # create matrix into which parent combinations are placed
        parents_matrix <- matrix(nrow = population_size, ncol = 2)
        
        # select males and females that are mating with each other - males column 1 females column 2
        parents_matrix[, 1] <- sample(as.numeric(male_parents), size = population_size, replace = TRUE, prob = NULL)
        parents_matrix[, 2] <- sample(as.numeric(female_parents), size = population_size, replace = TRUE, prob = NULL)
        
        # have parents, now lets make number_offspring making number_offspring*number_offspring/pairing young
        offspring_matrix <- as.data.frame(matrix(nrow = population_size * number_offspring, ncol = (7 + number_loci * 2)))
        
        # now have to loop through the number_offspring to make offspring...
        for (parent in 1:dim(parents_matrix)[1]) {
          offspring_matrix[((parent - 1) * number_offspring + 1):(parent * number_offspring), 1] <-
            sample(c("Male", "Female"), size = number_offspring, replace = TRUE) # sex
          offspring_matrix[((parent - 1) * number_offspring + 1):(parent * number_offspring), 2] <- 1 # source population
          offspring_matrix[((parent - 1) * number_offspring + 1):(parent * number_offspring), 3] <- g # generation number
          offspring_matrix[((parent - 1) * number_offspring + 1):(parent * number_offspring), 4] <- parents_matrix[parent, 1] # male parent
          offspring_matrix[((parent - 1) * number_offspring + 1):(parent * number_offspring), 5] <- parents_matrix[parent, 2] # female parent
          
          # 6 is reserved for ind # and 7 is reserved for alive/dead status
          
          # getting genotype of each number_offspring for a pairing...
          
          # parent chromosomes
          mchrome1 <- c(starting_population[parents_matrix[parent, 1], 7], starting_population[parents_matrix[parent, 1], 8])
          mchrome2 <- c(starting_population[parents_matrix[parent, 1], 9], starting_population[parents_matrix[parent, 1], 10])
          fchrome1 <- c(starting_population[parents_matrix[parent, 2], 7], starting_population[parents_matrix[parent, 2], 8])
          fchrome2 <- c(starting_population[parents_matrix[parent, 2], 9], starting_population[parents_matrix[parent, 2], 10])
          
          for (n in 1:number_offspring) {
            # getting the genotype from dad
            randad <- runif(1) # getting a random deviate to see if crossing over happens
            mctemp1 <- rep(NA, 2)
            mctemp2 <- rep(NA, 2)
            if (randad < recombination_rate & male_recombination) {
              mctemp1 <- c(mchrome1[1], mchrome2[2])
              mctemp2 <- c(mchrome2[1], mchrome1[2])
            } else {
              mctemp1 <- mchrome1 # store 1st male gamete from dad
              mctemp2 <- mchrome2 # store 2nd male gamete from dad
            }
            malechromes <- list(mctemp1, mctemp2)
            
            # pick which male gamete this number_offspring is getting...
            offspring_matrix[(parent - 1) * number_offspring + n, 8:9] <- malechromes[[sample(c(1, 2), 1)]]
            
            # getting the genotype from mom
            randmom <-runif(1) # getting a random deviate to see if crossing over happens
            fctemp1 <- rep(NA, 2)
            fctemp2 <- rep(NA, 2)
            if (randmom < recombination_rate) {
              fctemp1 <- c(fchrome1[1], fchrome2[2])
              fctemp2 <- c(fchrome2[1], fchrome1[2])
            } else {
              fctemp1 <- fchrome1 # store 1st female gamete from mom
              fctemp2 <- fchrome2 # store 2nd female gamete from mom
            }
            femalechromes <- list(fctemp1, fctemp2)
            # pick which male gamete this number_offspring is getting...
            offspring_matrix[(parent - 1) * number_offspring + n, 10:11] <-
              femalechromes[[sample(c(1, 2), 1)]]
          }
        }
        offspring_matrix[, 6] <- as.numeric(rownames(offspring_matrix))
        offspring_matrix[, 7] <- 1 # "1 = alive"
        offspring_matrix$V12 <- runif(dim(offspring_matrix)[1])
#       for (ii in 1:dim(offspring_matrix)[1]) {
          # apply mortality if an individual is homozygous for the second loci on each chromosome
#          if (offspring_matrix[ii, "V8"] == offspring_matrix[ii, "V10"]) {
            # test to find homozygotes
#            if (offspring_matrix[ii, "V12"] < m) {
              # if the random number assigned to that ind is less than mort probability
#              offspring_matrix$V7[ii] <- 0   # set the individual to dead
#            }
#          }
#        }
        offspring_matrix[which(offspring_matrix$V8==offspring_matrix$V10 & offspring_matrix$V12<m),7]<-0
        # get the male and female number_offsprings that are alive and then select the ones we'll use for
        # the next generation
        maleoffs <- offspring_matrix[which(offspring_matrix$V1 == "Male" & offspring_matrix$V7 == 1),]
        femaleoffs <- offspring_matrix[which(offspring_matrix$V1 == "Female" & offspring_matrix$V7 == 1),]
        starting_population <- rbind(
          maleoffs[sample(1:dim(maleoffs)[1], size = population_size / 2, replace = FALSE),], 
          femaleoffs[sample(1:dim(femaleoffs)[1], size = population_size / 2, replace = FALSE),])
        starting_population <- starting_population[, c("V1", "V2", "V3","V4", "V5", "V6", "V8", "V9", "V10","V11")]
        colnames(starting_population) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10")
        pop1 <- rbind(pop1, starting_population)
        
        #calculate heterozygosity
        count_allele_A1 <- sum(length(which(starting_population$V7 == "A1")), length(which(starting_population$V9 == "A1")))
        count_allele_A2 <- (population_size * 2) - count_allele_A1
        count_allele_B1 <- sum(length(which(starting_population$V8 == "B1")), length(which(starting_population$V10 == "B1")))
        count_allele_B2 <- (population_size * 2) - count_allele_B1
        alle_freq_A1 <- count_allele_A1 / (population_size * 2)
        alle_freq_A2 <- count_allele_A2 / (population_size * 2)
        alle_freq_B1 <- count_allele_B1 / (population_size * 2)
        alle_freq_B2 <- count_allele_B2 / (population_size * 2)
        het_A <- 1 - sum((alle_freq_A1 ^ 2), (alle_freq_A2 ^ 2))
        het_B <- 1 - sum((alle_freq_B1 ^ 2), (alle_freq_B2 ^ 2))
        het_tot_A <- rbind(het_tot_A, het_A)
        het_tot_B <- rbind(het_tot_B, het_B)
        generation_tot <- rbind(generation_tot, g)
        alle_freq_tot_A1 <- rbind(alle_freq_tot_A1, alle_freq_A1)
      }
      # storing after loop generation
      colnames(pop1) <- c("sex","source","gen","malep","femalep","indnum","c1-l1","c1-l2","c2-l1","c2-l2")
      write.csv(pop1, file = paste( scenarioname,"/", scenarioname, "-cr","-mort",m * 100,"-rep-", i, ".csv",sep = ""),
        row.names = FALSE)
      write.csv(sum(het_tot_B>0),file=paste(scenarioname,"/endOfTime-",i,".csv",sep=""))
      plot(generation_tot, het_tot_A, type = "l")
      plot(generation_tot, het_tot_B, type = "l")
      plot(generation_tot, alle_freq_tot_A1,
        type = "l",
        ylab = "Frequency allele A1",
        xlab = "Number of generations",
        main = paste("Ne", population_size, "number_offspring", number_offspring, "s", m, "freq", alleles_frequencies[1])
      )
    }
  }
}
et=proc.time()[1]
message(paste("Stop the clock ",et-st))
