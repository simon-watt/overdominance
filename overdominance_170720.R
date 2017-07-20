# script to simulate overdominance without migration...

# Richmond change...
#setwd("C:/Users/z5122213/Desktop/r_overdominance")
setwd("~/R")

# how many loci are being modeled

nloci<-2

# how many metapopulations
meta<-1

# how many generations are we simulating
gens<-100

# what is the total pop size
popsize<-50

# what years are individuals being transferred

transfer<-400 # setting this to a high number in order to skip transfers...

# how many individuals are being transferred per transfer (currently restricted to 1 or 2)

numtrans<-0 # not doing transfers

# how many offspring does each pairing produce?
pairings<-10

# does recombination happen in males? TRUE = yes, FALSE = no
mrecomb<-TRUE

# how many replicates to run
replicates<-1

# file name info...
scenarioname<-"scenario1"
dir.create(scenarioname)

# need to loop over a range of crossing over rates...
for(c in seq(0,0,0)){
 
 message("c = ",c)
  
  # need to loop over a range of mortality rates due to homozygosity
  for(m in seq(0.15,0.15)){
    message("m = ",m)

    # need to loop through a series of replicate runs...
    for(i in 1:replicates){
      message("i = ",i)
      locilist<-list()
      allelefreqp1<-list()
      # first need to generate a base population...
      # Name the alleles at loci 1 (the potentially lethal one...)
      locilist[[1]]<-c(1,2)
      # Allele frequencies for loci 1
      allelefreqp1[[1]]<-c(0.2,0.8)
      # Name the alleles at loci 2 (the neutral one...)
      locilist[[2]]<-c(10,20)
      # Allele frequencies for loci 2
      allelefreqp1[[2]]<-c(0.5,0.5)

      # Calculating the number of individuals/meta population
      metasz<-c(popsize/meta)
  
      # Create a data.frame to track details...
      popgenome1<-as.data.frame(matrix(ncol=nloci*2+6,nrow=metasz))
  
      # first column is sex
      popgenome1[1:(metasz/2),1]<-"Male"
      popgenome1[(metasz/2+1):metasz,1]<-"Female"
  
      # second column is population
      popgenome1[,2]<-1
  
      # third column is generation
      popgenome1[,3]<-0
  
      # fourth column is male parent
      popgenome1[,4]<-0
  
      # fifth column is female parent
      popgenome1[,5]<-0
  
      # sixth column is generation individual number
      popgenome1[,6]<-rownames(popgenome1)
  
      
      # generate initial genotypes for population 1
      
      # going to put the two loci on a single stretch of chromosome one after the other in the data table
      for(k in 1:nloci){
        popgenome1[1:metasz,(k+6)]<-sample(locilist[[k]], size=metasz, replace=TRUE, prob=allelefreqp1[[k]])
        popgenome1[1:metasz,(k+8)]<-sample(locilist[[k]], size=metasz, replace=TRUE, prob=allelefreqp1[[k]])
      }
  
      # We have our population genomes, now we have to simulate the population over time...
  
      pop1<-popgenome1
     
      #heterozygosity variables
      het_tot_1=list()
      het_tot_2=list()
      generation_tot=list()
      alle_freq__tot_1=list()
  
      for (g in 1:gens){
        #message("g = ",g)
        # males and females from pop 1
        malep1<-which(popgenome1$V1=="Male")
        femalep1<-which(popgenome1$V1=="Female")
    
        # create matrix into which parent combinations are placed
        parents1<-matrix(nrow=metasz,ncol=2)
    
        # select males and females that are mating with each other - males column 1 females column 2
        parents1[,1]<-sample(as.numeric(malep1), size=metasz, replace=TRUE, prob=NULL)
        parents1[,2]<-sample(as.numeric(femalep1), size=metasz, replace=TRUE, prob=NULL)
        
        # have parents, now lets make offspring making pairings*offspring/pairing young
        offspring1<-as.data.frame(matrix(nrow=metasz*pairings,ncol=(7+nloci*2)))
        
        # now have to loop through the pairings to make offspring...
        for(p in 1:dim(parents1)[1]){
          # give the offspring a sex
          offspring1[((p-1)*pairings+1):(p*pairings),1]<-sample(c("Male","Female"),size=pairings, replace=TRUE)
          # mark the source population of the individual
          offspring1[((p-1)*pairings+1):(p*pairings),2]<-1
          # mark the generation of the offspring
          offspring1[((p-1)*pairings+1):(p*pairings),3]<-g
          # give the male parent
          offspring1[((p-1)*pairings+1):(p*pairings),4]<-parents1[p,1]
          # give the female parent
          offspring1[((p-1)*pairings+1):(p*pairings),5]<-parents1[p,2]
          
          # 6 is reserved for ind # and 7 is reserved for alive/dead status
          
          # getting genotype of each offspring for a pairing...
          
          # parent chromosomes
          mchrome1<-c(popgenome1[parents1[p,1],7],popgenome1[parents1[p,1],8])
          mchrome2<-c(popgenome1[parents1[p,1],9],popgenome1[parents1[p,1],10])
          fchrome1<-c(popgenome1[parents1[p,2],7],popgenome1[parents1[p,2],8])
          fchrome2<-c(popgenome1[parents1[p,2],9],popgenome1[parents1[p,2],10])
          for(n in 1:pairings){
            # getting the genotype from dad
            randad<-runif(1) # getting a random deviate to see if crossing over happens
            mctemp1<-rep(NA,2)
            mctemp2<-rep(NA,2)
            if(randad<c & mrecomb){
            mctemp1<-c(mchrome1[1],mchrome2[2])
            mctemp2<-c(mchrome2[1],mchrome1[2])
            } else {
              mctemp1<-mchrome1 # store 1st male gamete from dad
              mctemp2<-mchrome2 # store 2nd male gamete from dad
            }
            malechromes<-list(mctemp1,mctemp2)
            # pick which male gamete this offspring is getting...
            offspring1[(p-1)*pairings+n,8:9]<-malechromes[[sample(c(1,2),1)]]
            
            # getting the genotype from mom
            ranmom<-runif(1) # getting a random deviate to see if crossing over happens
            fctemp1<-rep(NA,2)
            fctemp2<-rep(NA,2)
            if(ranmom<c){
             fctemp1<-c(fchrome1[1],fchrome2[2])
            fctemp2<-c(fchrome2[1],fchrome1[2])
            } else {
              fctemp1<-fchrome1 # store 1st female gamete from mom
              fctemp2<-fchrome2 # store 2nd female gamete from mom
            }
            femalechromes<-list(fctemp1,fctemp2)
            # pick which male gamete this offspring is getting...
            offspring1[(p-1)*pairings+n,10:11]<-femalechromes[[sample(c(1,2),1)]]
          }
        }
        offspring1[,6]<-as.numeric(rownames(offspring1))
        offspring1[,7]<-1 # "1 = alive"
        offspring1$V12<-runif(dim(offspring1)[1])
        for(ii in 1:dim(offspring1)[1]){
          # apply mortality if an individual is homozygous for the second loci on each chromosome
          if(offspring1[ii,"V8"]==offspring1[ii,"V10"]){ # test to find homozygotes
            if(offspring1[ii,"V12"]<m){  # if the random number assigned to that ind is less than mort probability
              offspring1$V7[ii]<-0   # set the individual to dead
            }
          }
        }
        # get the male and female offsprings that are alive and then select the ones we'll use for 
        # the next generations
        maleoffs<-offspring1[which(offspring1$V1=="Male" & offspring1$V7==1),]
        femaleoffs<-offspring1[which(offspring1$V1=="Female" & offspring1$V7==1),]
        popgenome1<-rbind(maleoffs[sample(1:dim(maleoffs)[1],size=metasz/2,replace=FALSE),],femaleoffs[sample(1:dim(femaleoffs)[1],size=metasz/2,replace=FALSE),])
        popgenome1<-popgenome1[,c("V1","V2","V3","V4","V5","V6","V8","V9","V10","V11")]
        colnames(popgenome1)<-c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10")
        pop1<-rbind(pop1,popgenome1)
        
        #calculate heterozygosity
        count_allele_10<-sum(length(which(popgenome1$V8==10)),length(which(popgenome1$V10==10))) 
        count_allele_20<-100-count_allele_10 
        het_1<-1-sum(((count_allele_10/100)^2),((count_allele_20/100)^2))
        count_allele_1<-sum(length(which(popgenome1$V7==1)),length(which(popgenome1$V9==1)))
        count_allele_2<-100-count_allele_1
        het_2<-1-sum(((count_allele_1/100)^2),((count_allele_2/100)^2))
        het_tot_1<-rbind(het_tot_1,het_1)
        het_tot_2<-rbind(het_tot_2,het_2)
        generation_tot<-rbind(generation_tot,g)
        # frequency allele 1
        alle_freq_1<-(count_allele_1*(sum(count_allele_1,count_allele_2))/10000)
        alle_freq__tot_1<-rbind(alle_freq__tot_1, alle_freq_1)
      }
      colnames(pop1)<-c("sex","source","gen","malep","femalep","indnum","c1-l1","c1-l2","c2-l1","c2-l2")
      write.csv(pop1,file=paste(scenarioname,"\\",scenarioname,"-cr","-mort",m*100,"-rep-",i,".csv",sep=""),row.names=FALSE)
      #plot(generation_tot,het_tot_1)
      plot(generation_tot,het_tot_2)
      plot(generation_tot, alle_freq__tot_1,type="l", ylab = "Allele frequency allele 1",xlab = "Number of generations",main = "case: No recombination, m=0.15 ")
    }
  }
}

