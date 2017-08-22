import numpy as np
import os
import time

pid=os.getpid()
print("pid = {}".format(pid))

# set random seed so each run different

np.random.seed(pid)

popsize=100
popsize_half=round(popsize/2)
noff=5
ngen=200
mort_A1=0
mort_A2=0
nloci=2
male_recomb=False
reps=1
recomb_rate=0

A1=1
A2=2
B1=10
B2=20


alle_names=list()
alle_names.append([A1,A2])
alle_names.append([B1,B2])
alle_freq=list()
alle_freq.append([0.5,0.5])
alle_freq.append([0.2,0.8])

het_tot_A=list()
het_tot_B=list()
gen_tot=list()
alle_freq_tot_A1=list()
alle_freq_tot_A2=list()

def onePop():
    starting_pop=np.zeros(shape=(popsize,2*nloci+6+1),dtype=np.int) 
    # extra col to coincide with R

    starting_pop[0:popsize_half,1]=0 # males
    starting_pop[popsize_half:popsize,1]=1 # females
    starting_pop[:,2]=1 # population
    starting_pop[:,3]=0 # initial generation
    starting_pop[:,4]=0 # male parent
    starting_pop[:,5]=0 # female parent
    starting_pop[:,6]=range(popsize) # id

    for k in range(nloci):
            starting_pop[:,6+k+1]=np.random.choice(alle_names[k],p=alle_freq[k],size=popsize,replace=True)
            starting_pop[:,8+k+1]=np.random.choice(alle_names[k],p=alle_freq[k],size=popsize,replace=True)

    nextgen=starting_pop

    g=0
    while (True):
        g=g+1
        male_parents=np.where(starting_pop[:,1]==0)[0]
        female_parents=np.where(starting_pop[:,1]==1)[0]
        parent_matrix=np.zeros(shape=(popsize,2+1),dtype=np.int)
        parent_matrix[:,1]=np.random.choice(male_parents,size=popsize,replace=True)
        parent_matrix[:,2]=np.random.choice(female_parents,size=popsize,replace=True)
        offspring_matrix=np.zeros(shape=(popsize*noff,8+2*nloci+1),dtype=np.int)
        for p in range(popsize):
            offspring_matrix[p*noff:(p+1)*noff,1]=np.random.choice([0,1],size=noff,replace=True)
            offspring_matrix[p*noff:(p+1)*noff,2]=1
            offspring_matrix[p*noff:(p+1)*noff,3]=parent_matrix[p,1]
            offspring_matrix[p*noff:(p+1)*noff,4]=parent_matrix[p,2]
        
            mchrome1=list((starting_pop[parent_matrix[p,1],7],starting_pop[parent_matrix[p,1],8]))
            mchrome2=list((starting_pop[parent_matrix[p,1],9],starting_pop[parent_matrix[p,1],10]))
            fchrome1=list((starting_pop[parent_matrix[p,2],7],starting_pop[parent_matrix[p,2],8]))
            fchrome2=list((starting_pop[parent_matrix[p,2],9],starting_pop[parent_matrix[p,2],10]))
        
            for n in range(noff):
                randad=np.random.rand(1)
                if randad<recomb_rate and male_recomb:
                    mctemp1=list((mchrome1[0],mchrome2[1]))
                    mctemp2=list((mchrome2[0],mchrome1[1]))
                else:
                    mctemp1=list((mchrome1[0],mchrome1[1]))
                    mctemp2=list((mchrome2[0],mchrome2[1]))
                malechrome=list((mctemp1,mctemp2))
                offspring_matrix[p*noff+n,8:10]=malechrome[np.random.choice([0,1])]
            
                ranmom=np.random.rand(1)
                if ranmom<recomb_rate:
                    fctemp1=list((fchrome1[0],fchrome2[1]))
                    fctemp2=list((fchrome2[0],fchrome1[1]))
                else:
                    fctemp1=list((fchrome1[0],fchrome1[1]))
                    fctemp2=list((fchrome2[0],fchrome2[1]))
                femalechrome=list((fctemp1,fctemp2))
                offspring_matrix[p*noff+n,10:12]=femalechrome[np.random.choice([0,1])]
                
        offspring_matrix[:,6]=range(popsize*noff) # id
        offspring_matrix[:,7]=1 # alive
        
        offspring_matrix[:,12]=np.random.rand(popsize*noff)<mort_A1
        deadA1=(offspring_matrix[:,8]==A1) & (offspring_matrix[:,10]==A1) & offspring_matrix[:,12]
        offspring_matrix[deadA1,7]=0
        offspring_matrix[:,12]=np.random.rand(popsize*noff)<mort_A2
        deadA2=(offspring_matrix[:,8]==A2) & (offspring_matrix[:,10]==A2) & offspring_matrix[:,12]
        offspring_matrix[deadA2,7]=0 
        
        male_off=np.where((offspring_matrix[:,1]==0) & offspring_matrix[:,7])[0]
        female_off=np.where((offspring_matrix[:,1]==1) & offspring_matrix[:,7])[0]
    
        collist=[0,1,2,3,4,5,6,8,9,10,11]
        starting_pop[0:popsize_half]=offspring_matrix[np.random.choice(male_off,size=popsize_half)][:,collist]    
        starting_pop[popsize_half:]=offspring_matrix[np.random.choice(female_off,size=popsize_half)][:,collist]
    
        nextgen=np.vstack((nextgen,starting_pop))
    
        count_A1=sum(starting_pop[:,7]==A1)+sum(starting_pop[:,9]==A1)
        count_A2=sum(starting_pop[:,7]==A2)+sum(starting_pop[:,9]==A2)
        count_B1=sum(starting_pop[:,8]==B1)+sum(starting_pop[:,10]==B1)
        count_B2=sum(starting_pop[:,8]==B2)+sum(starting_pop[:,10]==B2)
        alle_freq_A1=count_A1/(2.0*popsize)
        alle_freq_A2=count_A2/(2.0*popsize)
        alle_freq_B1=count_B1/(2.0*popsize)
        alle_freq_B2=count_B2/(2.0*popsize)
        het_A=1-(np.power(alle_freq_A1,2)+np.power(alle_freq_A2,2))
        het_B=1-(np.power(alle_freq_B1,2)+np.power(alle_freq_B2,2))
        het_tot_A.append(het_A)
        het_tot_B.append(het_B)
        gen_tot.append(g)
        alle_freq_tot_A1.append(alle_freq_A1)
        alle_freq_tot_A2.append(alle_freq_A2)
        if alle_freq_A1==0:
            return ('A1',g)

        if alle_freq_A2==0:
            return ('A2',g) 
   
    
fname1='fix-py-A1-'+str(pid)+'.dat'
fname2='fix-py-A2-'+str(pid)+'.dat' 
fh1=open(fname1,'w')
fh2=open(fname2,'w')
nreps=10
st=time.clock()
for rep in range(nreps):
    ret=onePop()
    print('No more {} after generation {}'.format(ret[0],ret[1]))
    if ret[0]=='A1':
        fh1.write(str(ret[1])+'\n')
    else:
        fh2.write(str(ret[1])+'\n')
et=time.clock()

print("Stop the clock {}".format(et-st))

fh1.close()
fh2.close()
