import numpy as np
import itertools
import re


def checkDuplicates(geneList):
    '''Returns False if their are multiple of the same allele in the np.array'''
    return np.array([len(np.unique(np.char.lower(geneList[i]))) == len(geneList[i]) for i in range(len(geneList))])

    
def deleteDuplicates(geneList,findAlleles):
    '''Returns np.array of alleles np.arrays w/o duplicates'''
    newGeneList=[]
    for i in range(len(geneList)):
        if findAlleles[i] == True:
            newGeneList.append(geneList[i])
    return np.array(newGeneList)

        
def determineChildGenotypes(alleles1, alleles2, nHyb1):
    '''Returns all possible genotypes of the child and their likelihoods'''
    
    tempJ=[]
    tempFlattenJ=[]
    allChild = []
    for i in range(len(alleles1)):
        for j in range(len(alleles2)):
            for k in range(nHyb1):
                tempJ.append(np.array(sorted(np.array([alleles1[i],alleles2[j]])[:,k])))
            tempFlattenJ = np.array(tempJ).flatten()
            allChild.append('-'.join(re.findall('..',''.join(tempFlattenJ))))
            tempFlattenJ=[]
            tempJ=[]
  
    childGenDict = {i:(allChild.count(i))/len(allChild) for i in allChild}
        
    return childGenDict
        

def punnettSquare(gen1=0,gen2=0):
    '''Returns all possible genotypes of the child'''
    if gen1==0:
        genotype1 = input("Please enter genotype of parent 1 (eg. rr-yy-Ww-ss): ") # Parent 1
        genotype2 = input("Please enter genotype of parent 2: ") # Parent 2
    else:
        genotype1 = gen1
        genotype2 = gen2
    
    
    # Determine punnett square alleles for parent 1  
    genotype1ListTemp=np.array(genotype1.split('-'))
    nHybrid1 = genotype1ListTemp.size
    genotype1List=np.array([np.array(list(genotype1ListTemp[i])) for i in range(nHybrid1)])
    
    genotype1CombList = np.array(list(itertools.combinations(genotype1List.flatten(),nHybrid1)))
    genotype1FindAlleles = checkDuplicates(genotype1CombList)
    genotype1AllelesList = deleteDuplicates(genotype1CombList,genotype1FindAlleles)
    
    
    # Determine punnett square alleles for parent 2 
    genotype2ListTemp=np.array(genotype2.split('-'))
    nHybrid2 = genotype2ListTemp.size
    genotype2List=np.array([np.array(list(genotype2ListTemp[i])) for i in range(nHybrid2)])
    
    genotype2CombList = np.array(list(itertools.combinations(genotype2List.flatten(),nHybrid2)))
    genotype2FindAlleles = checkDuplicates(genotype2CombList)
    genotype2AllelesList = deleteDuplicates(genotype2CombList,genotype2FindAlleles)
    

    childGenAndProb = determineChildGenotypes(genotype1AllelesList, genotype2AllelesList, nHybrid1)
    
    print("The potential genotypes and corresponding probabilities of the child are:")
    for key, value in childGenAndProb.items():
         print('{} ({:.2%})'.format(key, value))


# Can either enter parent genotypes after running code or add them as arguments to the function below
# eg. punnettSquare('rr-yy-Ww-ss','RR-yy-WW-Ss')

punnettSquare()