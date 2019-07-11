#########################################
#Büşra YAŞAR-150114063
#Emine Feyza MEMİŞ-150114077
#Hale ŞAHİN-150116841
#########################################
import random

def indexToBase(i):
    return "ACGT"[i] #return a represent base for values 0->A, 1->C, 2->G, 3->T

def baseToIndex(a):
    return "ACGT".index(a) #return a represent value for bases A->0, C->1, G->2, T->3

def readInput( file, t ):
    dna = []
    i=0
    for i in range(t):
        line=file.readline().strip() #keep dna strings in dna array
        mutatedline=mutation(str(line))
        dna.append(mutatedline)
    return dna
##################################### Randomized Motif Search #########################################################
def randomizedMotifSearch( k, dna ):
    bestMotifs = randomMotifs(dna, k)  #choose first motifs as random
    bestScore = Score(bestMotifs,k) #find best score of motifs
    while True:
        profile = profileFromMotifs(bestMotifs,k,1) #calulate profile matrix
        motifs = motifsFromProfile(dna,profile) #find motifs from dna strings by profile matrix
        score = Score(motifs,k) #find score of new motifs
        if score < bestScore: #keep lowest score as best score
            bestMotifs = motifs
            bestScore = score
        else:
            return bestMotifs #return best motifs

def motifsFromProfile( dna, profile ):
    return [bestKmerForProfile(seq, profile) for seq in dna] #send dna seq to find best motif for each seq

def bestKmerForProfile( seq, profile ):
    k = len(profile)
    bestProb = -1
    bestKmer = ''
    # new motif calculation
    for start in range(len(seq)-k+1): #search dna 0 to 500-k+1
        kmer = seq[start:start+k] #take string with k size
        prob = probKmerInProfile(kmer, profile) #calculate probability of kmer by using profile matrix
        if prob > bestProb: #choose kmer with higher probability for the seq dna
            bestProb = prob
            bestKmer = kmer
    return bestKmer

def probKmerInProfile( kmer, profile ): #probability of kmer
    prob = 1.0
    for i in range(len(kmer)): #loop k times
        prob *= profile[i][baseToIndex(kmer[i])] #take ith base of kmer and find its prob from profile,
    return prob                                  #multiply probs of each base in the kmer to find total prob


def mut(base): #do the mutation
    if base=='A': #if related base is A then it can mutated into C or G or T
        rnd = random.randint(0, 2) #select a random value
        if rnd==0: #if zero mutate into C
            mutated='C'
        elif rnd==1: #if 1 mutate into G
            mutated='G'
        elif rnd==2: #if 2 mutate into T
            mutated='T'
    elif base=='C': #if related base is C then it can mutated into A or G or T
        rnd = random.randint(0, 2)
        if rnd==0:
            mutated='A'
        elif rnd==1:
            mutated='G'
        elif rnd==2:
            mutated='T'
    elif base=='G': #if related base is C then it can mutated into A or C or T
        rnd = random.randint(0, 2)
        if rnd==0:
            mutated='A'
        elif rnd==1:
            mutated='C'
        elif rnd==2:
            mutated='T'
    elif base=='T': #if related base is C then it can mutated into A or C or G
        rnd = random.randint(0, 2)
        if rnd==0:
            mutated='A'
        elif rnd==1:
            mutated='C'
        elif rnd==2:
            mutated='G'
    return mutated

def mutation(line):
    rndpos = random.randint(0, 490)  # select a random position in dna string
    tenmer=line[rndpos:rndpos + 10]  # take tenmer from that position
    l=list(range(10)) # to get random mutation points in ten-mer we used 0 to 9 numbers list
    random.shuffle(l) # we shuffled the list to get first four of the list as random positions different than each other
    j=0 #is for to keep tennmer position
    mutatedtenmer=[]
    #print(rndpos)
    #print(tenmer)
    for i in tenmer: #loop elements of tenmer
        if j==l[0] or j==l[1] or j==l[2] or j==l[3]: #if we arrive mut points do the mutation
            mutatedtenmer.append(mut(i)) #take mutated point from return value of mut function and append it to mutated string tenmer
            #print(str(j),i,mutbas)
        else:
            mutatedtenmer.append(i) #if there is no mutation append that base directly to the new tenmer
        j=j+1
    muttenmer = ''.join(mutatedtenmer) #convert array into string
    #print (muttenmer)
    mutatedLine=line[0:rndpos]+muttenmer+line[rndpos+10:len(line)] #create new mutated dna string
    #print (mutatedLine)
    return mutatedLine

def Score( motifs,k ):
    score = 0
    for count in countsFromMotifs(motifs,k,0): #find ACGT counts of motifs
        score += sum(count) - max(count) #calculate score by taking max count as a reference
    return score
############################################Gibbs Sampler #############################################################
def GibbsSampler(k,t,dna): #make iterative randomized motif search until condition is satisfied
    s_best = float('inf') #give starting value  as big number
    curr_motifs = randomMotifs(dna,k)  # choose first motifs as random
    motifs_best = []
    i=0
    while True:
        x = random.randint(0, t - 1)  # select a random x point from motifs
        curr_motifs.pop(x)  # take xth motif off from motifs
        profile = profileFromMotifs(curr_motifs,k,1)  # calculate profile matrix after deleting xth motif
        motif_select = ProfileRandomGenerator(profile, dna, k, x)  # select random motif from the xth dna with using probs
        curr_motifs.insert(x, motif_select) #insert selected motif into current motifs
        s_curr = Score(curr_motifs,k) #find score of motifs
        if s_curr < s_best: #if the new score is less than(better than) best score we kept, assign it as best
            motifs_best = curr_motifs.copy() #copy current motifs to best motifs
            s_best = s_curr
            i = 0
        else:
            i = i + 1  # if best score didn't change count i increased by 1
        if i > 150:  # loop until last 50 best score doesn't change,then break when there is no longer improve
            break
    return [s_best,motifs_best]

def ProfileRandomGenerator(profile, dna, k, index):
    probs_list = []
    for x in range(len(dna[index]) - k + 1):#loop over the chosen dna
        probability = 1 #give prob starting value 1 to multiply probs
        ex_kmer = dna[index][x : k + x]  #take sequentially kmers from the chosen dna
        for y in range(k): #loop kmer with index y
            probability *= profile[y][baseToIndex(ex_kmer[y])] #prob of ex_kmer occur is calculated
        probs_list.append(probability) #keep a score list for all kmer probabilities
   # rnd = random.uniform(0, sum(probs_list)) #take a random value in range sum of the all probs
    #current = 0
    randomm=0
    for z, bias in enumerate(probs_list): #loop for finding random value interval in the probs list
        if bias > randomm:
            randomm=bias
            point=z
    return dna[index][point: k + point]
     #enumerate probs starting from 0 increase by 1, z keeps index bias keeps prob value
     #   current += bias
     #   if rnd <= current: #when find random number lied interval current
     #       return dna[index][z : k + z] #return  current kmer starting from z point

########################################### Common functions ###########################################################
def randomMotifs( dna, k ):
    return [randomKmer(seq,k) for seq in dna] #find random motifs in dna strings for the start

def randomKmer( seq, k ):
    num=int(len(seq)-k)
    #print (num)
    start = random.randint(0, num) #select a random position in dna string
    return seq[start:start+k] #return k sized string motif

def countsFromMotifs( motifs, k, initCount ): #find counts of motifs
    counts = []
    for i in range(k): #loop k times
        currCount = [initCount] * 4 #make each value in the count array 0
        for motif in motifs:
            currCount[baseToIndex(motif[i])] += 1 #find base counts as for A currCount[0]=how many times A exist
        counts.append(currCount)                  #currCount[1] is how many times C exist in motifs ith location of kmer
    return counts

def profileFromMotifs( motifs,k, initCount ): #profile matrix calculation
    counts = countsFromMotifs(motifs,k,initCount) #take base counts from motifs
    profile = []
    for count in counts:
        total = float(sum(count)) #find total sum of counts and cast it to float to calculate prob as float
        probs = [c/total for c in count] #find count/total count as profile write it into c value in count array
        profile.append(probs) #add probs to profile array
    return profile #return profile values of given motifs


def findConsensus(bestMotifs,k):
    profile = profileFromMotifs(bestMotifs,k, 1)  # calculate profile matrix of best motifs
    maxProb=0.0
    consensus=[]
    for i in range(k):  # loop k times
        for j in range(4):
            if(profile[i][j]>maxProb): #take the base with max probability from the profile
                maxProb=profile[i][j]
                base=j
        maxProb=0.0
        consensus.append(indexToBase(base)) # convert base index to base,and add chosen base to the consensus
    strconsensus = ''.join(consensus) #convert consensus to string
    return strconsensus



print("Please enter the file path: " ) #take dna file name as an input
file = input()
#file= "C://Users//feyza//Desktop//input.txt"
print("please enter the k value") #take k-mer size value as input
k = int(input())
#k=9
t = 10 #number of input dna strings
dna = readInput(open(file), t) #call read file method
bestMotifs1 = randomizedMotifSearch(k, dna) #call randomized motif search for dna array from the input file and k value
print('Results of Randomized Motif Search')
print('###Best Motifs###')
for motif in bestMotifs1: #to print best motifs at the end of execution
    print (motif)
sc1=Score(bestMotifs1,k)
consensus1=findConsensus(bestMotifs1,k)   #call method to calculate consensus according to the best motifs
print('Consensus String:    '+consensus1)
print('Score:               ', str(sc1)) #print consensus
#########################Gibbs#########################################################################################
print('Results of Gibbs Sampler')
bestMotifs2 = GibbsSampler(k, t, dna) #call gibbs sampler for dna array from the input file and k value
print('###Best Motifs###')
for motif in bestMotifs2[1]: #to print best motifs at the end of execution
    print (motif)
sc2=bestMotifs2[0]
consensus2=findConsensus(bestMotifs2[1],k)   #call method to calculate consensus according to the best motifs
print('Consensus String:    '+consensus2)
print('Score:               ', str(sc2)) #print consensus
