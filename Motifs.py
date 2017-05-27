#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 16:07:26 2017

@author: igor
"""

# Input:  A set of kmers Motifs
# Output: Count(Motifs)
def Count(Motifs):
    count = {} # initializing the count dictionary
    # Filling the dictionary with adequate number of zeroes
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {}
    profile = Count(Motifs)
    for l in profile:
        for i in range(k):
            profile[l][i] = profile[l][i]/t
    return profile

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus  

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    consensus = Consensus(Motifs)
    score = 0
    m = len(Motifs)
    n = len(Motifs[0])
    for j in range(n):
        for i in range(m):
            if Motifs[i][j] != consensus[j]:
                score += 1
            else:
                score = score
    return score

# Input:  String Text and profile matrix Profile (where Profile is a dictionary)
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
   p = 1
   for i in range(len(Text)):
        p *= Profile[Text[i]][i]
   return p

# Input:  String Text, an integer k, and profile matrix Profile
# Output: A Profile-most probable k-mer in Text,
# i.e. ProfileMostProbablePattern(Text, k, Profile)
def ProfileMostProbablePattern(Text, k, Profile):
    most_prob = Text[0:k] 
    p_max = Pr(Text[0:k], Profile)
    for i in range(1, len(Text) - k + 1):
         if Pr(Text[i:i+k], Profile) > p_max:
                p_max = Pr(Text[i:i+k], Profile)
                most_prob = Text[i:i+k]        
    return most_prob       
                                   
# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    #Initialise best_motifs as a list featuring the 
    #first k-mer from each dna string:
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
        
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {} # initializing the count dictionary
    # Filling the dictionary with adequate number of zeroes
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(0)
    # filling in the dictionary with pseudocunts
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    #add 1 to every position to obtain pseudo count
    for symbol in 'ACGT':
        for j in range(k):
            count[symbol][j] += 1
    return count

# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {} # output variable
    # Filling the dictionary with adequate number of zeroes
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
             profile[symbol].append(0)
    count = CountWithPseudocounts(Motifs)
    sum = 0
    for i in 'ACGT':
        sum = sum + count[i][0]
    for i in 'ACGT':
        for j in range(k):
            profile[i][j] = count[i][j]/sum
    
    return profile

