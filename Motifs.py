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
