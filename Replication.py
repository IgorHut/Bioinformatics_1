# Input:  Strings Pattern and Text
# Output: The number of times Pattern appears in Text
def PatternCount(Pattern, Text):
    count = 0 # output variable
    # your code here
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count


# Input:  A string Text and an integer k
# Output: CountDict(Text, k)
def CountDict(Text, k):
    Count = {}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i] = PatternCount(Pattern, Text)
    return Count


# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    FrequentPatterns = []
    Count = CountDict(Text, k)
    m = max(Count.values())
    for i in Count:
        if Count[i] == m:
            FrequentPatterns.append(Text[i:i+k])
    FrequentPatternsNoDuplicates = remove_duplicates(FrequentPatterns)
    return FrequentPatternsNoDuplicates

# Input:  A list Items
# Output: A list containing all objects from Items without duplicates
def remove_duplicates(Items):
    ItemsNoDuplicates = [] # output variable

    for i in Items:
        if i not in ItemsNoDuplicates:
            ItemsNoDuplicates.append(i)

    return ItemsNoDuplicates

# Input:  A DNA string Pattern
# Output: The reverse complement of Pattern
def ReverseComplement(Pattern):
    revComp = '' # output variable

    for i in range(1, len(Pattern)+1):
        revComp = revComp + Pattern[-i]

    revComp = complement(revComp)
    return revComp

# Input:  A character Nucleotide
# Output: The complement of Nucleotide
def complement(Nucleotide):

    comp = '' # output variable

    for i in Nucleotide:
        if i == 'A':
            comp = comp + 'T'
        elif i == 'C':
            comp = comp + 'G'
        elif i == 'G':
            comp = comp + 'C'
        elif i == 'T':
            comp = comp + 'A'

        else:
           print ('Something is wrong, check whether the sequence cotains only valid nucleotides, i.e. ACGT!')
           break

    return comp

# Input:  Two strings, Pattern and Genome
# Output: A list containing all starting positions where Pattern appears as a substring of Genome
def PatternMatching(Pattern, Genome):
    positions = [] # output variable
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions

#########################################################################
# I need this to print out the solution for Bioinformatics Hacker Track
# The code is in Python 3
#data = PatternMatching(Pattern, Genome)
#print(*data, sep=' ')
#########################################################################
# Compute the Hamming distance between two strings p and q
# Input:  Two strings p and q of the same length
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    hamm_dist = 0
    for index, a in enumerate(p):
        if a!= q[index]:
            hamm_dist = hamm_dist + 1
        else:
            hamm_dist = hamm_dist
    return hamm_dist


# Input:  Strings Genome and symbol
# Output: SymbolArray(Genome, symbol)
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

# Input:  Strings Genome and symbol
# Output: FasterSymbolArray(Genome, symbol)
def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])
    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i] - 1
        elif ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i] + 1
    return array

# Input:  A String Genome
# Output: Skew(Genome)
def Skew(Genome):
    skew = {} #initializing the dictionary
    skew[0] = 0
    for i, k in enumerate(Genome):
        i = i+1
        skew[i] = skew[i-1]
        if k == 'G':
            skew[i] = skew[i-1] + 1
        elif k == 'C':
            skew[i] = skew[i] -1
    return skew

# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Genome):
    positions = [] # output variable
    min_skew = min(Skew(Genome).values())
    for i, k in Skew(Genome).items():
        if k == min_skew:
            positions.append(i)
    # your code here
    return positions

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Pattern, Text, d):
    positions = [] # initializing list of positions
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Pattern, Text[i:i+len(Pattern)]) <=d:
            positions.append(i)
    return positions

# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Pattern, Text, d):
    count = 0 # initialize count variable
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d:
            count = count+1
    return count