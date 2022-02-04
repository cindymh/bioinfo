import math
from pprint import pprint

def count(motifs):
    """
    Finds the count of corresponding nucleotides of parallel motifs

            Parameters:
                    motifs (list of strings): DNA sequences

            Returns:
                    count (dict): dictionary of nucleotides counts
    """
    count = {} # initializing the count dictionary
    k = len(motifs[0]) # motif length
    for symbol in "ACGT":
        count[symbol] = [0] * len(range(k))

    t = len(motifs) # number of motifs
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count[symbol][j] += 1

    return count

def profile(motifs):
    """
    Finds the profile matrix of list of motifs

            Parameters:
                    motifs (list of strings): DNA sequences

            Returns:
                    profile (dict): dictionary of nucleotides profile/proportion
    """
    profile = {}
    counts = count(motifs)
    t = len(motifs)

    for nuc in counts.keys():
        profile[nuc] = [val/t for val in counts[nuc]]

    return profile

def consensus(motifs):
    """
    Finds the consensus string of list of motifs

            Parameters:
                    motifs (list of strings): DNA sequences

            Returns:
                    consensus (string): the most common nucleotide appearance at
                    every column
    """
    k = len(motifs[0])
    counts = count(motifs)

    consensus = ""
    # iterate through every column
    for j in range(k):
        m = 0
        frequent_symbol = ""
        # compare frequency of symbol at each column
        for symbol in "ACGT":
            if counts[symbol][j] > m:
                m = counts[symbol][j]
                frequent_symbol = symbol
        consensus += frequent_symbol

    return consensus

def score(motifs):
    """
    Finds the score of the list of motifs

            Parameters:
                    motifs (list of strings): DNA sequences

            Returns:
                    score (int): number of nucleotides that do not match the
                    consensus string
    """
    cons = consensus(motifs)
    k = len(motifs[0])
    t = len(motifs)

    score = 0
    # adds 1 when everytime a motif element does not match consensus
    for col in range(k):
        score += sum(1 for row in range(t) if motifs[row][col] != cons[col])

    return score

def pr(text, profile):
    """
    Finds the probability of text occurring based on profile matrix

            Parameters:
                    text (string): DNA sequence
                    profile (dictionary): maps nucleotides to its probability of
                    occurrence at each location in text

            Returns:
                    pr (float): probability of k-mer text occurring
    """
    p = 1
    k = len(text)

    for i in range(k):
        breakpoint()
        p *= profile[text[i]][i]

    return p

def profile_most_probable_kmer(text, k, profile):
    """
    Finds the most probable kmer given a profile matrix

            Parameters:
                    text (string): DNA sequence
                    k (int): length of motif/k-mer
                    profile (dictionary): maps nucleotides to its probability of
                    occurrence at each location in text

            Returns:
                    most probable k-mer (string)
    """
    most_p = -1
    txt_len = len(text)

    # iterate through every possible k-mer
    for i in range(txt_len-k+1):
        kmer = text[i:i+k]
        # if current k-mer is more probable, then update most probable k-mer
        if  (p_kmer := pr(kmer, profile)) > most_p:
            most_p = p_kmer
            most_kmer = kmer

    return most_kmer

def greedy_motif_search(dna, k, t):
    """
    Finds the most probable kmer

           Parameters:
                    dna (list of string): DNA sequences
                    k (int): length of motif/k-mer
                    t (int): number of dna sequences

            Returns:
                    most probable k-mer/motif (string)
    """
    best_motifs = []
    # initializing best motifs as the first k-mer in list of DNAs
    for i in range(0, t):
        best_motifs.append(dna[i][0:k])

    n = len(dna[0])
    # iterate through every possible k-mer motif in the first DNA string
    for i in range(n-k+1):
        motifs = []
        motifs.append(dna[0][i:i+k])
        # build a profile matrix based on existing motifs
        for j in range(1, t):
            p = profile(motifs)
            # appends most probable motifs in the subsequent DNA string using 
            # profile matrix provided by 
            motifs.append(profile_most_probable_kmer(dna[j], k, p))

        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    
    return best_motifs

def entropy(motifs):
    p = profile(motifs)
    n = len(motifs[0])
    symbols = p.keys()

    h=0
    for i in range(n):
        for symbol in symbols:
            prob = p[symbol][i]

            try:
                h += prob * math.log(prob, 2)
            except ValueError:
                continue
    return -h