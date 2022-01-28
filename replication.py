from frequent_words import pattern_count

def symbol_array(genome, symbol):
    """
    Finds the symbol count along a sliding window of length half the genome

            Parameters:
                    genome (str): DNA string
                    symbol (str): pattern in DNA string

            Returns:
                    array (int): count of symbol along windows of length genome//2
    """
    array = {}
    n = len(genome)
    Extendedgenome = genome + genome[0:n//2]

    # look at the first half of genome to compute first array value
    array[0] = pattern_count(genome[0:n//2], symbol)

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if Extendedgenome[i-1] == symbol:
            array[i] = array[i]-1
        if Extendedgenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def skew_array(genome):
    """
    This function produces a skew array which indicates direction of increasing
    or decreasing value #G - #C as it traverses through a genome

            Parameters:
                    genome (str): DNA string

            Returns:
                    skew (list of int): value of #G-#C at each position in the
                    genome
    """

    skew = [0] * (len(genome) + 1)

    for idx, nuc in enumerate(genome):
        if nuc == 'C':
            skew[idx+1] = skew[idx] - 1
        elif nuc == 'G':
            skew[idx+1] = skew[idx] + 1
        else:
            skew[idx+1] = skew[idx] 

    return skew

def minimum_skew(genome):
    """
    This function produces a list of positions in which skew array gives minimum
    values

            Parameters:
                    genome (str): DNA string

            Returns:
                    list (int): indices of minimum values in an array
    """
    # set a variable equal to SkewArray(Genome)
    skew = skew_array(genome)
    # find the minimum value of all values in the skew array
    min_val = min(skew)

    # range over the length of the skew array and add all positions achieving the min to positions
    return [idx for idx, val in enumerate(skew) if val == min_val]

def hamming_distance(p, q):
    """
    This function returns the hamming distance between two genome sequences

            Parameters:
                    p (str): DNA string
                    q (str): DNA string

            Returns:
                    count (int): hamming distance between the two sequences
    """


    count = 0
    # increments counts when there are differences between the two genomes
    for nuc1, nuc2 in zip(p, q):
        if nuc1 != nuc2:
            count += 1

    return count

    # list comprehension
    #return sum([1 for (n1, n2) in zip(p, q) if n1 != n2])

def approximate_pattern_matching(genome, pattern, d):
    """
    This function returns the location where pattern approximately appears in
    genome

            Parameters:
                    genome (str): DNA string
                    pattern (str): pattern DNA string
                    d (int): hamming distance

            Returns:
                    locations (list): list of indices where pattern approximately
                    occurred in genome
    """
    pattern_length = len(pattern)
    locations = []
    # look for positions in genome where pattern approximately matches
    for i in range(0, len(genome)-pattern_length+1):
        if(hamming_distance(genome[i:i+pattern_length], pattern) <= d):
            locations.append(i)
    return locations

def approximate_pattern_count(pattern, genome, d):
    """
    This function approximates the number of times a pattern occurs in a genome
    using hamming distance

            Parameters:
                    genome (str): DNA string
                    pattern (str): pattern DNA string
                    d (int): hamming distance

            Returns:
                    int: approximate number of times a pattern occurs in genome
    """
    return len(approximate_pattern_matching(genome, pattern, d))