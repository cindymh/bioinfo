def pattern_matching(genome, pattern):
    """
    Finds the position of a pattern present in a genome string.

            Parameters:
                    genome (str): string of AGCT
                    pattern (str): a specific pattern to be matched in genome

            Returns:
                    idx (list of int): list of indices where pattern starts in
                    genome
    """
    idx = []
    for i in range(len(genome)-len(pattern)+1):
        if genome[i:i+len(pattern)] == pattern:
            idx.append(i)
    return idx