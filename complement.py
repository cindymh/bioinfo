def complement(pattern):
    """
    Finds the complement of a DNA string.

            Parameters:
                    pattern (str): string which complement is to be found

            Returns:
                    complement pattern (str): complement DNA pattern
    """
    # dictionary of DNA complements
    comp = {
        'A': 'T', 
        'T': 'A', 
        'C': 'G', 
        'G': 'C'
    }
    # creates translation table to create complements of DNA pattern
    table = str.maketrans(comp)
    return pattern.translate(table)
    # can use comp.get(char in pattern)