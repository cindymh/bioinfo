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

def reverse(pattern):
    """
    Reverses a string.

            Parameters:
                    pattern (str): string to be reversed

            Returns:
                    rev (str): reversed string
    """
    rev = ''
    for char in pattern:
        rev = char + rev
    # or just return pattern[::-1]
    return rev

def reverse_complement(pattern):
    """
    Finds the reverse complement of a DNA string/sequence.

            Parameters:
                    pattern (str): string which complement is to be found

            Returns:
                    complement pattern (str): reverse complement DNA pattern
    """
    return complement(reverse(pattern))