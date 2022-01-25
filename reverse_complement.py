from reverse import reverse
from complement import complement

def reverse_complement(pattern):
    """
    Finds the complement of a DNA string.

            Parameters:
                    pattern (str): string which complement is to be found

            Returns:
                    complement pattern (str): complement DNA pattern
    """
    return complement(reverse(pattern))

print(reverse_complement("CCAGATC"))