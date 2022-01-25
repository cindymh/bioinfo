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