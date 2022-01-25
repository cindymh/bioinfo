"""
Counts the occurrence of pattern in a text.

        Parameters:
                text (str)
                pattern (str): potential substring within text

        Returns:
                count (int): number of pattern occured in text
"""
def pattern_count(text, pattern):
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if text[i:i+len(pattern)] == pattern:
            count = count+1
    return count 