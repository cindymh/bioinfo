def pattern_count(text, pattern):
    """
    Counts the occurrence of pattern in a text.

            Parameters:
                    text (str)
                    pattern (str): potential substring within text

            Returns:
                    count (int): number of pattern occured in text
    """
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if text[i:i+len(pattern)] == pattern:
            count = count+1
    return count 

def frequency_map(text, k):
    """
    Creates pattern-count pairs in a dictionary.

            Parameters:
                    text (str)
                    k (int): length of patterns to be extracted from text

            Returns:
                    freq (dict): pattern-count pairs/frequency map
    """
    freq = {}
    n = len(text)
    for i in range(n-k+1):
        pattern = text[i:i+k]
        freq[pattern] = 0

    # hint: your code goes here!
    for pattern in freq.keys():
        freq[pattern] = pattern_count(text, pattern)

    return freq

def frequent_words(text, k):
    """
    Finds the most frequent k-mer pattern(s) which occurred in text

            Parameters:
                    text (str)
                    k (int): length of patterns to be extracted from text

            Returns:
                    words (list): the most frequent k-letter pattern(s)/word(s)
    """
    words = []
    freq = frequency_map(text, k)
    m = max(freq.values())
    for key in freq:
        # add each key to words whose corresponding frequency value is equal to m
        if freq[key] == m:
            words.append(key)
    return words