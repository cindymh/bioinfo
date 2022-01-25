from frequency_map import frequency_map

"""
Finds the most frequent k-mer pattern(s) which occurred in text

        Parameters:
                text (str)
                k (int): length of patterns to be extracted from text

        Returns:
                words (list): the most frequent k-letter pattern(s)/word(s)
"""

def frequent_words(text, k):
    words = []
    freq = frequency_map(text, k)
    m = max(freq.values())
    for key in freq:
        # add each key to words whose corresponding frequency value is equal to m
        if freq[key] == m:
            words.append(key)
    return words