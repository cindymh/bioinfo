from pattern_count import pattern_count

def frequency_map(text, k):
    freq = {}
    n = len(text)
    for i in range(n-k+1):
        pattern = text[i:i+k]
        freq[pattern] = 0

    # hint: your code goes here!
    for pattern in freq.keys():
        freq[pattern] = pattern_count(text, pattern)

    return freq