# This function returns counts of the occurence of k-mer pattern in a text
def pattern_count(text, pattern):
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if text[i:i+len(pattern)] == pattern:
            count = count+1
    return count 