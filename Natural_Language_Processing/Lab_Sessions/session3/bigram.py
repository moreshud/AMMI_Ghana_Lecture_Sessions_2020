import io, sys, math, re
from collections import defaultdict
import numpy as np

# GOAL: build a stupid backoff bigram model

def load_data(filename):
    fin = io.open(filename, 'r', encoding='utf-8')
    data = []
    vocab = defaultdict(lambda:0)
    for line in fin:
        sentence = line.split()
        data.append(sentence)
        for word in sentence:
            vocab[word] += 1
    return data, vocab


def remove_rare_words(data, vocab, mincount=0):
    ## FILL CODE
    # replace words in data that are not in the vocab 
    # or have a count that is below mincount
    return data_with_unk


def build_bigram(data):
    unigram_counts = defaultdict(lambda:0)
    bigram_counts  = defaultdict(lambda: defaultdict(lambda: 0.0))
    total_number_words = 0

    ## FILL CODE
    # Store the unigram and bigram counts as well as the total 
    # number of words in the dataset

    unigram_prob = defaultdict(lambda:0)
    bigram_prob = defaultdict(lambda: defaultdict(lambda: 0.0))

    ## FILL CODE
    # Build unigram and bigram probabilities from counts

    return {'bigram': bigram_prob, 'unigram': unigram_prob}

def get_prob(model, w1, w2):
    assert model["unigram"][w2] != 0, "Out of Vocabulary word!"
    ## FILL CODE
    # Should return the probability of the bigram (w1w2) if it exists
    # Else it return the probility of unigram (w2) multiply by 0.4

def perplexity(model, data):
    ## FILL CODE
    # follow the formula in the slides
    # call the function get_prob to get P(w2 | w1)
    return perp

def generate(model):
    sentence = ["<s>"]
    ## FILL CODE
    # generate a sentence. A sentence starts with a <s> and ends with a </s>
    # Possiblly a use function is:
    #   np.random.choice(x, 1, p = y)
    # where x is a list of things to sample from
    # and y is a list of probability (of the same length as x)
    return sentence

###### MAIN #######

print("load training set")
train_data, vocab = load_data("train2.txt")
## FILL CODE 
# If you have a Out of Vocabulary error (OOV) 
# call the function "remove_rare_words" to replace 
# rare words with <unk> in the dataset

print("build bigram model")
model = build_bigram(train_data)

print("load validation set")
valid_data, _ = load_data("valid2.txt")
## FILL CODE 
# If you have a Out of Vocabulary error (OOV) 
# call the function "remove_rare_words" to replace 
# OOV with <unk> in the dataset

print("The perplexity is", perplexity(model, valid_data))

print("Generated sentence: ",generate(model))
