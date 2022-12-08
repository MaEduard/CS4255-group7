from sys import stdin
import math


def letter_to_index(let):
    if let == "A":
        return 0 
    elif let == "C":
        return 1
    elif let == "G":
        return 2
    elif let == "T":
        return 3

def jukes_cantor(prof):

    dist = -(3/4)*math.log(1-(4/3))

def total_profile():
    pass


def make_profile(seq):
    """
    Makes a new profile based on one sequence. 
    """
    prof = make_empty_profile(len(seq))
    for i, let in enumerate(seq):
        prof[letter_to_index(let)][i] = 1
    return prof


def make_empty_profile(size):
    profile = [[],[],[],[]]
    for row in profile:
        for i in range(size):
            row.append(0)
    
    return profile

def profile_join(a, b):
    new_prof = make_empty_profile(len(a[0])) # Size of first row from a (sequence length)
    for col in range(len(a)): # Assume same size for a and b
        for row in range(4): # Only 4 nucleotides.
            new_prof[row][col] = (a[row][col] + b[row][col])/2

def prof_dist(a, b):
    """
    Calculates and returns the distance between two profiles.
    """
    dist = 0
    for i in range(len(a)):
        for x in range(4): # Only 4 nucleotides.
            for y in range(4):
                if x !=y:
                    dist += a[x][i] * b[y][i]

    return dist/len(a)

def find_min_dist(profile, profiles):
    dist = 0
    min_dist = float('inf')
    for prof in profiles:
        dist = prof_dist(profile, prof)
        if dist == 0:
            continue 
        if dist < min_dist:
            min_dist = dist 
            


seqs = []
for line in stdin:
    line = line.strip()
    seqs.append(line)

profiles = []
for seq in seqs:
    profile = make_profile(seq)
    profiles.append(profile)

print(profiles[0])


class Node:
    profile = []
    up_distance = 0
    is_active = True
    children = []