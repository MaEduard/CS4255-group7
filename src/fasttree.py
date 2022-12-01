from sys import stdin
import math


def read_file(file_name):
    f = open(file_name, "r")
    sequences = []
    for line in f.readlines():
        if len(line) > 4:  # hard-coded value just so we can skip the lines that have only a number
            sequences.append(line[:len(line)-1])  # to remove the \n

    return sequences


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
    profile = [[], [], [], []]
    for row in profile:
        for i in range(size):
            row.append(0)

    return profile


def profile_join(a, b):
    # Size of first row from a (sequence length)
    new_prof = make_empty_profile(len(a[0]))
    for col in range(len(a)):  # Assume same size for a and b
        for row in range(4):  # Only 4 nucleotides.
            new_prof[row][col] = (a[row][col] + b[row][col])/2


def prof_dist(a, b):
    """
    Calculates and returns the distance between two profiles.
    """
    dist = 0
    for i in range(len(a)):
        for x in range(4):  # Only 4 nucleotides.
            for y in range(4):
                if x != y:
                    dist += a[x][i] * b[y][i]

    return dist/len(a)
<<<<<<< HEAD


# seqs = []
# for line in stdin:
#     seqs.append()

# profiles = []
# profile = []
# nucs = ["A", "C", "G", "T"]
# for i in range(len(seqs)):
#     profile = []
#     for x in seqs[i]:
#         freq = [0, 0, 0, 0]
#         pos = nucs.index(x)
#         freq[pos] = 1
#         profile.append(freq)
#     profiles.append(profile)

# print(prof_dist(profiles[0], profiles[1]))
=======
    
seqs = []
for line in stdin:
    line = line.strip()
    seqs.append(line)

profiles = []
for seq in seqs:
    profile = make_profile(seq)
    profiles.append(profile)

print(profiles[0])
>>>>>>> afedf5a364fa212c6d2d0a8b6feb0f45da445bb0


class Node:
    profile = []
    up_distance = 0
    is_active = True
    children = []


def main():
    read_file('data/test-small.aln')


if __name__ == '__main__':
    main()
