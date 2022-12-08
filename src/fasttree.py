from sys import stdin
import math
from node import Node


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
    """
    Creates an empty profile filled with zeros. 
    """
    profile = [[], [], [], []]
    for row in profile:
        for i in range(size):
            row.append(0)

    return profile


def profile_join(a, b):
    """
    Joins two profile together and returns the result. 
    """
    new_prof = make_empty_profile(len(a[0]))
    for col in range(len(a)):  # Assume same size for a and b
        for row in range(4):  # Only 4 nucleotides.
            new_prof[row][col] = (a[row][col] + b[row][col])/2
    return new_prof

def prof_dist(a, b):
    """
    Calculates and returns the distance between two profiles.
    """
    dist = 0
    for i in range(len(a[0])):
        for x in range(4):  # Only 4 nucleotides.
            for y in range(4):
                if x != y:
                    dist += a[x][i] * b[y][i]

    return dist/len(a[0])

def find_min_dist_nodes(nodes):
    """
    Finds the minimum distance between all possible joins, and returns the index of both nodes. 

    Using exhaustive search. 
    """
    min_dist = float('inf')
    for i in range(len(nodes)):
        if not nodes[i].is_active: # Skip nodes that are inactive. 
            continue
        for j in range(i+1, len(nodes)):
            if not nodes[j].is_active: # Skip nodes that are inactive. 
                continue
            dist = prof_dist(nodes[i].profile, nodes[j].profile) - nodes[i].up_distance - nodes[j].up_distance
            if dist < min_dist:
                min_dist = dist
                ind1 = i
                ind2 = j
    return ind1, ind2, min_dist

def initialize_leaf_nodes(seqs) -> list:
    """ 
    Creates the initial profiles for the input sequences. 
    """
    nodes = []
    for seq in seqs:
        profile = make_profile(seq)
        node = Node(profile, up_distance=0, is_active=True)
        nodes.append(node)    
    return nodes

def print_profile(a):
    for row in a:
        print(row)

def find_joins(nodes:list, active_nodes):
    """
    Loops through all active nodes and joins the nodes with the minimum distance. 
    Stops when the number of active nodes is just one. 

    TODO Maybe removing the nodes from the list can be used instead of using a bool flag?
    TODO Node.is_active is currently not working
    """
    while active_nodes > 1:
        a, b, mindist = find_min_dist_nodes(nodes)
        nodes[a].is_active = False 
        nodes[b].is_active = False
        prof_a = nodes[a].profile
        prof_b = nodes[b].profile
        c = profile_join(prof_a,prof_b)
        up_dist = (prof_dist(prof_a,c) + prof_dist(prof_b,c))/2
        nodes.append(Node(profile=c, up_distance=up_dist, is_active=True))

        active_nodes -=1

seqs = read_file('../data/test-small.aln')
nodes = initialize_leaf_nodes(seqs)
active_nodes = len(nodes)
print(nodes[0].is_active)

find_joins(nodes, active_nodes)

def main():
    read_file('../data/test-small.aln')


if __name__ == '__main__':
    main()
