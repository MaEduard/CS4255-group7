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

def r_equation(node_i, total_profile, n, up_sum):
    ans = n*prof_dist(node_i.profile, total_profile) - prof_dist(node_i, node_i) - (n-1)*node_i.up_distance + node_i.up_distance - up_sum
    return ans/(n-2)

def find_min_dist_nodes(nodes, total_profile):
    """
    Finds the minimum distance between all possible joins, and returns the index of both nodes. 

    Using exhaustive search. 

    TODO include the r_equation to compute r(i) and r(j) to eventually compute d'(i,j) = d(i,j) - r(i) - r(j)
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
        node = Node(profile=profile, up_distance=0, is_active=True)
        node.add_sequence(seq)
        nodes.append(node)    
    return nodes

def print_profile(a):
    for row in a:
        print(row)

def compute_total_profile(seqs):
    """Computes the profile matrix of a seqs matrix
    Args:
        seqs (str[]): matrix of dna sequences

    Returns:
        float[]: profile matrix of size 4 x k
    """
    profile = make_empty_profile(len(seqs[0]))
    for column in range(len(seqs[0])):
        count = [0,0,0,0]
        # go over every motif per column
        for row in range(len(seqs)):
            if seqs[row][column] == 'A':
                count[0] += 1
            elif seqs[row][column] == 'C':
                count[1] += 1
            elif seqs[row][column] == 'G':
                count[2] += 1
            else:
                count[3] +=1
        for profile_row in range(4):
            profile[profile_row][column] = count[profile_row]/(len(seqs))
    return profile

def find_total_profile(prof_i, nodes):
    """Computes the total distance of node i to other nodes in the old fashion way. We should not do this but
    actually compute this value by comparing node i to the total profile of all active nodes.
    
    * I created this function to double check whether the formula on the right bottom of page 3 of old paper actually holds.
    Args:
        prof_i (float[][]): 
        nodes (str[][]): 

    Returns:
        float: 
    """
    total_profile = 0
    for node in nodes:
        if node.is_active:
            total_profile += prof_dist(prof_i, node.profile)
    return total_profile

def find_joins(nodes:list, active_nodes, seqs):
    """
    Loops through all active nodes and joins the nodes with the minimum distance. 
    Stops when the number of active nodes is just one. 

    TODO Maybe removing the nodes from the list can be used instead of using a bool flag?
    TODO Double check what equation to use for up_dist. In the supplementary they give a different equation?
    """

    total_profile = compute_total_profile(seqs)
    while active_nodes > 1:
        a, b, mindist = find_min_dist_nodes(nodes, total_profile)
        nodes[a].is_active = False 
        nodes[b].is_active = False
        prof_a = nodes[a].profile
        prof_b = nodes[b].profile
        c = profile_join(prof_a,prof_b)
        # up_dist = (prof_dist(prof_a,c) + prof_dist(prof_b,c))/2 
        up_dist = prof_dist(prof_a, prof_b)/2
        new_node = Node(profile=c, up_distance=up_dist, is_active=True)
        for seq in nodes[a].sequences:
            new_node.add_sequence(seq)
        for seq in nodes[b].sequences:
            new_node.add_sequence(seq)
        nodes.append(new_node)
        active_nodes -=1
    

def main():
    seqs = read_file('./data/test-small.aln')
    nodes = initialize_leaf_nodes(seqs)
    active_nodes = len(nodes)

    ### PROOF that equation on the right bottom of page 3 in old paper equals comparison to 'total profile'. 
    # equation is needed for computing r(i) and r(j)
    i = make_profile(seqs[0])
    print(find_total_profile(i, nodes))
    print(prof_dist(i, compute_total_profile(seqs))*len(seqs))
    ###
    
    find_joins(nodes, active_nodes, seqs)

if __name__ == '__main__':
    main()