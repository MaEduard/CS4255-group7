from sys import stdin
import math
from node import Node
import random

tot_up_dist = 0

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


def profile_join(a, j):
    """
    Joins two profile together and returns the result. 
    """
    new_prof = make_empty_profile(len(a[0]))
    for row in range(4):  # Only 4 nucleotides.
        for col in range(len(a[0])):  # Assume same size for a and j
            new_prof[row][col] = (a[row][col] + j[row][col])/2
    return new_prof

def node_dist(i:Node, j:Node):
    """
    Calculates and returns the distance between the profiles of node a and node b.
    Including subtracting with up_distances. 
    """
    dist = 0
    prof_i = i.profile
    prof_j = j.profile
    k = len(prof_i[0])
    for a in range(k):
        for x in range(4):  # Only 4 nucleotides.
            for y in range(4):
                if x != y:
                    dist += prof_i[x][a] * prof_j[y][a]

    return dist/k - i.up_distance - j.up_distance

def prof_dist(i, j):
    """
    Calculates and returns the distance between the profiles i and j. 
    """
    dist = 0
    k = len(i[0])
    for a in range(k):
        for x in range(4):  # Only 4 nucleotides.
            for y in range(4):
                if x != y:
                    dist += i[x][a] * j[y][a]

    return dist/k

def r_equation(node_i, total_profile, n):
    ans = (n*prof_dist(node_i.profile, total_profile)) - prof_dist(node_i.profile, node_i.profile) - (n-1)*node_i.up_distance - (tot_up_dist - node_i.up_distance)
    return ans/(n-2)

def find_min_dist_nodes(nodes, total_profile, n):
    """
    Finds the minimum distance between all possible joins, and returns the index of both nodes. 

    Using exhaustive search. 

    Parameters: 
        n = number of active nodes

    TODO include the r_equation to compute r(i) and r(j) to eventually compute d'(i,j) = d(i,j) - r(i) - r(j)
    """
    min_dist = float('inf')
    ri = 0
    rj = 0
    for i in range(len(nodes)):
        if nodes[i].is_active: # Skip nodes that are inactive. 
            for j in range(i+1, len(nodes)):
                if nodes[j].is_active: # Skip nodes that are inactive. 
                    dist_i_j = node_dist(nodes[i], nodes[j])
                    r_i = out_dist(nodes[i], total_profile, n) # shouldn't we use r_equation?
                    r_j = out_dist(nodes[j], total_profile, n) # shouldn't we use r_equation?
                    
                    dist_prime_i_j = dist_i_j - r_i - r_j

                    print(dist_prime_i_j)

                    if dist_prime_i_j < min_dist:
                        min_dist = dist_prime_i_j
                        ind1 = i
                        ind2 = j

    return ind1, ind2, min_dist

def initialize_leaf_nodes(seqs) -> list:
    """ 
    Creates the initial profiles for the input sequences. 
    """
    nodes = []
    for index, seq in enumerate(seqs):
        profile = make_profile(seq)
        node = Node(profile=profile, up_distance=0, is_active=True)
        node.indexes = str(index)
        nodes.append(node)    
    return nodes

def print_profile(a):
    for row in a:
        print(row)

def set_total_up_dist(val):
    global tot_up_dist
    tot_up_dist = val

def compute_total_profile(nodes, active_nodes):
    """
    Computes the profile matrix of all active nodes.

    Also sets the global variable for the total up distance. 

    Args:
        nodes: all the nodes of the tree
        active_nodes: the number of active nodes in the tree. 

    Returns:
        float[]: profile matrix of size 4 x k
    """
    k = len(nodes[0].profile[0])
    profile = make_empty_profile(k)
    dist = 0
    for node in nodes:
        if node.is_active:
            dist += node.up_distance
            for row in range(4):
                for col in range(k):
                    profile[row][col] += node.profile[row][col]/active_nodes

    set_total_up_dist(dist)
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

def out_dist(node, total_profile, n):
    """
    Calculates the out distance from one node. 
    """
    # return n*(prof_dist(node.profile, total_profile) + tot_up_dist - node.up_distance) / (n-2) 
    return (n*(prof_dist(node.profile, total_profile)) - prof_dist(node.profile, node.profile) - tot_up_dist - (n-1)*node.up_distance) / (n-2) 

def calc_branch_len(i:Node, j:Node, n, total_profile) -> tuple[float, float]:
    """
    Calculates the branch length from daughter nodes i and j. 

        Parameters:
            i = node i 
            j = node j
            n = number of active nodes
            total_profile = distance to all other nodes
            tot_up_dist = total up distance
    """
    dist_i_j = node_dist(i, j)

    r_i = out_dist(i, total_profile, n)
    r_j = out_dist(j, total_profile, n)
    res_i = (dist_i_j + r_i - r_j)/2
    res_j = (dist_i_j + r_j - r_i)/2

    return res_i, res_j

def find_joins(nodes:list, seqs):
    """
    Loops through all active nodes and joins the nodes with the minimum distance. 
    Stops when the number of active nodes is just two. 
    """

    branch_lengths = []
    active_nodes = len(nodes)
    global total_profile
    total_profile = compute_total_profile(nodes, active_nodes)

    seed_nodes = [node for node in nodes]
    # create_top_hits(seed_nodes, nodes, len(nodes))
    while active_nodes > 2:
        total_profile = compute_total_profile(nodes, active_nodes)
        i, j, mindist = find_min_dist_nodes(nodes, total_profile, active_nodes)
        # print(mindist)

        print(nodes[i].indexes, "<-->", nodes[j].indexes)
        print("----------------")

        prof_i = nodes[i].profile
        prof_j = nodes[j].profile

        nodes[i].is_active = False 
        nodes[j].is_active = False

        k = profile_join(prof_i, prof_j)

        up_dist = prof_dist(nodes[i].profile, nodes[j].profile)/2 # node_dist(nodes[i], nodes[j])/2 
        new_node = Node(profile=k, up_distance=up_dist, is_active=True)
        nodes.append(new_node)

        new_node.indexes = "("+ nodes[i].indexes + "," + nodes[j].indexes + ")"

        
        i_len, j_len = calc_branch_len(nodes[i],nodes[j], active_nodes, total_profile)
        branch_lengths.append(i_len)
        branch_lengths.append(j_len)
        active_nodes-=1
        # break
    
    return(branch_lengths)

# def create_top_hits(seed_nodes, nodes, n):
#     m = int(math.sqrt(len(nodes)))

#     while(len(seed_nodes) > 0):
#         seed_idx = random.randint(0, len(seed_nodes)-1) #pick random seed node. Note that randint picks random int from inclusive range, hence do len(seed_nodes) - 1
#         seed = seed_nodes[seed_idx]
#         top_hits = []
#         r_seed = out_dist(seed, total_profile, n)

#         for node in nodes:
#             if node != seed:
#                 dist_i_j = node_dist(node, seed)
#                 r_j = out_dist(node, total_profile, n) # shouldn't we use r_equation?
#                 dist_prime_i_j = dist_i_j - r_seed - r_j
#                 top_hits.append((dist_prime_i_j, node))

#         top_hits.sort()

#         nodes[seed_idx].top_hits = top_hits[0:m]
#         seed_nodes.remove(seed)

#         print(seed_nodes)

#         for neighbor in top_hits[0:m]:
#             neighbor_node = neighbor[1]
#             top_hits_neighbor = []
#             for node_tuple in top_hits[0:2*m]:
#                 node = node_tuple[1]
#                 if node != neighbor and len(node.top_hits) == 0:
#                     dist_i_j = node_dist(neighbor_node, node)
#                     r_j = out_dist(neighbor_node, total_profile, n) 
#                     r_i = out_dist(node, total_profile, n)
#                     dist_prime_i_j = dist_i_j - r_i - r_j
#                     top_hits_neighbor.append((dist_prime_i_j, node))
#             top_hits_neighbor.sort()
#             print("Top hits length:")
#             print(len(top_hits))
#             neighbor_node.top_hits = top_hits_neighbor
#             print(neighbor_node)
#             seed_nodes.remove(neighbor_node)


def main():
    seqs = read_file('data/test-small.aln')
    print(seqs)
    nodes = initialize_leaf_nodes(seqs)

    ### PROOF that equation on the right bottom of page 3 in old paper equals comparison to 'total profile'. 
    # equation is needed for computing r(i) and r(j)
    # i = make_profile(seqs[0])
    # print(find_total_profile(i, nodes))
    # print(prof_dist(i, compute_total_profile(seqs))*len(seqs))
    ###
    
    branch_lengths = find_joins(nodes, seqs)
    print(branch_lengths)

if __name__ == '__main__':
    main()