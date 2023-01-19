import random
from sys import stdin
import math
from node import Node

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

def make_empty_profile(size):
    """
    Creates an empty profile filled with zeros. 
    """
    profile = [[], [], [], []]
    for row in profile:
        for i in range(size):
            row.append(0)

    return profile

def make_profile(seq):
    """
    Makes a new profile based on one sequence. 
    """
    prof = make_empty_profile(len(seq))
    for i, let in enumerate(seq):
        prof[letter_to_index(let)][i] = 1
    return prof

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

def profile_join(a, j):
    """
    Joins two profile together and returns the result. 
    """
    new_prof = make_empty_profile(len(a[0]))
    for row in range(4):  # Only 4 nucleotides.
        for col in range(len(a[0])):  # Assume same size for a and j
            new_prof[row][col] = (a[row][col] + j[row][col])/2
    return new_prof

def out_dist(node, total_profile, n):
    """
    Calculates the out distance from one node. 
    """
    # return n*(prof_dist(node.profile, total_profile) + tot_up_dist - node.up_distance) / (n-2) 
    return (n*(prof_dist(node.profile, total_profile)) - prof_dist(node.profile, node.profile) - tot_up_dist - (n-1)*node.up_distance) / (n-2) 

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

                    if dist_prime_i_j < min_dist:
                        min_dist = dist_prime_i_j
                        ind1 = i
                        ind2 = j

    return ind1, ind2, min_dist

def find_nodes_to_be_joined(nodes):
    """
    Finds the minimum distance between all possible joins, and returns the index of both nodes. 

    Using exhaustive search. 

    Parameters: 
        n = number of active nodes

    TODO include the r_equation to compute r(i) and r(j) to eventually compute d'(i,j) = d(i,j) - r(i) - r(j)
    """
    min_dist = float('inf')
    node1 = None
    node2 = None
    ri = 0
    rj = 0
    for i in range(len(nodes)):
        if nodes[i].is_active: # Skip nodes that are inactive. 
            current_dist = nodes[i].top_hits[0][0] 
            if current_dist < min_dist:
                min_dist = current_dist
                node1 = nodes[i]
                node2 = nodes[i].top_hits[0][2]

    return node1, node2 

def merge_top_hits(top_hits1, top_hits2):
    top_hits_merged = top_hits1 + top_hits2
    top_hits_merged = list(set(top_hits_merged))

    return top_hits_merged

def filter_top_hits(merged_node, top_hits, active_nodes):
    new_top_hits = []
    for node_tuple in top_hits:
        node = node_tuple[2]
        dist_i_j = node_dist(merged_node, node)
        r_j = out_dist(merged_node, total_profile, active_nodes) 
        r_i = out_dist(node, total_profile, active_nodes) 
        dist_prime_i_j = dist_i_j - r_i - r_j 
        new_top_hits.append((dist_prime_i_j, node.index, node))
    
    new_top_hits.sort()
    return new_top_hits[0:top_hits_length]


def join_nodes(node1, node2, active_nodes):
    # create new profile
    # create new top hits list
    # set name of new node
    # create new up_distance
    # update total_profile?
    # set node1 and node2 to inactive

    new_profile = profile_join(node1.profile, node2.profile)
    new_updistance = prof_dist(node1.profile, node2.profile)
    merged_node = Node(profile=new_profile, up_distance=new_updistance, is_active=True)
    top_hits_merged = merge_top_hits(node1.top_hits, node2.top_hits)
    filtered_top_hits = filter_top_hits(merged_node, top_hits_merged, active_nodes)
    merged_node.top_hits = filtered_top_hits

    node1.is_active = False
    node2.is_active = False

    return merged_node

def initialize_leaf_nodes(seqs) -> list:
    """ 
    Creates the initial profiles for the input sequences. 
    """
    nodes = []
    for index, seq in enumerate(seqs):
        profile = make_profile(seq)
        node = Node(profile=profile, up_distance=0, is_active=True)
        node.value = str(index)
        node.index = index
        node.main_sequence = seq
        nodes.append(node)    
    return nodes

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
    global top_hits_length
    global node_counter
    node_counter = active_nodes-1
    top_hits_length = int(math.sqrt(len(nodes)))
    total_profile = compute_total_profile(nodes, active_nodes)
    create_top_hits(nodes, len(nodes), total_profile)
    while active_nodes > 2:
        total_profile = compute_total_profile(nodes, active_nodes)
        node1, node2 = find_nodes_to_be_joined(nodes)
        node12 = join_nodes(node1, node2, active_nodes)
        # i, j, mindist = find_min_dist_nodes(nodes, total_profile, active_nodes)

        print(node1.value, "<-->", node2.value)
        print("----------------")

        # prof_i = nodes[i].profile
        # prof_j = nodes[j].profile

        # nodes[i].is_active = False 
        # nodes[j].is_active = False

        # k = profile_join(prof_i, prof_j)

        # up_dist = prof_dist(nodes[i].profile, nodes[j].profile)/2 # node_dist(nodes[i], nodes[j])/2 
        # new_node = Node(profile=k, up_distance=up_dist, is_active=True)
        nodes.append(node12)

        node12.value = "("+ node1.value + "," + node2.value + ")"

        ## REDO THIS? 
        # i_len, j_len = calc_branch_len(nodes[i],nodes[j], active_nodes, total_profile)
        # branch_lengths.append(i_len)
        # branch_lengths.append(j_len)
        active_nodes-=1
        # break
    return(branch_lengths)

def create_top_hits(nodes, n, total_profile):
    m = int(math.sqrt(len(nodes))) # maybe use math.ceil? 
    unused_nodes = [i for i in range(len(nodes))]

    while(len(unused_nodes) > 0):
        seed_idx = random.sample(unused_nodes, 1)[0] #pick random seed node. Note that randint picks random int from inclusive range, hence do len(seed_nodes) - 1
        seed = nodes[seed_idx]
        top_hits = []
        r_seed = out_dist(seed, total_profile, n)

        for node in nodes:
            if node.main_sequence != seed.main_sequence:
                dist_i_j = node_dist(node, seed)
                r_j = out_dist(node, total_profile, n) # shouldn't we use r_equation?
                dist_prime_i_j = dist_i_j - r_seed - r_j
                top_hits.append((dist_prime_i_j, node.index, node))

        top_hits.sort() # the distances for some top hits are apparently exactly the same and thus it starts to compare nodes, which is impossible. 

        nodes[seed_idx].top_hits = top_hits[0:m]
        unused_nodes.remove(seed_idx)

        for neighbor in top_hits[0:m]:
            neighbor_node = neighbor[2]
            if len(neighbor_node.top_hits) != 0:
                continue
            top_hits_neighbor = []
            for node_tuple in top_hits[0:2*m]:
                node = node_tuple[2]
                if node.main_sequence != neighbor_node.main_sequence:
                    dist_i_j = node_dist(neighbor_node, node) 
                    r_j = out_dist(neighbor_node, total_profile, n) 
                    r_i = out_dist(node, total_profile, n) 
                    dist_prime_i_j = dist_i_j - r_i - r_j 
                    top_hits_neighbor.append((dist_prime_i_j, node.value, node)) # add node.value to make sure that if two nodes have exactly the same distance_prime_i_j, pick one randomly
                    
            top_hits_neighbor.sort()
            neighbor_node.top_hits = top_hits_neighbor[0:m]
            if neighbor_node in unused_nodes:
                unused_nodes.remove(neighbor_node)

def main():
    seqs = read_file('data/test-small.aln')


    ### TODO ###
    # perhaps change the top_hits list to contain the index of the node instad of the node object itself.

    nodes = initialize_leaf_nodes(seqs)

    branch_lengths = find_joins(nodes, seqs)
    print(branch_lengths)

if __name__ == '__main__':
    main()