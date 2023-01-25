# import math
import math
from node import Node
import random
import numpy

tot_up_dist = 0
tot_prof_size = 0


def read_file(file_name):
    """

    Reads the input file and appends the sequences.

        Parameters:
            file_name: the name of the input file.

        Returns:
            The sequences from the input file.
    """
    f = open(file_name, "r")
    sequences = []
    real_line = False
    for line in f.readlines():
        if real_line:
            sequences.append(line[:len(line)-1])  # to remove the \n
        real_line = not real_line

    return sequences


def letter_to_index(let):
    """
    Maps each nucleotide base to an index.

        Returns:
            The index of the given nucleotide base.
    """
    if let == "A":
        return 0
    elif let == "C":
        return 1
    elif let == "G":
        return 2
    elif let == "T":
        return 3


def jukes_cantor(d_u):

    dist = -(3/4)*math.log(1-(4/3)*d_u)

    return dist


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

        Returns:
            profile: A profile filled with zeros.
    """
    profile = [[], [], [], []]
    for row in profile:
        for i in range(size):
            row.append(0)

    return profile


def profile_join(a, b):
    """
    Joins two profile together and returns the resulting profile.

        Parameters:
            a,b: the profiles to be joined.
        Returns:
            the resulting profile.
    """
    new_prof = make_empty_profile(len(a[0]))
    for row in range(4):  # Only 4 nucleotides.
        for col in range(len(a[0])):  # Assume same size for a and j
            new_prof[row][col] = (a[row][col] + b[row][col])/2
    return new_prof


def prof_dist(i, j):
    """
    Calculates and returns the distance between the profiles i and j.

        returns:
            The profile distance between i and j.
    """
    dist = 0
    k = len(i[0])
    for a in range(k):
        for x in range(4):  # Only 4 nucleotides.
            for y in range(4):
                if x != y:
                    dist += i[x][a] * j[y][a]
    return dist/k


def find_min_dist_nodes(nodes, total_profile, n):
    """
    Finds the minimum distance between all possible joins,
    and returns the index of both nodes.

    Using exhaustive search.

    Parameters:
        n = number of active nodes
    """
    min_dist = float('inf')
    for i in range(len(nodes)):
        if nodes[i].is_active:  # Only look at active nodes.
            for j in range(i+1, len(nodes)):
                if nodes[j].is_active:  # Only look at active nodes.
                    dist_i_j = calc_d(nodes[i], nodes[j])
                    r_i = out_dist(nodes[i], total_profile, n)
                    r_j = out_dist(nodes[j], total_profile, n)

                    # r_i = dist_to_all_other_nodes(i, nodes)/(n-2)
                    # r_j = dist_to_all_other_nodes(j, nodes)/(n-2)

                    # print("r_i : ", r_i, "distance to all other nodes: ", r_i2)

                    # The more negative, the better (in general)
                    dist_prime_i_j = dist_i_j - r_i - r_j

                    if dist_prime_i_j < min_dist:
                        min_dist = dist_prime_i_j
                        ind1 = i
                        ind2 = j

    return ind1, ind2, min_dist


def initialize_leaf_nodes(seqs) -> list:
    """
    Creates the initial profiles for the input sequences.

        Parameters:
            seqs: the list of sequences to be analyzed.

        Returns:
            nodes: each sequence turned into a leaf node.
    """
    nodes = []
    for index, seq in enumerate(seqs):
        profile = make_profile(seq)
        node = Node(profile=profile, up_distance=0, is_active=True)
        node.value = str(index)
        node.index = index
        nodes.append(node)

    return nodes


def print_profile(a):
    """
    Used for testing.
    """
    dist = 0
    for col in range(len(a[0])):
        for row in range(len(a)):
            dist += a[row][col]
        print(dist)
        dist = 0


def set_total_up_dist(val):
    """
    Sets the global variable for total up distance.
    """
    global tot_up_dist
    tot_up_dist = val


def update_total_profile(i: Node, j: Node, k: Node, total_profile, n):
    """
    Updates the total profile by subtracting the joined nodes,
    and adding the newly created node. More efficient than iterating
    over all nodes after each join.

    TODO check if this returns the correct total profile
    """

    new_up = tot_up_dist - i.up_distance - j.up_distance + k.up_distance
    set_total_up_dist(new_up)
    iprof = i.profile
    jprof = j.profile
    kprof = k.profile

    for row in range(4):
        for col in range(len(i.profile[0])):  # Number of columns in a profile.
            total_profile[row][col] += kprof[row][col]/(n-1)\
                - iprof[row][col]/n-jprof[row][col]/n

    return total_profile


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


def out_dist(node, total_profile, n):
    """
    Calculates the out distance from one node.
    """
    return (n*(prof_dist(node.profile, total_profile)) -
            prof_dist(node.profile, node.profile) -
            (n-1)*node.up_distance + node.up_distance - tot_up_dist) / (n-2)


def calc_branch_len(i: Node, j: Node, n, total_profile):
    """
    Calculates the branch length from daughter nodes i and j.

        Parameters:
            i = node i
            j = node j
            n = number of active nodes
            total_profile = distance to all other nodes
        Returns:
            res_i, res_j: the branch lengths for the parent nodes.
    """
    dist_i_j = calc_d(i, j)

    r_i = out_dist(i, total_profile, n)
    r_j = out_dist(j, total_profile, n)
    res_i = (dist_i_j + r_i - r_j)/2
    res_j = (dist_i_j + r_j - r_i)/2

    if res_i < 0:
        res_j -= res_i
        res_i = 0
    elif res_j < 0:
        res_i -= res_j
        res_j = 0

    return round(res_i, 3), round(res_j, 3)


def calc_branch_len_without_totprof(i, j, n, nodes):
    """
    Calculates the branch length from daughter nodes i and j.

        Parameters:
            i = node i
            j = node j
            n = number of active nodes
            total_profile = distance to all other nodes
        Returns:
            res_i, res_j: the branch lengths for the parent nodes.
    """
    dist_i_j = calc_d(nodes[i], nodes[j])

    r_i = dist_to_all_other_nodes(i, nodes)/(n-2)
    r_j = dist_to_all_other_nodes(j, nodes)/(n-2)
    res_i = (dist_i_j + r_i - r_j)/2
    res_j = (dist_i_j + r_j - r_i)/2

    if res_i < 0:
        res_j -= res_i
        res_i = 0
    elif res_j < 0:
        res_i -= res_j
        res_j = 0

    return round(res_i, 3), round(res_j, 3)


def recalculate_profiles(node):
    """ After a nearest neighbor interchange, this method recalculates the profiles
    of the nodes influenced by the change by propagating it up to the root of the tree.
    This is done recursively. """
    if node.parent == None:
        return
    new_profile = profile_join(node.left.profile, node.right.profile)
    node.profile = new_profile

    recalculate_profiles(node.parent)


def nearest_neighbor_interchanges(root):
    print("----------NNI-----------------")
    queue = []
    stack = []

    queue.append(root)

    while len(queue) > 0:
        node = queue.pop(0)
        if not node.left == None:
            queue.append(node.left)
            stack.append(node.left)
        if not node.right == None:
            queue.append(node.right)
            stack.append(node.right)

    while len(stack) > 0:
        node = stack.pop()

        # find nodes that have both of their children not leaves. Following the illustration about NNI in Fig. 1 in the paper
        if node.has_two_children():
            if node.left.has_two_children() and node.right.has_two_children():
                node_a = node.left.left
                node_b = node.left.right
                node_c = node.right.left
                node_d = node.right.right

                # Use log-corrected profile distance (without up-distance)
                d_a_b = jukes_cantor(prof_dist(node_a.profile, node_b.profile))
                d_c_d = jukes_cantor(prof_dist(node_c.profile, node_d.profile))
                d_a_c = jukes_cantor(prof_dist(node_a.profile, node_c.profile))
                d_b_d = jukes_cantor(prof_dist(node_b.profile, node_d.profile))
                d_a_d = jukes_cantor(prof_dist(node_a.profile, node_d.profile))
                d_b_c = jukes_cantor(prof_dist(node_b.profile, node_c.profile))

                # now there are 2 possible cases in which we need to rearrange the nodes
                switch_flag = False
                if (d_a_c + d_b_d) < min((d_a_b + d_c_d), (d_a_d + d_b_c)):
                    # rearrange so (A, C) and (B, D) are together (new profiles?)
                    node.left.right = node_c
                    node.right.left = node_b
                    switch_flag = True

                    print("switch " + node_c.indexes + " " + node_b.indexes)

                elif (d_a_d + d_b_c) < min((d_a_b + d_c_d), (d_a_c + d_b_d)):
                    # rearrange so (A, D) and (B, C) are together(new profiles?)
                    switch_flag = True
                    node.left.right = node_d
                    node.right.left = node_b
                    print("switch " + node_d.indexes + " " + node_b.indexes)

                if switch_flag:
                    node.left.profile = profile_join(
                        node_a.profile, node_b.profile)
                    node.right.profile = profile_join(
                        node_c.profile, node_d.profile)
                    recalculate_profiles(node)

        # print(node.indexes)


def main():
    seqs = read_file('data\\test-small.aln')


def calc_d(i: Node, j: Node):
    """
    Calculates the d(i,j) equation from the paper.

        Parameters:
            i, j: the nodes to be compared.
        Returns:
            The distance between two nodes.
    """
    res = prof_dist(i.profile, j.profile) - i.up_distance - j.up_distance
    return res


def dist_to_all_other_nodes(i: int, nodes):
    """
    Calculates the distance to all nodes.

    Equivalent to total matrix, but only used for testing.
    """
    dist = 0
    for j in range(len(nodes)-1):
        if nodes[j].is_active and j != i:
            dist += calc_d(nodes[i], nodes[j])
    return dist


def branch_len(i, j, n, nodes):
    """
    Helper function for get_branch_lengths.

    Based on the algorithm from Wikipedia.
    """
    return 0.5 * calc_d(nodes[i], nodes[j]) + 1/(2*(n-2)) *\
        (dist_to_all_other_nodes(i, nodes) - dist_to_all_other_nodes(j, nodes))


def get_branch_lengths(i, j, active_nodes, nodes):
    """
    Calculates branch length based on the Wikipedia article
    on neighbor joining. Can be used to test the other branch
    length calculation function.
    """
    i_len = branch_len(i, j, active_nodes, nodes)
    i_len = round(i_len, 3)
    j_len = calc_d(nodes[i], nodes[j]) - i_len
    j_len = round(j_len, 3)

    if i_len < 0:
        j_len += -i_len
        i_len = 0
    elif j_len < 0:
        i_len += -j_len
        j_len = 0
    return i_len, j_len


def join_two_nodes(i: int, j: int, nodes: list) -> Node:
    """
    Joins two nodes into a new node.

        Parameters:
            i, j: the indices of the nodes to be joined
            nodes: the list of all nodes in the tree.

        Returns:
            The newly created node.
    """
    new_prof = profile_join(nodes[i].profile, nodes[j].profile)
    nodes[i].is_active = False
    nodes[j].is_active = False

    up_dist = prof_dist(nodes[i].profile, nodes[j].profile)/2

    # Make a new node.
    new_node = Node(profile=new_prof, up_distance=up_dist, is_active=True)
    new_node.left = nodes[i]    # we add the children of the node
    new_node.right = nodes[j]

    nodes[i].parent = new_node
    nodes[j].parent = new_node

    return new_node


def get_node_value(i: int, j: int, nodes: list, n: int, total_profile):
    """
    Updates each newly created node with the values of its parent nodes.
    This way, the last node in the tree will contain the whole tree.

        Returns:
            The name of the newly created node.
    """
    # TODO calc branch lengths correctly
    i_len, j_len = calc_branch_len_without_totprof(
        i, j, n, nodes)
    # print("branch length without total profile: ", i_len, j_len, " using nodes: ", i, j)
    #i_len, j_len = calc_branch_len(nodes[i], nodes[j], n, total_profile)
    # print("branch with total profile: ", i_len, j_len, " using nodes: ", i, j)

    # Building up the phylogenetic tree in Newick format.
    return "(" + nodes[i].value + ":" + str(i_len) +\
        "," + nodes[j].value + ":" + str(j_len) + ")"


def join_last_nodes(nodes, total_profile):
    last_nodes = []
    for i in range(len(nodes)):
        if nodes[i].is_active:
            last_nodes.append(i)
    last_node = join_two_nodes(last_nodes[0], last_nodes[1], nodes)
    last_node.value = get_node_value(
        last_nodes[0], last_nodes[1], nodes, 3, total_profile)

    return last_node


def create_phylogenetic_tree(nodes: list):
    """
    Loops through all active nodes and joins the nodes
    with the minimum distance. Stops when the number
    of active nodes is 3.

        Parameters:
            nodes: the initial leaf nodes created from the input sequences

        Returns:
            The phylogenetic tree in Newick format.

    This will complete the tree and make sure it is in correct format.
    """
    for active_nodes in range(len(nodes), 2, -1):
        total_profile = compute_total_profile(nodes, active_nodes)
        i, j, join_criterion = find_min_dist_nodes(
            nodes, total_profile, active_nodes)
        new_node = join_two_nodes(i, j, nodes)

        new_node.value = get_node_value(
            i, j, nodes, active_nodes, total_profile)
        new_node.index = len(nodes) - 1
        nodes.append(new_node)

    last_node = join_last_nodes(nodes, total_profile)
    print(last_node.value)
    return last_node


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

def dfs_search(node, val):
    stack = []
    stack.append(node)

    while len(stack) != 0:
        curr_node = stack.pop()
        if curr_node.value == val:
            return curr_node
        else:
            if curr_node.left != None:
                stack.append(curr_node.left)
            if curr_node.right != None:
                stack.append(curr_node.right)


def test_nearest_neighbor_interchange(root):
    node_7 = dfs_search(root, "7")
    node_1 = dfs_search(root, "1")
    print("found node: " + node_1.value)

    node_1.parent.right = node_7
    node_7.parent.left = node_1

    node_1.parent = node_7.parent
    node_7.parent = node_1.parent

    recalculate_profiles(node_1.parent)
    recalculate_profiles(node_7.parent)

    nearest_neighbor_interchanges(root)

    print("in test")

    # node_1.parent


def update_values(root, nodes):
    num_active_nodes = len(nodes)
    queue = []
    stack = []

    queue.append(root)

    while len(queue) > 0:
        node = queue.pop(0)
        if not node.left == None:
            queue.append(node.left)
            stack.append(node.left)
        if not node.right == None:
            queue.append(node.right)
            stack.append(node.right)

    while len(stack) > 0:
        node = stack.pop()
        print(node.value)
        if node.left != None and node.right != None:
            total_profile = compute_total_profile(nodes, num_active_nodes)
            new_value = get_node_value(
                node.left.index, node.right.index, nodes, num_active_nodes, total_profile)
            node.value = new_value
            # update num_active_nodes


def main():
    seqs = read_file('data\\test-small.aln')
    print(seqs)
    nodes = initialize_leaf_nodes(seqs)

    # PROOF that equation on the right bottom of page 3 in old paper equals comparison to 'total profile'.
    # equation is needed for computing r(i) and r(j)
    # i = make_profile(seqs[0])
    # print(find_total_profile(i, nodes))
    # print(prof_dist(i, compute_total_profile(seqs))*len(seqs))
    ###

    root = create_phylogenetic_tree(nodes)
    # update_values(root, nodes)
    print("before NNI")
    # test_nearest_neighbor_interchange(root)
    print("here")


if __name__ == '__main__':
    main()
