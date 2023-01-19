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


def node_dist(i: Node, j: Node):
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


def r_equation(node_i, total_profile, n, up_sum):
    ans = n*prof_dist(node_i.profile, total_profile) - prof_dist(node_i,
                                                                 node_i) - (n-1)*node_i.up_distance + node_i.up_distance - up_sum
    return ans/(n-2)


def help_function(i, nodes, n):
    dist = 0
    for node in nodes:
        if node.indexes != i.indexes:
            if node.is_active:
                dist += prof_dist(i.profile, node.profile) - \
                    i.up_distance - node.up_distance

    return dist/(n-2)


def find_min_dist_nodes(nodes, total_profile, n):
    """
    Finds the minimum distance between all possible joins, and returns the index of both nodes. 

    Using exhaustive search. 

    Parameters: 
        n = number of active nodes

    TODO include the r_equation to compute r(i) and r(j) to eventually compute d'(i,j) = d(i,j) - r(i) - r(j)
    """
    min_dist = float('inf')
    for i in range(len(nodes)):
        if nodes[i].is_active:  # Skip nodes that are inactive.
            for j in range(i+1, len(nodes)):
                if nodes[j].is_active:  # Skip nodes that are inactive.
                    dist_i_j = node_dist(nodes[i], nodes[j])
                    r_i = out_dist(nodes[i], total_profile, n)
                    r_j = out_dist(nodes[j], total_profile, n)
                    dist_prime_i_j = dist_i_j - r_i - r_j

                    if dist_prime_i_j < min_dist:
                        min_dist = dist_prime_i_j
                        ind1 = i
                        ind2 = j

    print(out_dist(nodes[ind1], total_profile, n))
    print(help_function(nodes[ind1], nodes, n))

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


def calc_up_dist(nodes):
    dist = 0
    for node in nodes:
        if node.is_active:
            dist += node.up_distance
    return dist


def out_dist(node, total_profile, n):
    """
      return (n*(prof_dist(node.profile, total_profile)) + tot_up_dist - node.up_distance) / (n-2)
    Calculates the out distance from one node. 
    """
    return (n*(prof_dist(node.profile, total_profile)) - prof_dist(node.profile, node.profile) - (n-1)*node.up_distance + node.up_distance - tot_up_dist) / (n-2)


def calc_branch_len(i: Node, j: Node, n, total_profile):
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


def set_last_join(nodes):
    active_joins = 3


def find_joins(nodes: list, seqs):
    """
    Loops through all active nodes and joins the nodes with the minimum distance. 
    Stops when the number of active nodes is just one. 

    """
    branch_lengths = []
    active_nodes = len(nodes)
    while active_nodes > 2:
        total_profile = compute_total_profile(nodes, active_nodes)

        i, j, mindist = find_min_dist_nodes(nodes, total_profile, active_nodes)

        print(nodes[i].indexes, "<-->", nodes[j].indexes)
        print("----------------")

        prof_i = nodes[i].profile
        prof_j = nodes[j].profile

        nodes[i].is_active = False
        nodes[j].is_active = False

        k = profile_join(prof_i, prof_j)

        up_dist = prof_dist(nodes[i].profile, nodes[j].profile)/2
        new_node = Node(profile=k, up_distance=up_dist, is_active=True)
        new_node.left = nodes[i]    # we add the children of the node
        new_node.right = nodes[j]
        nodes.append(new_node)

        # Set the new node to be the parent of the two joined nodes.
        nodes[i].parent = new_node
        nodes[j].parent = new_node

        new_node.indexes = "(" + nodes[i].indexes + \
            "," + nodes[j].indexes + ")"

        i_len, j_len = calc_branch_len(
            nodes[i], nodes[j], active_nodes, total_profile)  # TODO Remove this?
        branch_lengths.append(str(i) + " " + str(i_len))
        branch_lengths.append(str(j) + " " + str(j_len))
        active_nodes -= 1
        # break

    last_node1 = None
    last_node2 = None

    for node in nodes:
        if node.is_active:
            node.is_active = False
            if last_node1 is None:
                last_node1 = node
            else:
                last_node2 = node

    print(last_node1.indexes, "<-->", last_node2.indexes)
    print("----------------")

    root_profile = profile_join(last_node1.profile, last_node2.profile)
    root_up_dist = prof_dist(last_node1.profile, last_node2.profile)/2
    root = Node(profile=root_profile, up_distance=root_up_dist, is_active=True)
    root.left = last_node1
    root.right = last_node2

    return branch_lengths, root


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

    nodes = initialize_leaf_nodes(seqs)

    # PROOF that equation on the right bottom of page 3 in old paper equals comparison to 'total profile'.
    # equation is needed for computing r(i) and r(j)
    # i = make_profile(seqs[0])
    # print(find_total_profile(i, nodes))
    # print(prof_dist(i, compute_total_profile(seqs))*len(seqs))
    ###

    branch_length, root = find_joins(nodes, seqs)
    print(branch_length)
    nearest_neighbor_interchanges(root)


if __name__ == '__main__':
    main()
