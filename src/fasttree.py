from node import Node

tot_up_dist = 0


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

                    dist_prime_i_j = dist_i_j - r_i - r_j

                    if dist_prime_i_j < min_dist:
                        #print(dist_i_j, r_i, r_j, "nodes used: ", i, j)
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
        nodes.append(node)

    return nodes


def print_profile(a):
    """
    Used for testing.
    """
    for row in a:
        print(row[0:4])


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


def find_total_profile(prof_i, nodes):
    """
    Computes the total distance of node i to other nodes in the old fashion
    way. We should not do this but actually compute this value by comparing
    node i to the total profile of all active nodes.

    * I created this function to double check whether the formula on the right
    bottom of page 3 of old paper actually holds.
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
    return (n*(prof_dist(node.profile, total_profile)) -
            prof_dist(node.profile, node.profile) -
            (n-1)*node.up_distance + node.up_distance - tot_up_dist) / (n-2)


def calc_branch_len(i: Node, j: Node, n, total_profile) -> tuple[float, float]:
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
            dist += calc_d(nodes[i], nodes[j]) - \
                nodes[i].up_distance - nodes[j].up_distance
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
    up_dist = calc_d(nodes[i], nodes[j])/2

    # Make a new node.
    new_node = Node(profile=new_prof, up_distance=up_dist, is_active=True)

    return new_node


def get_node_value(i: int, j: int, nodes: list, n: int, total_profile):
    """
    Updates each newly created node with the values of its parent nodes.
    This way, the last node in the tree will contain the whole tree.

        Returns:
            The name of the newly created node.
    """
    i_len, j_len = calc_branch_len(
        nodes[i], nodes[j], n, total_profile)
    # i_len, j_len = get_branch_lengths(i, j, n, nodes)

    # Building up the phylogenetic tree in Newick format.
    return "(" + nodes[i].value + ":" + str(i_len) +\
        "," + nodes[j].value + ":" + str(j_len) + ")"


def create_phylogenetic_tree(nodes: list):
    """
    Loops through all active nodes and joins the nodes
    with the minimum distance. Stops when the number
    of active nodes is 3.

        Parameters:
            nodes: the initial leaf nodes created from the input sequences

        Returns:
            The phylogenetic tree in Newick format.

    TODO join all the active nodes. Perhaps in a different function.
    This will complete the tree and make sure it is in correct format.
    """
    for active_nodes in range(len(nodes), 2, -1):
        total_profile = compute_total_profile(nodes, len(nodes))
        i, j, mindist = find_min_dist_nodes(nodes, total_profile, active_nodes)

        new_node = join_two_nodes(i, j, nodes)
        new_node.value = get_node_value(
            i, j, nodes, active_nodes, total_profile)
        nodes.append(new_node)


    print(new_node.value)


def main():
    seqs = read_file('../data/test-small.aln')

    nodes = initialize_leaf_nodes(seqs)

    create_phylogenetic_tree(nodes)


if __name__ == '__main__':
    main()
