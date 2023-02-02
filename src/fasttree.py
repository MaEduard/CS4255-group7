import math
from node import Node
import random


global tot_up_dist
global m


def read_file(file_name):
    """Reads the input file and outputs the strings.

    Args:
        file_name (str): string formatted file name location

    Returns:
        str[]: list of sequence strings
    """
    f = open(file_name, "r")
    sequences = []
    real_line = False
    for line in f.readlines():
        if real_line:
            sequences.append(line[:len(line)-1])  # to remove the \n
        real_line = not real_line

    return sequences


def letter_to_index(letter):
    """Converts a char letter to integer index.

    Args:
        letter (char): input letter (A/C/G/T)

    Returns:
        int: index
    """
    if letter == "A":
        return 0
    elif letter == "C":
        return 1
    elif letter == "G":
        return 2
    elif letter == "T":
        return 3


def jukes_cantor(d_u):
    """Computes the Jukes Cantor distance for the NNI heuristic.

    Args:
        d_u (float): uncorrected distance

    Returns:
        float: Jukes Cantor distance
    """
    dist = -(3/4)*math.log(1-(4/3)*d_u)

    return dist


def make_profile(seq):
    """Produces frequency profile based on given string sequence.

    Args:
        seq (str): DNA sequence

    Returns:
        float[][]: 2D (4 x len(seq)) frequency matrix
    """
    prof = make_empty_profile(len(seq))
    for i, let in enumerate(seq):
        prof[letter_to_index(let)][i] = 1
    return prof


def make_empty_profile(size):
    """Creates empty frequency profile.

    Args:
        size (int): number of columns for the frequency profile

    Returns:
        int[][]: 4 x size empty (i.e. 0) frequency profile
    """
    profile = [[], [], [], []]
    for row in profile:
        for i in range(size):
            row.append(0)

    return profile


def profile_join(a, b):
    """Joins two profile together and returns the resulting profile.

    Args:
        a (float[][]): profile of first node to be merged
        b (float[][]): profile of second node to be merged

    Returns:
        float[][]: newly merged profile
    """
    new_prof = make_empty_profile(len(a[0]))
    for row in range(4):  # Only 4 nucleotides.
        for col in range(len(a[0])):  # Assume same size for a and j
            new_prof[row][col] = (a[row][col] + b[row][col])/2
    return new_prof


def prof_dist(prof1, prof2):
    """Computes the distance between profiles.

    Args:
        prof1 (float[][]): _description_
        prof2 (float[][]): _description_

    Returns:
        float: distance between two profiles
    """
    dist = 0
    k = len(prof1[0])
    for a in range(k):
        for x in range(4):  # Checking for A/C/G/T nucleotides in prof1.
            for y in range(4):  # Checking for A/C/G/T nucleotides in prof2.
                if x != y:
                    dist += prof1[x][a] * prof2[y][a]
    return dist/k


def profile_join(prof1, prof2):
    """Joins two profile matrices based by taking the average of the two input profiles.

    Args:
        prof1 (_type_): _description_
        prof2 (_type_): _description_

    Returns:
        float[][]: joined profile matrix
    """
    new_prof = make_empty_profile(len(prof1[0]))
    for row in range(4):  # Only 4 nucleotides.
        # Assume same size for prof1 and prof2
        for col in range(len(prof1[0])):
            new_prof[row][col] = (prof1[row][col] + prof2[row][col])/2
    return new_prof


def out_dist(node, total_profile, n):
    """Computes the out distance {'r(i)' in the paper}  based on the total profile.

    Args:
        node (Node): node to calculate the out distance of
        total_profile (float[][]): total profile so far.
        n (int): number of active nodes

    Returns:
        float: out distance
    """
    return (n*(prof_dist(node.profile, total_profile)) - prof_dist(node.profile, node.profile) - tot_up_dist - (n-1)*node.up_distance) / (n-2)


def join_criterion(a, b, n, total_profile):
    """Calculates the Neighbor Joining criterion value that we use to find the best nodes to join (the nodes for which this criterion is 
    minimized using top hits heuristic).

    Args:
        a (Node): node a
        b (Node): node b
        n (int): number of active nodes
        total_profile (float[][]): _description_

    Returns:
        float: join criterion value
    """
    dist_i_j = calc_d(a, b)
    r_j = out_dist(a, total_profile, n)
    r_i = out_dist(b, total_profile, n)
    return dist_i_j - r_i - r_j


def find_nodes_to_be_joined(nodes, total_profile, n):
    """Finds the best nodes to be joined based on their top hits list and the Neighbor Joining Criterion. Updates a node's top hits list in a
    a lazy way (on the fly) when it encounters a node that has already been joined.

    Args:
        nodes (Node[]): list of all nodes created so far
        total_profile (float[][]): profile matrix based on all nodes.
        n (int): number of active nodes.

    Returns:
        (Node), (Node): the best nodes to be joined according to the joining criteria.
    """

    min_dist = float('inf')
    node1 = None
    node2 = None
    for i in range(len(nodes)):
        if nodes[i].is_active:  # Skip nodes that are inactive.
            current_dist = nodes[i].top_hits[0][0]
            top_node_i = nodes[i].top_hits[0][1]
            top_node = nodes[top_node_i]
            if top_node.is_active == False:
                while top_node.is_active == False:
                    top_node = top_node.parent
                dist = join_criterion(nodes[i], top_node, n, total_profile)
                nodes[i].top_hits[0] = (dist, top_node.index)
                current_dist = dist

            if current_dist < min_dist:
                min_dist = current_dist
                node1 = i
                node2 = nodes[i].top_hits[0][1]

    return node1, node2


def merge_top_hits(top_hits1, top_hits2):
    """Merges the top hits lists of the nodes that were found to be best to join. Removes any duplicate top hits after merging.
    Notice that our top hits list stores both the neighbor joining criterion value and the index of the node that is a top hit for the 
    node to which the top hits list belongs.

    Args:
        top_hits1 ((dist, Node)[]): top hits list of fist node to be joined
        top_hits2 ((dist, Node[])): top hits list of second node to be joined

    Returns:
        ((dist, Node)[]): merged list of top hits.
    """
    top_hits_merged = top_hits1 + top_hits2
    top_hits_merged = list(set(top_hits_merged))

    return top_hits_merged


def filter_top_hits(merged_node, top_hits, active_nodes, total_profile, nodes):
    """Filters the merged top hits list (created by the top hits lists of the nodes that were joined) to represent the best 
    top hits for this newly created node called 'merge_node').

    Args:
        merged_node (Node): Newly created node after a join
        top_hits ((dist, Node)[]): list of top hits created after merging top hits lists of both nodes
        active_nodes (int): number of active nodes
        total_profile (float[][]): total profile matrix based on all nodes that are active
        nodes (Node[]): list of all nodes

    Returns:
        ((dist, Node)[]): newly created top hits list (of length m) for merged node 
    """
    new_top_hits = []
    for node_tuple in top_hits:
        node_index = node_tuple[1]
        curr = nodes[node_index]
        # If the best hit is inactive, iterate until active ancestor is found
        while curr.is_active == False:
            curr = curr.parent
            node_index = curr.index

        if node_index == merged_node.index:  # If pointing to self
            continue

        dist_i_j = calc_d(merged_node, nodes[node_index])
        r_j = out_dist(merged_node, total_profile, active_nodes)
        r_i = out_dist(nodes[node_index], total_profile, active_nodes)
        dist_prime_i_j = dist_i_j - r_i - r_j
        new_top_hits.append((dist_prime_i_j, node_index))

    new_top_hits.sort()
    return new_top_hits[0:m]


def update_top_hits(merged_node, nodes, n, total_profile):
    """Refreshes the top hits list of the merged node when the size of the list has shrunk below 0.8*m

    Args:
        merged_node (Node): node that has just been merged
        nodes (Node[]): list of all nodes
        n (int): number of active nodes
        total_profile (float[][]): total profile derived from all current active nodes
    """
    new_top_hits_list = []

    for node in nodes:
        if node.is_active and node.index != merged_node.index:
            dist = join_criterion(merged_node, node, n, total_profile)
            new_top_hits_list.append((dist, node.index))

    new_top_hits_list.sort()
    merged_node.top_hits = new_top_hits_list[0:m]

    # Updates neighbors top hits lists.
    for (_, i) in new_top_hits_list[0:m]:
        neighbor_node = nodes[i]
        new_neighbor_top_hits = []
        for (_, j) in new_top_hits_list[0:2*m]:
            other_node = nodes[j]
            if other_node.is_active and neighbor_node.index != other_node.index:
                dist = join_criterion(
                    neighbor_node, other_node, n, total_profile)
                new_neighbor_top_hits.append((dist, other_node.index))
        new_neighbor_top_hits = list(
            set(new_neighbor_top_hits + neighbor_node.top_hits))
        new_neighbor_top_hits.sort()
        neighbor_node.top_hits = new_neighbor_top_hits[0:m]


def join_nodes(node1, node2, active_nodes, nodes, total_profile):
    """Joins node1 and node2 by setting them to inactive and creates the merged node with a new profile based on the 
    profiles of node1 and node2 and the newly calculated up distance. Takes into account the old top hits lists update and
    updates it appropriately.

    Args:
        node1 (Node): first node to join
        node2 (Node): second node to join
        active_nodes (int): number of active nodes currently
        nodes (Node[]): list of all nodes
        total_profile (flaot[][]): profile based on all active nodes

    Returns:
        merged_node (Node): the newly created merged node
    """
    nodes[node1].is_active = False
    nodes[node2].is_active = False

    new_profile = profile_join(nodes[node1].profile, nodes[node2].profile)
    new_updistance = prof_dist(nodes[node1].profile, nodes[node2].profile)/2
    merged_node = Node(
        profile=new_profile, up_distance=new_updistance, is_active=True, index=len(nodes))

    nodes[node1].parent = merged_node
    nodes[node2].parent = merged_node
    merged_node.right = nodes[node1]
    merged_node.left = nodes[node2]

    top_hits_merged = merge_top_hits(
        nodes[node1].top_hits, nodes[node2].top_hits)
    filtered_top_hits = filter_top_hits(
        merged_node, top_hits_merged, active_nodes, total_profile, nodes)
    merged_node.top_hits = filtered_top_hits

    if len(filtered_top_hits) < ((0.8*m)):
        update_top_hits(merged_node, nodes, active_nodes, total_profile)

    merged_node.value = get_node_value(
        node1, node2, nodes, active_nodes, total_profile)

    return merged_node


def initialize_leaf_nodes(seqs):
    """Creates a leaf node for all sequences in the input at the start of the program.

    Args:
        seqs (str[]): the list of sequences to be analyzed.

    Returns:
        nodes (Node[]): list of nodes that represent leaves in the tree.
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
    """Sets the global varibale tot_up_dist

    Args:
        val (float): value to update the tot_up_dist
    """
    global tot_up_dist
    tot_up_dist = val


def compute_total_profile(nodes, active_nodes):
    """Computes the profile matrix of all active nodes.

    Also sets the global variable for the total up distance.

    Args:
        nodes (Node[]): all the nodes of the tree
        active_nodes (int): the number of active nodes in the tree.

    Returns:
        profile (float[]): profile matrix of size 4 x k
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
    """Computes the out distance {'r(i)' in the paper}  based on the total profile.

    Args:
        node (Node): node to calculate the out distance of
        total_profile (float[][]): total profile so far.
        n (int): number of active nodes

    Returns:
        float: out distance
    """
    return (n*(prof_dist(node.profile, total_profile)) -
            prof_dist(node.profile, node.profile) -
            (n-1)*node.up_distance + node.up_distance - tot_up_dist) / (n-2)


def calc_branch_len(i, j, n, total_profile):
    """Calculates the branch length from daughter nodes i and j.

    Args:
        i (Node): node i
        j (Node): node j
        n (int): number of active nodes
        total_profile (float[][]): profile based on all active nodes

    Returns:
        (float, float): length of the branch between node i and j
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


def recalculate_profiles(node):
    """After a nearest neighbor interchange, this method recalculates the profiles
    of the nodes influenced by the change by propagating it up to the root of the tree.
    This is done recursively.

    Args:
        node (Node): node to be propagated
    """
    if node.parent is None:
        return
    new_profile = profile_join(node.left.profile, node.right.profile)
    node.profile = new_profile

    recalculate_profiles(node.parent)


def nearest_neighbor_interchanges(root):
    """Calculates log-corrected distances and based on that switches nodes in the topology.
    Every time there is a switch the profiles of the nodes are also recalculated.

    Args:
        root (Node): the node at the root of the create tree
    """
    queue = []
    stack = []

    queue.append(root)

    while len(queue) > 0:
        node = queue.pop(0)
        if node.left is not None:
            queue.append(node.left)
            stack.append(node.left)
        if node.right is not None:
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

                elif (d_a_d + d_b_c) < min((d_a_b + d_c_d), (d_a_c + d_b_d)):
                    # rearrange so (A, D) and (B, C) are together(new profiles?)
                    switch_flag = True
                    node.left.right = node_d
                    node.right.left = node_b

                if switch_flag:
                    node.left.profile = profile_join(
                        node_a.profile, node_b.profile)
                    node.right.profile = profile_join(
                        node_c.profile, node_d.profile)
                    recalculate_profiles(node)


def calc_d(i: Node, j: Node):
    """Calculates the distance d(i,j) between nodes i and j.

    Args:
        i (Node): node i
        j (Node): node j

    Returns:
        res (float): node distance
    """
    res = prof_dist(i.profile, j.profile) - i.up_distance - j.up_distance
    return res


def dist_to_all_other_nodes(i, nodes):
    """Calculates the distance to all nodes. Equivalent to total profile, but only used for testing.

    Args:
        i (int): node index
        nodes (Node[]): list of all nodes

    Returns:
        dist (int): distance between node i and all other nodes.
    """
    dist = 0
    for j in range(len(nodes)-1):
        if nodes[j].is_active and j != i:
            dist += calc_d(nodes[i], nodes[j])
    return dist


def branch_len(i, j, n, nodes):
    """Helper function for get_branch_lengths. Based on the algorithm from Wikipedia
    (https://en.wikipedia.org/wiki/Neighbor_joining#Distance_from_the_pair_members_to_the_new_node).
    Calculates the branch length of to newly created node.

    Args:
        i (int): node index for node i
        j (int): node index for node j
        n (int): number of active nodes currently
        nodes (Node[]): list all nodes

    Returns:
        float: distance
    """
    return 0.5 * calc_d(nodes[i], nodes[j]) + 1/(2*(n-2)) *\
        (dist_to_all_other_nodes(i, nodes) - dist_to_all_other_nodes(j, nodes))


def get_branch_lengths(i, j, active_nodes, nodes):
    """Calculates branch length based on the Wikipedia article
    on neighbor joining. Can be used to test the other branch
    length calculation function.

    Args:
        i (_type_): _description_
        j (_type_): _description_
        active_nodes (_type_): _description_
        nodes (_type_): _description_

    Returns:
        _type_: _description_
    """
    i_len = branch_len(i, j, active_nodes, nodes)
    i_len = round(i_len, 3)
    j_len = calc_d(nodes[i], nodes[j]) - i_len
    j_len = round(j_len, 3)

    # If branch lengths are negative, make them zero and subtract the value
    # from the other branch.
    if i_len < 0:
        j_len += -i_len
        i_len = 0
    elif j_len < 0:
        i_len += -j_len
        j_len = 0
    return i_len, j_len


def join_last_nodes(i, j, nodes):
    """Joins node i and j and outputs the newly created node.

    Args:
        i (int): index of first node
        j (int): index of second node
        nodes (Node[]): list of all nodes

    Returns:
        Node: newly created node based on nodes[i] and nodes[j]
    """
    # Create a new profile based on the daughter nodes.
    new_prof = profile_join(nodes[i].profile, nodes[j].profile)
    nodes[i].is_active = False
    nodes[j].is_active = False

    up_dist = prof_dist(nodes[i].profile, nodes[j].profile)/2

    # Make a new node.
    new_node = Node(profile=new_prof, up_distance=up_dist, is_active=True)
    new_node.right = nodes[i]    # we add the children of the node
    new_node.left = nodes[j]

    nodes[i].parent = new_node
    nodes[j].parent = new_node

    return new_node


def get_node_value(i, j, nodes, n, total_profile):
    """Updates each newly created node with the values of its parent nodes.
    This way, the last node in the tree will contain the whole tree.

    Args:
        i (int): index of first daughter node
        j (int): index of second daughter node
        nodes (Node[]): list of all nodes
        n (int): number of active nodes currently
        total_profile (float[][]): profile based on all active nodes

    Returns:
        str: The value of the newly created node.
    """
    i_len, j_len = calc_branch_len(
        nodes[i], nodes[j], n, total_profile)

    # Building up the phylogenetic tree in Newick format.
    return "(" + nodes[i].value + ":" + str(i_len) +\
        "," + nodes[j].value + ":" + str(j_len) + ")"


def find_last_join(nodes, total_profile):
    """Joins the last two active nodes, and returns the root of the tree.

    Parameters:
        nodes (Node[]): the list containing the nodes
        total_profile (float[][]): profile based on all active nodes

    Returns:
        last_node (Node): root node of the tree
    """
    last_nodes = []
    for i in range(len(nodes)):
        if nodes[i].is_active:
            last_nodes.append(i)

    last_node = join_last_nodes(last_nodes[0], last_nodes[1], nodes)
    last_node.value = get_node_value(
        last_nodes[0], last_nodes[1], nodes, 3, total_profile)

    return last_node


def create_phylogenetic_tree(nodes):
    """Loops through all active nodes and joins the nodes
    with the minimum distance, based on top hits heuristics.

    Stops when only two remaining nodes are active.

    Parameters:
        nodes (Node[]): the initial leaf nodes created from the input sequences

    Returns:
        last_node (Node): The root of the tree that contains the full tree information in last_node.value

    """
    create_top_hits(nodes, len(nodes))
    initial_nodes = len(nodes)

    # Loop until only two of the remaining leaves are left.
    for active_nodes in range(initial_nodes, 2, -1):
        total_profile = compute_total_profile(nodes, active_nodes)
        node1, node2 = find_nodes_to_be_joined(
            nodes, total_profile, active_nodes)
        new_node = join_nodes(node1, node2, active_nodes, nodes, total_profile)
        nodes.append(new_node)

    # 2 active nodes left. Merge these.
    last_node = find_last_join(nodes, total_profile)
    return last_node


def create_top_hits(nodes, n):
    """
    Method that will calculate top hits list of all initial leaf nodes. Starts of with a seed node and computes distance 
    between side node and all other nodes. Sorts the list and picks the top m hits as the top hits list for the seed node. 
    Furthermore, for these m hits, finds their top m hits from the top 2*m hits of the seed. Exhaustively repeats this full process until all 
    nodes have their top hits list created.

    Args:
        nodes (Node[]): list all nodes
        n (int): number of active nodes
    """
    m = int(math.sqrt(len(nodes)))  # maybe use math.ceil?
    unused_nodes = [i for i in range(len(nodes))]  # list of index integers
    total_profile = compute_total_profile(nodes, n)

    while (len(unused_nodes) > 0):
        # pick random seed node. Note that randint picks random int from inclusive range, hence do len(seed_nodes) - 1
        seed_idx = random.sample(unused_nodes, 1)[0]
        seed = nodes[seed_idx]
        top_hits = []

        for node in nodes:
            if node.index != seed.index:
                dist_prime = join_criterion(seed, node, n, total_profile)

                top_hits.append((dist_prime, node.index))

        # the distances for some top hits are apparently exactly the same and thus it starts to compare nodes, which is impossible.
        top_hits.sort()

        nodes[seed_idx].top_hits = top_hits[0:m]
        unused_nodes.remove(seed_idx)

        for neighbor in top_hits[0:m]:
            neighbor_node = nodes[neighbor[1]]
            if len(neighbor_node.top_hits) != 0:
                continue
            top_hits_neighbor = []
            for node_tuple in top_hits[0:2*m]:
                node = nodes[node_tuple[1]]
                if node.index != neighbor_node.index:
                    dist = join_criterion(
                        node, neighbor_node, n, total_profile)
                    top_hits_neighbor.append((dist, node.index))

            top_hits_neighbor.sort()
            neighbor_node.top_hits = top_hits_neighbor[0:m]
            if neighbor_node in unused_nodes:
                unused_nodes.remove(neighbor_node)


def dfs_search(node, val):
    """Performs a depth-first search of the tree and returns the node
    that corresponds to the value passed as a parameter.

    Args:
        node (Node): the starting node of the search (in our case - root)
        val (String): the value the search is looking for in the nodes
    """
    stack = []
    stack.append(node)

    while len(stack) != 0:
        curr_node = stack.pop()
        if curr_node.value == val:
            return curr_node
        else:
            if curr_node.left is not None:
                stack.append(curr_node.left)
            if curr_node.right is not None:
                stack.append(curr_node.right)


def test_nearest_neighbor_interchange(root, nodes):
    """This test shows the efficacy of NNI. It switches two nodes in the topology
    which our phylogenic tree algorithm came up with and performs NNI on that
    "false" topology. The result after the NNI is that the topology is brought
    back to the original one before the manual interchange. This proves that NNI
    at least in this particular case executes the right interchange of nodes.
    By running the test the topologies are printed in the console in Newick
    format and easy to understand.

    Args:
        root (Node): the node at the root of the create tree
        nodes (Node[]): list of all nodes
    """
    node_7 = dfs_search(root, "7")
    node_1 = dfs_search(root, "1")

    node_1.parent.left = node_7
    node_7.parent.right = node_1

    node_1.parent = node_7.parent
    node_7.parent = node_1.parent

    recalculate_profiles(node_1.parent)
    recalculate_profiles(node_7.parent)

    print("--------Test---------")

    update_values(root, nodes)
    print("After node switching:")
    print(root.value)

    nearest_neighbor_interchanges(root)
    update_values(root, nodes)

    print("After NNI:")
    print(root.value)
    print("---------End Test-------")

    # node_1.parent


def update_values(root, nodes):
    """In case any switching of nodes has happened (in NNI or in testing),
    the values of the nodes which are the representation of the joins and lengths
    still need to be updated because the switching in NNI happens only at the
    references in the nodes. This method handles the values update.

    Args:
        root (Node): the node at the root of the create tree
        nodes (Node[]): list of all nodes
    """
    # Here I caclculate the number of nodes in the final tree because we don't keep track of that number anywhere else
    num_active_nodes = 1
    queue = []
    stack = []

    queue.append(root)

    # traverse the nodes in an order starting from the leaves
    while len(queue) > 0:
        node = queue.pop(0)
        if node.left is not None:
            queue.append(node.left)
            stack.append(node.left)
            num_active_nodes += 1
        if node.right is not None:
            queue.append(node.right)
            stack.append(node.right)
            num_active_nodes += 1

    stack.insert(0, root)
    while len(stack) > 0:
        node = stack.pop()
        if node.left is not None and node.right is not None:
            total_profile = compute_total_profile(nodes, num_active_nodes)
            new_value = get_node_value(
                node.left.index, node.right.index, nodes, num_active_nodes, total_profile)
            node.value = new_value
            num_active_nodes -= 1


def main():
    """
    Main method to start the algorithm. Loads the data and starts FastTree.
    """
    try:
        seqs = read_file('data/test-small.aln')
    except FileNotFoundError:
        seqs = read_file('../data/test-small.aln')

    nodes = initialize_leaf_nodes(seqs)

    global m  # Used globally for top hits heuristics.
    m = round(math.sqrt(len(nodes)))

    # Create initial topology.
    root = create_phylogenetic_tree(nodes)

    nni_iterations = math.ceil(math.log2(len(nodes))) + 1
    for _ in range(nni_iterations):
        nearest_neighbor_interchanges(root)

    update_values(root, nodes)
    test_nearest_neighbor_interchange(root, nodes)
    print(root.value)


if __name__ == '__main__':
    main()
