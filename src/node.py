class Node:
    """
    Class representation of a node in the phylogenetic tree.
    """
    profile: list
    up_distance: float
    is_active: bool
    sequences: list
    top_hits = list
    indexes: str
    index: int # index value pointing to the node in the nodes list
    value: str #name of the node itself --> e.g. node 1 merged with node 2 is called '12'
    main_sequence = ""
    value: str
    index: int

    def __init__(self, profile=[], up_distance=0, is_active=False, main_sequence="", top_hits = [], index=0):
        self.profile = profile
        self.up_distance = up_distance
        self.is_active = is_active
        self.sequences = []
        self.indexes = ""
        self.left = None
        self.right = None
        self.parent = None
        self.top_hits = []
        self.main_sequence = ""
        self.value = ""
        self.index = index

    def set_active_status(self, val: bool):
        self.is_active = val

    def has_two_children(self):
        return (self.left != None and self.right != None)

    def set_active_status(self, val: bool):
        self.is_active = val
