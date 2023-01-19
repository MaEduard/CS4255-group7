class Node:
    """
    Class representation of a node in the phylogenetic tree. 
    """
    # profile: list
    # up_distance: float
    # is_active: bool
    # sequences: list
    # indexes: str
    # left: Node

    def __init__(self, profile=[], up_distance=0, is_active=False):
        self.profile = profile
        self.up_distance = up_distance
        self.is_active = is_active
        self.sequences = []
        self.indexes = ""
        self.left = None
        self.right = None
        self.parent = None

    def set_active_status(self, val: bool):
        self.is_active = val

    def has_two_children(self):
        return (self.left != None and self.right != None)
