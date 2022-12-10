class Node:
    """
    Class representation of a node in the phylogenetic tree. 
    """
    profile: list
    up_distance: float 
    is_active: bool
    sequences: list
    indexes: str

    def __init__(self, profile=[], up_distance=0, is_active=False):
        self.profile = profile
        self.up_distance = up_distance
        self.is_active = is_active
        self.sequences = []
        self.indexes = ""
    
    def set_active_status(self, val:bool):
        self.is_active = val