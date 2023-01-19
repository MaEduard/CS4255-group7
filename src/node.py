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
    index: int
    value: str #name of the node itself --> e.g. node 1 merged with node 2 is called '12'
    main_sequence = ""

    def __init__(self, profile=[], up_distance=0, is_active=False, main_sequence="", top_hits = []):
        self.profile = profile
        self.up_distance = up_distance
        self.is_active = is_active
        self.top_hits = []
        self.indexes = None
        self.main_sequence = ""
        self.index = -1
        self.value = ""

    def set_active_status(self, val:bool):
        self.is_active = val