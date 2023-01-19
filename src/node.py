class Node:
    """
    Class representation of a node in the phylogenetic tree.
    """
    profile: list
    up_distance: float
    is_active: bool
    sequences: list
    top_hits = list
    main_sequence = ""
    value: str
    index: int

    def __init__(self, profile=[], up_distance=0, is_active=False, main_sequence=""):
        self.profile = profile
        self.up_distance = up_distance
        self.is_active = is_active
        self.top_hits = []
        self.main_sequence = ""
        self.value = ""
        self.index = None
    
    def set_active_status(self, val: bool):
        self.is_active = val
