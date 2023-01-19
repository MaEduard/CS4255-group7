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
    main_sequence = ""

    def __init__(self, profile=[], up_distance=0, is_active=False, main_sequence="", top_hits = []):
        self.profile = profile
        self.up_distance = up_distance
        self.is_active = is_active
        self.top_hits = []
        self.indexes = None
        self.main_sequence = ""
    
    def set_active_status(self, val:bool):
        self.is_active = val