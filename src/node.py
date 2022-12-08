class Node:
    def __init__(self, sequence):
        self.sequences = [sequence]
        self.profile = self.createProfile()
        self.up_distance = 0

    def createProfile(self):
        """Computes the profile based on the sequences that are in the node, WITHOUT using pseudocounts

        Returns:
            float[][]: 4 x k profile matrix
        """
        profile_temp = [[0 for i in range(len(self.sequences[0]))] for j in range(4)]
        for column in range(len(self.sequences[0])):
            count = [0,0,0,0]
            # go over every motif per column
            for row in range(len(self.sequences)):
                if self.sequences[row][column] == 'A':
                    count[0] += 1
                elif self.sequences[row][column] == 'C':
                    count[1] += 1
                elif self.sequences[row][column] == 'G':
                    count[2] += 1
                else:
                    count[3] +=1
            for profile_row in range(4):
                profile_temp[profile_row][column] = count[profile_row]/(len(self.sequences))
        return profile_temp

    def __str__(self):
        return f"{self.sequences}({self.profile})"
