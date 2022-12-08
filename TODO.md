# TO DO list

 - [ ] Create a function that can read the aln file and outputs a 2D list / matrix with the DNA sequence entries. Make sure the entries are parseable as a string.
 - [ ] Create a function that computes D(P_l(i,a), P_l(j, b))
 - [ ] Create a function that computes the profile distance P(i,j)
 - [ ] Create a function that computes “total profile” (P(i, T)) which is the average over all **active** nodes where T is the total profile
 - [ ] Create a function that computes u(i).
 - [ ] Create a function that computes d(i,j)
 - [ ] How to implement r(i) and r(j)? If anyone knows let us know. 
 - [ ] Implement d'(i,j) 

## Notes

The distance bewtween two profiles can be computed in O(La) time,instead of the naive O(La^2) time, by using the eigen decomposition of the alphabet’s distance matrix (see Methods)