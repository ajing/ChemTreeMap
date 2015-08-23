# A simple version for tree building

# Architecture
## Similarity Calculation
input: smile
output: distance matrix

## Preparing input for RapidNJ
input: distance matrix
output: RapidNJ input ([PHYLIP format](http://www.mothur.org/wiki/Phylip-formatted_distance_matrix))

## Run RapidNJ
input: filename for RapidNJ input
