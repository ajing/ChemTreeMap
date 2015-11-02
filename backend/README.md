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
output: newick tree format

## Newick to Dot
input: newick file
output: dot file

## run graphviz on dot
input: dot filename
output: dot filename

## convert dot file to json with combination of other info
input: dot filename, other info
output: json
