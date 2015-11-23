# A simple version for tree building

# Architecture
## Similarity Calculation
file:     DistanceMeasure.py
function: GenerateDistanceFile()
input:    smile
output:   distance matrix

#o# Preparing input for RapidNJ
file:     DistanceMeasure.py
function: GenerateDistanceFile()
input:    distance matrix
output:   RapidNJ input ([PHYLIP format](http://www.mothur.org/wiki/Phylip-formatted_distance_matrix))

## Run RapidNJ
file:     RunRapidNJ.py
function: RunRapidNJ()
input:    filename for RapidNJ input
output:   newick tree format

## Newick to Dot
file:     RunGraphViz.py
function: WriteDotFile()
input:    newick file
output:   dot file

## run graphviz on dot
file:     RunGraphViz.py
function: SFDPonDot()
input:    dot filename
output:   dot filename

## convert dot file to json with combination of other info
file:     RunGraphViz.py
function: Dot2Dict()
input:    dot filename, other info
output:   json

## For adding other physical properties
file:     DistanceMeasure.py
function: AddLigEff, AddSLogP
input:    a json for server input, ligand profile dictionary
output:   a new json

## If you have additional columns in your input file, you can change INTEREST_COLUMN in Model.py
