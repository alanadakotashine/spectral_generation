# spectral_generation
Generate random graphs subject to similar eigenvalues of the symmetric normalized laplacian

Will need the following:

to read gml
https://www.mathworks.com/matlabcentral/fileexchange/45741-read-gml

to compute largest component
https://www.mathworks.com/matlabcentral/fileexchange/30926-largest-component

heap and hash table implementations
https://github.com/armon/c-minheap-indirect
https://troydhanson.github.io/uthash/

To compile mex files run 

mex -v CC=gcc LD=gcc COPTIMFLAGS='-O3 -DNDEBUG' genEqualityConstraints.c
mex -v CC=gcc LD=gcc COPTIMFLAGS='-O3 -DNDEBUG' fiedRoundEigScore.c

To run:

spectralGeneration(graphGML,directoryForOutputPath,numOptSteps,numGraphs,OptFlag)

To run optimization, set OptFlag = 1. 

All graphs will be written as csv files and placed in directoryForOutput along with a text file containing the spectrum of the normalized laplacian of a graph. Input graphs must be in gml format. 

