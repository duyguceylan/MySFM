Below is an explanation of the arguments for the 'solveLinear' matlab function. You have to call this function once you have built your image matching graph with candidate alignments and want to optimize for the edge costs of this graph.

Input Arguments:

***noImages***: number of images

***alignmentsFile***: the current estimate of pairwise image alignments based on feature
matches

format:
First line is the number of pairwise image alignment. Then for each image alignment we have the following format:
				
<index of first image> <index of second image> <weight> <isGrid boolean>
<dummy number 1><dummy number 2><dummy number 3>
<number of grid alignments>
	for each grid alignment:
	<template index><plane index 1><plane index 2><grid index 1><grid index 2><column shift><row shift>
	<rotation>
	
'weight' is between 0 and 1 and measures the quality of the alignment, in our case it is computed as described in Section of our paper. 'isGrid' is a number that is either a 0 or 1 encoding if the alignment is encoded as grid shifts (1) or not (0). For each grid alignment, 'template index' is the index of the user marked template (0-based). It is 0 if only one template has been marked. 'plane index 1' is the index of the facade plane that contains the grid (0-based), 'plane index 2' is defined similarly for the second grid. Plane indices are not used in the algorithm. 'grid index 1' is the index of the grid among all the grids detected in the first image (0-based) and 'grid index 2' is defined similarly for the second image. 'column shift' and 'row shift' encode the number of rows and columns the grid in the second image needs to be shifted over the first image as specified by the alignment.

Below is a simple example file:

2									//2 alignments
0 1 0.819033 1						//alignment between images 0 and 1
0.995718 293.365 -116.771
1									//1 grid alignment
0 0 0 0 0 0 -1
0.998718 0.0159908 -0.0480297		//rotation matrix
-0.015744 0.999861 0.00551244
0.0481112 -0.00474919 0.998831
0 2 0.776911 1						alignment between images 0 and 2
0.987741 -91.6972 -7.98977
1									//1 grid alignment
0 0 0 0 0 0 -1
0.998603 0.0128446 -0.0512509		//rotation matrix
-0.0132547 0.999883 -0.00766929
0.0511464 0.00833789 0.998656

				
***cycleFile***: list of three-cycles in the image matching graph

format: each line contains three integers which are the indices of the images				in a cycle, indices are zero based
				
for example the following line means the images 0, 1, and 3 are in a cycle:
		...
		0 1 3
		...

***newAlignmentsFile****: Contains the refined pairwise image alignments, has a format similar to the input alignments file.

First line is the number of pairwise image alignments. Then for each image alignment we have the following format:
				
<index of first image> <index of second image> <cost assigned by optimization> <weight of the alignment>
<number of grid alignments>
	for each grid alignment:
	<template index><plane index 1><plane index 2><grid index 1><grid index 2><column shift><row shift>
	<rotation>

Plane indices output by the algorithm are all 0 since this variable is not really used in the optimization. The costs assigned to the pairwise alignments are in the range (0.1, 1), 0.1 denoting a good alignment. The minimum spanning tree of the image matching graph based on the edge costs give the correct alignments that register all the images. The file 'parser.cpp' includes code to parse this file.

Below is an example output file:

2								//2 alignments
0 1 0.100000 0.449367			//alignment between images 0 and 1 has a cost 0.1
1								//1 grid alignment
1 0 0 0 0 0 0					//template id 1, grid indices are 0 and 0 row/column shift
0.999415 -0.006731 -0.033532	//rotation
0.008717 0.998193 0.059446
0.033071 -0.059704 0.997668
0 2 0.100000 0.278247			//alignment between images 0 and 2 has a cost 0.1
2								//2 grid alignments
0 0 0 0 0 1 0					//template id 0, grid indices are 0
1 0 0 0 0 1 0					//template id 1, grid indices are 0
0.996764 0.001849 -0.080367		//rotation
0.003860 0.997482 0.070821
0.080296 -0.070902 0.994246