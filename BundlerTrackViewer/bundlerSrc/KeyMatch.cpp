/* 
 *  Copyright (c) 2008-2010  Noah Snavely (snavely (at) cs.cornell.edu)
 *    and the University of Washington
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 */

/* KeyMatch.cpp */
/* Read in keys, match, write results to a file */

#include <time.h>

#include "keys.h"

class KeyMatch
{
	/*static void RunKeyMatch(char *keys1_in, char *keys2_in, char *file_out)
	{
		double ratio;
    
		ratio = 0.6; // atof(argv[3]);
		
		clock_t start = clock();

		short *keys1, *keys2;

		int num1 = ReadKeyFile(keys1_in, &keys1);
		int num2 = ReadKeyFile(keys2_in, &keys2);

		clock_t end = clock();    
		printf("Reading keys took %0.3fs\n", (end - start) / ((double) CLOCKS_PER_SEC));

		// Compute likely matches between two sets of keypoints
		std::vector<KeypointMatch> matches = 
		MatchKeys(num1, keys1, num2, keys2, ratio);

		int num_matches = (int) matches.size();
		// int num_matches_sym = (int) matches_sym.size();

		printf("num_matches = %d\n", num_matches);
		// printf("num_matches_sym = %d\n", num_matches_sym);

		if (num_matches >= 16) {
			FILE *f = fopen(file_out, "w");
		
			// Write the number of matches 
			fprintf(f, "%d\n", (int) matches.size());

			for (int i = 0; i < num_matches; i++) {
				fprintf(f, "%d %d\n", matches[i].m_idx1, matches[i].m_idx2);
			}

			fclose(f);
		}
	}*/
};
