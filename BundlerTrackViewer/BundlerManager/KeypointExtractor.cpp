/*
 *  KeyPointExtractor.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/16/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "KeyPointExtractor.h"
#include <cmath>

void KeypointExtractor::WriteKeyFile(const char *filename, vector<short> &keys, vector<keypt_t> &info)
{
	FILE *f = fopen(filename, "w");
	fprintf(f, "%d %d\n", info.size(), 128);
	for(int j=0; j<info.size(); j++)
	{
		fprintf(f, "%f %f %f %f %d\n", info[j].y, info[j].x, info[j].scale, info[j].orient, 0 );
		for(int l=0;l<7; l++)
		{
			if(l<6)
			{
				fprintf(f, "%hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd\n", 
						keys[j*128+l*20+0], keys[j*128+l*20+1], keys[j*128+l*20+2], keys[j*128+l*20+3], keys[j*128+l*20+4],
						keys[j*128+l*20+5], keys[j*128+l*20+6], keys[j*128+l*20+7], keys[j*128+l*20+8], keys[j*128+l*20+9],
						keys[j*128+l*20+10], keys[j*128+l*20+11], keys[j*128+l*20+12], keys[j*128+l*20+13], keys[j*128+l*20+14],
						keys[j*128+l*20+15], keys[j*128+l*20+16], keys[j*128+l*20+17], keys[j*128+l*20+18], keys[j*128+l*20+19]);
				
			}
			else
			{
				fprintf(f,"%hd %hd %hd %hd %hd %hd %hd %hd\n",
						keys[j*128+l*20+0], keys[j*128+l*20+1], keys[j*128+l*20+2], keys[j*128+l*20+3],
						keys[j*128+l*20+4], keys[j*128+l*20+5], keys[j*128+l*20+6], keys[j*128+l*20+7]);
			}
		}
	}
	fclose(f);
}

int KeypointExtractor::ReadKeyFile(const char *filename, vector<short> &keys, vector<keypt_t> &info)
{
    FILE *file;
	file = fopen (filename, "r");
    if (! file) 
	{
        // Try to file a gzipped keyfile
        char buf[1024];
        sprintf(buf, "%s.gz", filename);
        gzFile gzf = gzopen(buf, "rb");
		
        if (gzf == NULL) 
		{
            printf("Could not open file: %s\n", filename);
            return 0;
        } 
		else 
		{
            int n = ReadKeysGzip(gzf, keys, info);
            gzclose(gzf);
            return n;
        }
    }
	else
	{
		int n = ReadKeys(file, keys, info);
		fclose(file);
		return n;
	}
}

int KeypointExtractor::ReadKeysGzip(gzFile fp, vector<short> &keys, vector<keypt_t> &info)
{
    int i, num, len;
	
    char header[256];
    gzgets(fp, header, 256);
	
    if (sscanf(header, "%d %d", &num, &len) != 2) 
	{
        printf("Invalid keypoint file.\n");
        return 0;
    }
	
    if (len != 128) 
	{
        printf("Keypoint descriptor length invalid (should be 128).");
        return 0;
    }
	
    keys.resize(num*128, 0);
	info.resize(num);
	
    for (i = 0; i < num; i++) 
	{
        // Allocate memory for the keypoint.
        float x, y, scale, ori;
		int repeated;
        char buf[1024];
        gzgets(fp, buf, 1024);
		
        if (sscanf(buf, "%f %f %f %f\n", &y, &x, &scale, &ori/*, &repeated*/) != 4) 
		{
            printf("Invalid keypoint file format.");
            return 0;
        }
		
		info[i].x = x;
		info[i].y = y;
		info[i].scale = scale;
		info[i].orient = ori;
		if(repeated == 1)
			info[i].repeatedKeypt = true;
		else
			info[i].repeatedKeypt = false;

        for (int line = 0; line < 7; line++) 
		{
            char *str = gzgets(fp, buf, 1024);
			
            if (line < 6) 
			{
                sscanf(buf, "%hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd", 
					   &keys[i*128+line*20+0], &keys[i*128+line*20+1], &keys[i*128+line*20+2], &keys[i*128+line*20+3], &keys[i*128+line*20+4],
					   &keys[i*128+line*20+5], &keys[i*128+line*20+6], &keys[i*128+line*20+7], &keys[i*128+line*20+8], &keys[i*128+line*20+9],
					   &keys[i*128+line*20+10], &keys[i*128+line*20+11], &keys[i*128+line*20+12], &keys[i*128+line*20+13], &keys[i*128+line*20+14],
					   &keys[i*128+line*20+15], &keys[i*128+line*20+16], &keys[i*128+line*20+17], &keys[i*128+line*20+18], &keys[i*128+line*20+19]);
			} 
			else 
			{
                sscanf(buf, "%hd %hd %hd %hd %hd %hd %hd %hd",
					   &keys[i*128+line*20+0], &keys[i*128+line*20+1], &keys[i*128+line*20+2], &keys[i*128+line*20+3], 
					   &keys[i*128+line*20+4], &keys[i*128+line*20+5], &keys[i*128+line*20+6], &keys[i*128+line*20+7]);
            }
        }
    }
	
    return num;
}

int KeypointExtractor::ReadKeys(FILE *fp, vector<short> &keys, vector<keypt_t> &info)
{
    int i, num, len;
	
    if (fscanf(fp, "%d %d", &num, &len) != 2) 
	{
        printf("Invalid keypoint file\n");
        return 0;
    }
	
    if (len != 128) 
	{
        printf("Keypoint descriptor length invalid (should be 128).");
        return 0;
    }
	
    keys.resize(num*128, 0);
	info.resize(num);
	
    for (i = 0; i < num; i++) 
	{
        float x, y, scale, ori;
		int repeated;
        if (fscanf(fp, "%f %f %f %f %d\n", &y, &x, &scale, &ori, &repeated) != 5) 
		{
            printf("Invalid keypoint file format.");
            return 0;
        }
        
		
        info[i].x = x;
		info[i].y = y;
		info[i].scale = scale;
		info[i].orient = ori;
		if(repeated == 1)
			info[i].repeatedKeypt = true;
		else
			info[i].repeatedKeypt = false;
		
        char buf[1024];
        for (int line = 0; line < 7; line++)
		{
            fgets(buf, 1024, fp);
			
            if (line < 6) 
			{
                sscanf(buf, "%hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd %hd", 
					   &keys[i*128+line*20+0], &keys[i*128+line*20+1], &keys[i*128+line*20+2], &keys[i*128+line*20+3], &keys[i*128+line*20+4],
					   &keys[i*128+line*20+5], &keys[i*128+line*20+6], &keys[i*128+line*20+7], &keys[i*128+line*20+8], &keys[i*128+line*20+9],
					   &keys[i*128+line*20+10], &keys[i*128+line*20+11], &keys[i*128+line*20+12], &keys[i*128+line*20+13], &keys[i*128+line*20+14],
					   &keys[i*128+line*20+15], &keys[i*128+line*20+16], &keys[i*128+line*20+17], &keys[i*128+line*20+18], &keys[i*128+line*20+19]);
				
            } 
			else 
			{
                sscanf(buf, "%hd %hd %hd %hd %hd %hd %hd %hd",
					   &keys[i*128+line*20+0], &keys[i*128+line*20+1], &keys[i*128+line*20+2], &keys[i*128+line*20+3], 
					   &keys[i*128+line*20+4], &keys[i*128+line*20+5], &keys[i*128+line*20+6], &keys[i*128+line*20+7]);
				
            }
        }
    }
	
    return num;
}

void KeypointExtractor::formFeaturePairs(std::vector<keypt_t> &keyInfo, std::vector<short> &keys, std::vector<keyPair> &keyPairs, float minDist, float maxDist)
{
	int noKeys = keyInfo.size();
	Vec2f unitX(1.0, 0.0);
	
	for(int i=0; i<noKeys; i++)
	{
		for(int j=i; j<noKeys; j++)
		{
			Vec2f pos1(keyInfo[i].x, keyInfo[i].y);
			Vec2f pos2(keyInfo[j].x, keyInfo[j].y);
			float dist = (pos1-pos2).length();
			float angle = 0.0;
			
			if(dist<maxDist && dist>minDist)
			{
				angle = acos(dot(unitX, (pos2-pos1).normalize()));
				//if((pos2-pos1).normalize()[1] < 0)
				//	angle += M_PI;
				
				//printf("orient:%f %f\n", keyInfo[i].orient, keyInfo[j].orient);
				keyPair p1(keyInfo[i],keyInfo[j], i, j);
				p1.angle = angle;
				keyPairs.push_back(p1);
				
				keyPair p2(keyInfo[j],keyInfo[i], j, i);
				p2.angle = angle;
				//keyPairs.push_back(p2);
			}
		}
	}
}