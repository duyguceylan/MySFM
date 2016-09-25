/*
 *  BundlerOutputParser.cpp
 *  BundlerTrackViewer
 *
 *  Created by Duygu Ceylan on 11/8/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <fstream>
#include <sys/stat.h>
#include "matrix.h"
#include "BundlerOutputParser.h"
#include "PointCloud.h"

void BundlerOutputParser::readBundlerOutputFile(const char *filename, int width, int height)
{
	ifstream fin( filename, ios::in);
	if( fin.fail() ) 
	{
		printf("Cannot not open bundle file %s\n", filename);
		return;
	}
	
	printf("Reading bundle file:%s\n", filename);
	
	char buffer[1024];
	string tmp = "#";
	string tmp2 = "v";
	do 
	{
		fin.getline(buffer, 1024);
		if( fin.fail() ) 
		{
			printf("Broken bundle file.\n");
			return ;
		}
	} 
	while( strncmp(buffer,"#", tmp.length()) == 0  || strncmp(buffer,"v", tmp2.length()) == 0 ) ;
	
	int noCameras, noPoints;
	sscanf(buffer,"%d %d",&noCameras,&noPoints);
	
	printf("%d cameras and %d points found.\n", noCameras, noPoints);
	
	//read cameras
	float f, k1, k2;
	float r1, r2, r3, r4, r5, r6, r7, r8, r9;
	float t1, t2, t3;
	for(int i=0; i<noCameras; i++)
	{
		camera_params_t cam;
		fin.getline(buffer, 1024);
		sscanf(buffer,"%f %f %f",&f,&k1,&k2);
		fin.getline(buffer, 1024);
		sscanf(buffer,"%f %f %f",&r1,&r2,&r3);
		fin.getline(buffer, 1024);
		sscanf(buffer,"%f %f %f",&r4,&r5,&r6);
		fin.getline(buffer, 1024);
		sscanf(buffer,"%f %f %f",&r7,&r8,&r9);
		fin.getline(buffer, 1024);
		sscanf(buffer,"%f %f %f",&t1,&t2,&t3);
		cam.f = f;
        cam.k[0] = k1;
        cam.k[1] = k2;
		cam.R[0] = r1; cam.R[1] = r2; cam.R[2] = r3;
		cam.R[3] = r4; cam.R[4] = r5; cam.R[5] = r6;
		cam.R[6] = r7; cam.R[7] = r8; cam.R[8] = r9;
		cam.t[0] = t1; cam.t[1] = t2; cam.t[2] = t3;
		
		cameraParams.push_back(cam);
        
	}
	
	//read points
	float centerX = width/2.0;
	float centerY = height/2.0;
	
	Vec3f pos;
	Vec3d col;
	int noVisViews;
	int viewIndex;
	int keyIndex;
	Vec2f pos2D;
    
    PointCloud *pc = new PointCloud();
    
	for(int i=0; i<noPoints; i++)
	{
		CorrespondenceTrack track;
		
		fin >> pos[0] >> pos[1] >> pos[2];
		
		fin >> col[0] >> col[1] >> col[2];
        
        Vec3uc c;
        c[0] = col[0];
        c[1] = col[1];
        c[2] = col[2];
        
        //printf("****************\n");
        //printf("%f %f %f\n", pos[0], pos[1], pos[2]);
        pc->addPoint(pos, Vec3f(0,0,0), c);
		
		//view list
		fin >> noVisViews;
		
		track.pos3D = pos;
		track.color = col;
		
		for(int j=0; j<noVisViews; j++)
		{
			Correspondence c;
			fin >> viewIndex >> keyIndex >> pos2D[0] >> pos2D[1];
			pos2D[0] = (pos2D[0] + centerX);
			pos2D[1] = ((-pos2D[1]) + centerY);
			c.imageIndex = viewIndex;
			c.keyIndex = keyIndex;
			c.position = pos2D;
			track.correspondences.push_back(c);
		}
		
		tracks.push_back(track);
	}
    
    pc->writeMesh("/Users/ceylan/Desktop/sparse.ply");
}

void BundlerOutputParser::undistortImage(Img &img, int w, int h, const std::string &out, const camera_params_t &camera)
{ 
	//printf("Undistorting image\n");
	fflush(stdout);
	
	Img imgOut(w, h);
	
	double f2Inv = 1.0 / (camera.f * camera.f);
	
	for (int y = 0; y < h; y++) 
	{
		for (int x = 0; x < w; x++) 
		{
			double xC = x - 0.5 * w;
			double yC = y - 0.5 * h;
            
			double r2 = (xC * xC + yC * yC) * f2Inv;
			double factor = 1.0 + camera.k[0] * r2 + camera.k[1] * r2 * r2;
            
			xC *= factor;
			yC *= factor;
			
			xC += 0.5 * w;
			yC += 0.5 * h;
            
			Color c;
			if (xC >= 0 && xC <= w - 1 && yC >= 0 && yC <= h - 1) 
			{
				c = img.lerp(xC, yC);
			} 
			else 
			{
				c = Color(0.0, 0.0, 0.0);
			}
            
			imgOut.setColor(x, y, c);
		}
	}
	
	imgOut.write(out);
}

void BundlerOutputParser::writeCamAndImage(const char *outputPath, std::vector<string> &imagePaths, std::vector<PerspectiveCamera> &cameras)
{
	int numCameras = (int) cameraParams.size();
	
	//create the output directory
	mkdir(outputPath, 0770);
	
	//create the camera directory
	string camDirectory(outputPath);
	camDirectory += "/cam";
	mkdir(camDirectory.c_str(), 0770);
	
	//create the image directory
	string imgDirectory(outputPath);
	imgDirectory += "/images";
	mkdir(imgDirectory.c_str(), 0770);
	
	int count = 0;
	for (int i = 0; i < numCameras; i++) 
	{
		if (cameraParams[i].f == 0.0)
			continue;
		
		char buf[256];
		sprintf(buf, "%s/%08d.txt", camDirectory.c_str(), count);
		FILE *f = fopen(buf, "w");
		assert(f);
		
		char imgBuf[256];
		string extension = imagePaths[i].substr(imagePaths[i].find(".")+1);
		sprintf(imgBuf, "%s/%08d.%s", imgDirectory.c_str(), count, extension.c_str());
		
		Img img;
		img.read(imagePaths[i]);
		int w = img.width();
		int h = img.height();
		
		// Compute the projection matrix
		double focal = cameraParams[i].f;
		double *R = cameraParams[i].R;
		double *t = cameraParams[i].t;
	
		double *Rt = new double[9];
		double *c = new double[3];
		matrix_transpose(3, 3, R, Rt);
		matrix_scale(3,3,Rt, -1.0, Rt);
		matrix_product(3, 3, 3, 1, Rt, t, c);
		
		//printf("****cam%d******\n", i);
		//printf("focal:%f\n", focal);
		//printf("R:%f %f %f %f %f %f %f %f %f\n", R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8]);
		//printf("t:%f %f %f\n", t[0], t[1], t[2]);
		printf("c:%f %f %f\n", c[0], c[1], c[2]);
		
		double K[9] = { -focal, 0.0, 0.5 * w - 0.5,
			0.0, focal, 0.5 * h - 0.5,
			0.0, 0.0, 1.0 };
		
		double Ptmp[12] = { R[0], R[1], R[2], t[0],
			R[3], R[4], R[5], t[1],
			R[6], R[7], R[8], t[2] };
        
		double P[12];
		matrix_product(3, 3, 3, 4, K, Ptmp, P);
		matrix_scale(3, 4, P, -1.0, P);
		
		fprintf(f, "CONTOUR\n");
		fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[0], P[1], P[2],  P[3]);
		fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[4], P[5], P[6],  P[7]);
		fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[8], P[9], P[10], P[11]);
		
		Matrix4f projMatrix;
		projMatrix[0][0] = P[0]; projMatrix[1][0] = P[1]; projMatrix[2][0] = P[2]; projMatrix[3][0] = P[3];
		projMatrix[0][1] = P[4]; projMatrix[1][1] = P[5]; projMatrix[2][1] = P[6]; projMatrix[3][1] = P[7];
		projMatrix[0][2] = P[8]; projMatrix[1][2] = P[9]; projMatrix[2][2] = P[10]; projMatrix[3][2] = P[11];
		projMatrix[0][3] = 0.0; projMatrix[1][3] = 0.0; projMatrix[2][3] = 0.0; projMatrix[3][3] = 1.0;
		PerspectiveCamera cam;
		cam.setProjectionMatrix(projMatrix);
		cameras.push_back(cam);
		
		fclose(f);
		
		//undistortImage(img, w, h, imgBuf, cameraParams[i]);
		
		count++;
	}
	
}

void BundlerOutputParser::readCalibration(const char *outputPath, int noImages, std::vector<PerspectiveCamera> &cameras)
{
	int maxLevel = 3;
	string camDirectory(outputPath);
	camDirectory += "/cam";
    
	for(int i=0; i<noImages; i++)
	{
        PerspectiveCamera cam;
        
        char buf[256];
        sprintf(buf, "%s/%08d.txt", camDirectory.c_str(), i);
        string fname(buf);
        ifstream fin(fname.c_str());
		
        if(fin)
        {
            cam.init(fname, maxLevel);
        }
		
		cameras.push_back(cam);
	}
}

void BundlerOutputParser::convertSfMOutputToPMVS(const char *input)
{
	ifstream fin(input, ios::in);
	
	string filename(input);
	filename = filename.substr(0, filename.find_last_of("/")+1);
	filename += "cam";
	//create the output directory
	mkdir(filename.c_str(), 0770);
	
	int noCameras;
	fin >> noCameras;
	Matrix3f K;
	K.setIdentity();
	K[0][0] = 2945.455;
	K[1][1] = 2945.455;
	K[2][0] = 1440;
	K[2][1] = 1080;
	for(int i=0; i<noCameras; i++)
	{
		int index;
		Matrix3f R;
		Vec3f t;
		fin >> index;
		fin >> R[0][0] >> R[1][0] >> R[2][0] >> t[0];
		fin >> R[0][1] >> R[1][1] >> R[2][1] >> t[1];
		fin >> R[0][2] >> R[1][2] >> R[2][2] >> t[2];
		
		R = K*R;
		t = K*t;
		
		char buf[256];
		sprintf(buf, "%s/%08d.txt", filename.c_str(), i);
		
		ofstream fout(buf, ios::out);
		fout << "CONTOUR\n";
		fout << R[0][0] << " " << R[1][0] << " " << R[2][0] << " " << t[0] << endl;
		fout << R[0][1] << " " << R[1][1] << " " << R[2][1] << " " << t[1] << endl;
		fout << R[0][2] << " " << R[1][2] << " " << R[2][2] << " " << t[2] << endl;
		fout.close();
	}
}

/*void convertPMVSOutputToBundler(vector<PerspectiveCamera> &cameras, Matrix3f intrinsic)
{
	double focal = cameraParams[i].f;
	double *R = cameraParams[i].R;
	double *t = cameraParams[i].t;
	
	double *Rt = new double[9];
	double *c = new double[3];
	matrix_transpose(3, 3, R, Rt);
	matrix_scale(3,3,Rt, -1.0, Rt);
	matrix_product(3, 3, 3, 1, Rt, t, c);
	
	printf("****cam%d******\n", i);
	printf("focal:%f\n", focal);
	printf("R:%f %f %f %f %f %f %f %f %f\n", R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8]);
	printf("t:%f %f %f\n", t[0], t[1], t[2]);
	printf("c:%f %f %f\n", c[0], c[1], c[2]);
	
	double K[9] = { -focal, 0.0, 0.5 * w - 0.5,
		0.0, focal, 0.5 * h - 0.5,
		0.0, 0.0, 1.0 };
	
	double Ptmp[12] = { R[0], R[1], R[2], t[0],
		R[3], R[4], R[5], t[1],
		R[6], R[7], R[8], t[2] };
	
	double P[12];
	matrix_product(3, 3, 3, 4, K, Ptmp, P);
	matrix_scale(3, 4, P, -1.0, P);
	
	fprintf(f, "CONTOUR\n");
	fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[0], P[1], P[2],  P[3]);
	fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[4], P[5], P[6],  P[7]);
	fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[8], P[9], P[10], P[11]);
	
	Matrix4f projMatrix;
	projMatrix[0][0] = P[0]; projMatrix[1][0] = P[1]; projMatrix[2][0] = P[2]; projMatrix[3][0] = P[3];
	projMatrix[0][1] = P[4]; projMatrix[1][1] = P[5]; projMatrix[2][1] = P[6]; projMatrix[3][1] = P[7];
	projMatrix[0][2] = P[8]; projMatrix[1][2] = P[9]; projMatrix[2][2] = P[10]; projMatrix[3][2] = P[11];
	projMatrix[0][3] = 0.0; projMatrix[1][3] = 0.0; projMatrix[2][3] = 0.0; projMatrix[3][3] = 1.0;
	PerspectiveCamera cam;
	cam.setProjectionMatrix(projMatrix);
	cameras.push_back(cam);
	
}*/

