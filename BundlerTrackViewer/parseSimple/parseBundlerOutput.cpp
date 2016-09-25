#include <fstream>
#include "matrix.h"
#include "Image.h"
#include <vector>
#include <string>
#include <sys/stat.h>
#include <stdio.h> 
#include <assert.h> 

#define NUM_CAMERA_PARAMS 9
#define POLY_INVERSE_DEGREE 6

using namespace std;

typedef struct {
    double R[9];     /* Rotation */
    double t[3];     /* Translation */
    double f;        /* Focal length */
    double k[2];     /* Undistortion parameters */
    double k_inv[POLY_INVERSE_DEGREE]; /* Inverse undistortion parameters */
    char constrained[NUM_CAMERA_PARAMS];
    double constraints[NUM_CAMERA_PARAMS];  /* Constraints (if used) */
    double weights[NUM_CAMERA_PARAMS];      /* Weights on the constraints */
    double K_known[9];  /* Intrinsics (if known) */
    double k_known[5];  /* Distortion params (if known) */

    char fisheye;            /* Is this a fisheye image? */
    char known_intrinsics;   /* Are the intrinsics known? */
    double f_cx, f_cy;       /* Fisheye center */
    double f_rad, f_angle;   /* Other fisheye parameters */
    double f_focal;          /* Fisheye focal length */

    double f_scale, k_scale; /* Scale on focal length, distortion params */
} camera_params_t;

std::vector<camera_params_t> cameraParams;

void writeCam(const char *imageListFile, const char* inputPath, const char *outputPath)
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
	imgDirectory += "/imagesNumbered";
	mkdir(imgDirectory.c_str(), 0770);
	
	string cameras(outputPath);
	cameras += "/cameras.obj";
	
	ofstream camFout(cameras.c_str(), ios::out);
	
	ifstream fin( imageListFile, ios::in);
	
	if( fin.fail() ) 
	{
		printf("Cannot not open image list file %s\n", imageListFile);
		return;
	}
	
	char buffer[1024];
	int count = 0;
	for (int i = 0; i < numCameras; i++) 
	{
		fin.getline(buffer, 1024);
		
		if (cameraParams[i].f == 0.0)
			continue;
		
		string readLine(buffer);
		readLine = readLine.substr(0, readLine.find_first_of(" "));
		string imageFilename(inputPath);
		imageFilename += "/" + readLine;
		
		Img tmpImg;
		tmpImg.read(imageFilename);
		double w = tmpImg.width();
		double h = tmpImg.height();
		
		char buf[256];
		sprintf(buf, "%s/%08d.jpg", imgDirectory.c_str(), count);
		string imgNewFilename(buf);
		tmpImg.write(imgNewFilename);
		
		printf("Image %s has width:%d height:%d\n", imageFilename.c_str(), tmpImg.width(), tmpImg.height());
		
		sprintf(buf, "%s/%08d.txt", camDirectory.c_str(), count);
		FILE *f = fopen(buf, "w");
		assert(f);
		
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
		camFout << "v " << c[0] << " " << c[1] << " " << c[2] << "255 0 0" << endl;
		
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
		
		fclose(f);
		
		count++;
	}
	
	camFout.close();
}

void readBundlerOutputFile(const char *filename, const char *output3D)
{
	ifstream fin( filename, ios::in);
	if( fin.fail() ) 
	{
		printf("Cannot not open bundle file %s\n", filename);
		return;
	}
	
	printf("Reading bundle file:%s\n", filename);
	
	char buffer[1024];
	std::string tmp = "#";
	std::string tmp2 = "v";
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
	ofstream fout(output3D, ios::out);
	//fout << "OFF" << endl;
	//fout << noPoints << " 0 0" << endl;
	//fout << noPoints << endl;
	double x, y, z;
	int cx, cy, cz;
	char longBuffer[10000];
	for(int i=0; i<noPoints; i++)
	{
		fin.getline(buffer, 1024);
		
		fout << "v " << buffer;// << endl;
		//printf("pt:%f %f %f\n", x, y, z);
		//fout << x << " " << y << " " << z << endl;
		
		//fin >> cx >> cy >> cz;
        fin.getline(buffer, 10000);
        fout << " " << buffer << endl;
        
        fin.getline(longBuffer, 10000);
       	/*int noVisViews;
		fin >> noVisViews;
		
		int viewIndex, keyIndex;
		double kpx, kpy;
		for(int j=0; j<noVisViews; j++)
		{
			fin >> viewIndex >> keyIndex >> kpx >> kpy;
		}*/
	}
	fout.close();
}

int main(int argc, const char* argv[])
{
	const char* inputPath = argv[1];
	const char *outputPath = argv[2];
	
	string imageListFile(inputPath);
	imageListFile += "/output/list.txt";
	
	string bundlerFile(inputPath);
	bundlerFile += "/output/bundle.out";
	
	string pointsFile(inputPath);
	pointsFile += "/points3D.obj";
	readBundlerOutputFile(bundlerFile.c_str(), pointsFile.c_str());
	
	
	writeCam(imageListFile.c_str(), inputPath, outputPath);
	return 0;
}