/*
 *  PlyReader.cpp
 *  CathedralProject
 *
 *  Created by Duygu Ceylan on 9/12/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <OpenMesh/Core/IO/IOManager.hh>
#include <OpenMesh/Core/System/omstream.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/Utils/color_cast.hh>
#include "PlyReader.h"

namespace OpenMesh {
	namespace IO {
		
		
		_PLYReader_  PLYReaderInstance;
		_PLYReader_& PLYReader() { return PLYReaderInstance; }
		
		
		_PLYReader_::_PLYReader_() 
		{ 
			IOManager().register_module(this); 
		}
		
		
		bool _PLYReader_::read(const std::string& filename, BaseImporter& bi, Options& opt)
		{
			ifstream fin( filename.c_str(), ifstream::in);
			if( fin.fail() ) 
			{
				omerr() << "[PLYReader] : cannot not open file " << filename << std::endl;
				return false;
				
			}
			
			{
#if defined(WIN32)
				string::size_type dot = filename.find_last_of("\\/");
#else
				string::size_type dot = filename.rfind("/");
#endif
				path = (dot == string::npos) ? "./" : string(filename.substr(0,dot+1));
			}
			
			bool ok = read(fin, bi, opt);
			
			fin.close();
			return ok;
		}
		
		bool _PLYReader_::read(istream& fin, BaseImporter& bi, Options& opt)
		{
			omlog() << "[PLYReader] : read file\n";
			
			
			BaseImporter::VHandles vhandles, faceVertexHandles;
			std::vector<Vec3f>     normals;
			std::vector<Vec3uc>     colors;
			
			string tag, tmp;
			char buffer[1024];
			fin.getline(buffer,1024);
			fin.getline(buffer, 1024);
			
			int noListElements;
			
			if( strcmp(buffer, "format binary_little_endian 1.0") != 0 )
			{
				omerr() << "[PLYReader] doesn't know this format: " << buffer <<std::endl;
				return false;
			}
			
			tmp = "element";
			do 
			{
				fin.getline(buffer, 1024);
				if( fin.fail() ) 
				{
					omerr() << "[PLYReader] broken file" << std::endl;
					return false;
				}
			} 
			while( strncmp(buffer,"element", tmp.length()) != 0 ) ;
			
			int nver = 0;
			int nface = 0;
			tmp = "element vertex";
			if( strncmp(buffer, "element vertex", tmp.length()) == 0 )
				sscanf(buffer, "element vertex %d", &nver);
			omlog() << "[PLYReader] number of vertices: " << nver << std::endl;
			
			bool bver=false, bnrm=false, bcol=false, bAlpha = false;
			char str[1024], type[1024];
			tmp = "end_header";
			do 
			{
				fin.getline(buffer,1024);
				if( strncmp(buffer, "element face", tmp.length()) == 0 )
				{
					sscanf(buffer, "element face %d", &nface);
					continue;
				}
				
				sscanf(buffer,"property %s %s",type,str);
				if( !strcmp(str,"x")  || !strcmp(str,"y") || !strcmp(str,"z") ) 
					bver=true;
				if( !strcmp(str,"nx") || !strcmp(str,"ny") || !strcmp(str,"nz") ) 
				{
					bnrm=true;
					opt += Options::VertexNormal;
				}
				if(!strcmp(str,"red") ||
				   !strcmp(str,"green") ||
				   !strcmp(str,"blue") || !strcmp(str,"diffuse_red") ||
				   !strcmp(str,"diffuse_green") ||
				   !strcmp(str,"diffuse_blue") ) 
				{
					bcol=true;
					opt += Options::VertexColor;
				}
				if(!strcmp(str,"alpha"))
				{
					bAlpha = true;
				}
			} 
			while( strncmp(buffer,"end_header", tmp.length()) != 0 );
			
			if( !bver) 
			{
				omerr() <<"[PLYReader] nrm   aint exist wo ver" << std::endl;
				return false;
			}
			
			float tempPos[3], tempNormal[3];
			unsigned char tempChar[4];
			int tempInt;
			unsigned char tempUChar;
			for(int i=0; i<nver; i++)
			{
				if( bver ) 
				{
					fin.read((char*)tempPos, sizeof(float)*3 );
					vhandles.push_back( bi.add_vertex(OpenMesh::Vec3f(tempPos[0],tempPos[1],tempPos[2])) );
				}
				if( bnrm ) 
				{
					fin.read((char*)tempNormal, sizeof(float)*3 );
					bi.set_normal(vhandles.back(),OpenMesh::Vec3f(tempNormal[0],tempNormal[1],tempNormal[2]).normalize());
				}
				if( bcol ) 
				{
					if(bAlpha)
					{
						fin.read((char*)tempChar, sizeof(unsigned char)*4 );					
					}
					else
					{
						fin.read((char*)tempChar, sizeof(unsigned char)*3 );
					}
					bi.set_color(vhandles.back(),OpenMesh::Vec3uc(tempChar[0], tempChar[1], tempChar[2]));
				}
				
			}
			
			for(int i=0; i<nface; i++)
			{
				fin.read((char*)&tempUChar, sizeof(unsigned char));
				noListElements = (int)tempUChar;
				faceVertexHandles.clear();
				for(int j=0; j<noListElements; j++)
				{
					fin.read((char*)&tempInt, sizeof(int));
					faceVertexHandles.push_back(vhandles[tempInt]);
				}
				bi.add_face(faceVertexHandles);
			}
			
			return true;
		}
	}
}
