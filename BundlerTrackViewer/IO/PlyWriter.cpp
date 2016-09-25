/*
 *  PlyWriter.cpp
 *  CathedralProject
 *
 *  Created by Duygu Ceylan on 9/14/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "PlyWriter.h"

#include <fstream>

// OpenMesh
#include <OpenMesh/Core/System/config.hh>
#include <OpenMesh/Core/IO/BinaryHelper.hh>
#include <OpenMesh/Core/IO/writer/OBJWriter.hh>
#include <OpenMesh/Core/IO/IOManager.hh>
#include <OpenMesh/Core/System/omstream.hh>


namespace OpenMesh {
	namespace IO {
		
		
		_PLYWriter_  __PLYWriterinstance;
		_PLYWriter_& PLYWriter() { return __PLYWriterinstance; }
		
		
		_PLYWriter_::_PLYWriter_() { IOManager().register_module(this); }
		
		bool _PLYWriter_::write(const std::string& filename, BaseExporter& be, Options opt, std::streamsize _precision) const
		{
			ofstream fout( filename.c_str(), ifstream::out);
			if( fout.fail() )
			{
				omerr() << "[PLYWriter] : cannot not open file " << filename << std::endl;
				return false;
				
			}
			
			bool result = writeMaterial(fout, be, opt);
			
			fout.close();
			return result;
		}
		
		bool _PLYWriter_::write(std::ostream& os, BaseExporter& be, Options opt, std::streamsize _precision) const
        {
            bool result = writeMaterial(os, be, opt);
            return result;
        }
        
		bool _PLYWriter_::writeMaterial(ostream & fout, BaseExporter& be, Options opt) const
		{
			unsigned int i, nV;
			Vec3f v, n;
			Vec3uc c;
			VertexHandle vh;
			
			omlog() << "[PLYWriter] : write file\n";
			
			// check exporter features
			if (!check( be, opt))
				return false;
			
			nV=be.n_vertices();
			
			//header
			fout << "ply" << endl;
			fout << "format binary_little_endian 1.0" << endl;
			fout << "element vertex " << nV << endl;
			fout << "property float x" << endl;
			fout << "property float y" << endl;
			fout << "property float z" << endl;
			
			if(opt.check(Options::VertexNormal))
			{
				fout << "property float nx" << endl;
				fout << "property float ny" << endl;
				fout << "property float nz" << endl;
			}
			if(opt.check(Options::VertexColor))
			{
				fout << "property uchar red" << endl;
				fout << "property uchar green" << endl;
				fout << "property uchar blue" << endl;
			}
			fout << "end_header" << endl;
			
			
			if(opt.check(Options::VertexNormal) && opt.check(Options::VertexColor))
			{
				for (i=0; i<nV; ++i)
				{
					vh = VertexHandle(i);
					v  = be.point(vh);
					n = be.normal(vh);
					c  = be.color (vh);
					
					fout.write((char*)&v[0], sizeof(float)*3);
					fout.write((char*)&n[0], sizeof(float)*3);
					fout.write((char*)&c[0], sizeof(unsigned char)*3);
				}
			}
			else if(opt.check(Options::VertexNormal))
			{
				for (i=0; i<nV; ++i)
				{
					vh = VertexHandle(i);
					v  = be.point(vh);
					n = be.normal(vh);
					
					fout.write((char*)&v[0], sizeof(float)*3);
					fout.write((char*)&n[0], sizeof(float)*3);
				}
			}
			else if(opt.check(Options::VertexColor))
			{
				for (i=0; i<nV; ++i)
				{
					vh = VertexHandle(i);
					v  = be.point(vh);
					c  = be.color (vh);
					fout.write((char*)&v[0], sizeof(float)*3);
					fout.write((char*)&c[0], sizeof(unsigned char)*3);
				}
			}
			else
			{
				for (i=0; i<nV; ++i)
				{
					vh = VertexHandle(i);
					v  = be.point(vh);
					fout.write((char*)&v[0], sizeof(float)*3);
				}
			}
			
			return true;
		}
		
		
	}
}