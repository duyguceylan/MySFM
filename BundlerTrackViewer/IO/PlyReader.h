/*
 *  PlyReader.h
 *  CathedralProject
 *
 *  Created by Duygu Ceylan on 9/12/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef __PLYREADER_HH__
#define __PLYREADER_HH__

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <map>
#include <stdio.h>

#include <OpenMesh/Core/System/config.hh>
#include <OpenMesh/Core/Utils/SingletonT.hh>
#include <OpenMesh/Core/IO/importer/BaseImporter.hh>
#include <OpenMesh/Core/IO/reader/BaseReader.hh>

using namespace std;

namespace OpenMesh {
	namespace IO {
		
		
		class _PLYReader_ : public BaseReader
		{
		public:
			
			_PLYReader_();
			
			virtual ~_PLYReader_() { }
			
			string get_description() const { return "Stanford Polygon"; }
			string get_extensions()  const { return "ply"; }
			
			bool read(const string& filename, BaseImporter& bi, Options& opt);
			
		private:
			
			bool read(istream &fin, BaseImporter& bi, Options& opt);
			
			string path;
			
		};
		
		
		extern _PLYReader_  PLYReaderInstance;
		_PLYReader_& PLYReader();
		
	}
}
#endif
