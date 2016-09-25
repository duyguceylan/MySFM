/*
 *  PlyWriter.h
 *  CathedralProject
 *
 *  Created by Duygu Ceylan on 9/14/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __PLYWRITER_HH__
#define __PLYWRITER_HH__

#include <string>
#include <stdio.h>
#include <iostream>
#include <OpenMesh/Core/System/config.hh>
#include <OpenMesh/Core/Utils/SingletonT.hh>
#include <OpenMesh/Core/IO/exporter/BaseExporter.hh>
#include <OpenMesh/Core/IO/writer/BaseWriter.hh>

using namespace std;

namespace OpenMesh {
	namespace IO {
		
		
		class OPENMESHDLLEXPORT _PLYWriter_ : public BaseWriter
		{
		public:
			
            _PLYWriter_();
			
            /// Destructor
            virtual ~_PLYWriter_() {};
            
			std::string get_description() const  { return "Stanford Polygon"; }
			std::string get_extensions()  const  { return "ply"; }
			
			bool write(const std::string&, BaseExporter&, Options, std::streamsize _precision = 6) const;
			bool write(std::ostream&, BaseExporter&, Options, std::streamsize _precision = 6) const;
            
			size_t binary_size(BaseExporter&, Options) const { return 0; }
			
		private:
			
			bool writeMaterial(ostream&, BaseExporter&, Options) const;
		};
		
		
		extern _PLYWriter_  __PLYWriterinstance;
		OPENMESHDLLEXPORT _PLYWriter_& PLYWriter();
		
	}
}

#endif
