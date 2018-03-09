// ======================================================================== //
// Copyright 2009-2017 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#pragma once

#include <ospcommon/vec.h>
#include <ospcommon/box.h>
#include <stdexcept>
#include <sstream>
#include <type_traits>

/*! _everything_ in the ospray core universe should _always_ be in the
  'ospray' namespace. */
namespace ospray {

  /*! though not required, it is good practice to put any module into
    its own namespace (isnide of ospray:: ). Unlike for the naming of
    library and init function, the naming for this namespace doesn't
    particularlly matter. E.g., 'bilinearPatch', 'module_blp',
    'bilinar_patch' etc would all work equally well. */
  namespace impi {

    // use ospcommon for vec3f etc
    using namespace ospcommon;
    
    /*! helper class to parse command-line arguments */
    struct CommandLine {
      CommandLine() = default;
      CommandLine(int ac, const char **av) { Parse(ac, av); }
      void Parse(int ac, const char **av);
      std::vector<std::string> inputFiles;
    };

    inline void CommandLine::Parse(int ac, const char **av)
    {
      for (int i=1;i<ac;i++) {
        const std::string arg = av[i];
        if (arg[0] == '-') {
          throw std::runtime_error("un-handled cmdline argument '"+arg+"'");
        } else {
          // no arg: must be an input file
          inputFiles.push_back(arg);
        }
      }
    }
    
    template <class T1, class T2> T1 lexical_cast(const T2& t2) {
      std::stringstream s; s << t2; T1 t1;
      if(s >> t1 && s.eof()) { return t1; }
      else {
	throw std::runtime_error("bad conversion " + s.str());
	return T1();
      }
    }

    template<int N, typename T> T Parse(const int ac, const char** av, int &i, T& v) {
      const int init = i;
      if (init + N < ac) {	
	if constexpr (std::is_scalar<T>::value) {
	  v = lexical_cast<double, const char*>(av[i+1]);
	} else {	  
	  for (int k = 0; k < N; ++k) {
	    v[k] = lexical_cast<double, const char*>(av[i+1+k]);
	  }
	}
      } else {
	throw std::runtime_error(std::to_string(N) + " values required for " + av[init]);
      }
    }
    
  } // ::ospray::bilinearPatch
} // ::ospray
