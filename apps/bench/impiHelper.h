// ======================================================================== //
// Copyright SCI Institute, University of Utah, 2018
// ======================================================================== //

#pragma once

#include <stdint.h>
#include <stdio.h>
#include <errno.h>
#ifdef _WIN32
#  include <malloc.h>
#else
#  include <alloca.h>
#endif
#include <chrono>
#include <sstream>
#include <type_traits>

namespace ospray {
  namespace impi {

    // helper function to write the rendered image as PPM file
    inline void writePPM(const char *fileName,
			 const size_t sizex, 
			 const size_t sizey,
			 const uint32_t *pixel)
    {
      FILE *file = fopen(fileName, "wb");
      if (!file) {
	fprintf(stderr, "fopen('%s', 'wb') failed: %d", fileName, errno);
	return;
      }
      fprintf(file, "P6\n%i %i\n255\n", sizex, sizey);
      unsigned char *out = (unsigned char *)alloca(3 * sizex);
      for (int y = 0; y < sizey; y++) {
	const unsigned char *in = (const unsigned char *) &pixel[(sizey - 1 - y) * sizex];
	for (int x = 0; x < sizex; x++) {
	  out[3 * x + 0] = in[4 * x + 0];
	  out[3 * x + 1] = in[4 * x + 1];
	  out[3 * x + 2] = in[4 * x + 2];
	}
	fwrite(out, 3 * sizex, sizeof(char), file);
      }
      fprintf(file, "\n");
      fclose(file);
    }
    inline void writePPM(const std::string &fileName,
			 const size_t sizex, 
			 const size_t sizey,
			 const uint32_t *pixel)
    {
      writePPM(fileName.c_str(), sizex, sizey, pixel);
    }

    // timer
    typedef std::chrono::high_resolution_clock::time_point time_point;
    time_point Time() {
      return std::chrono::high_resolution_clock::now();
    }
    double Time(const time_point& t1) {
      time_point t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> et = 
	std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
      return et.count();  
    }

    template <class T1, class T2> T1 lexical_cast(const T2& t2) {
      std::stringstream s; s << t2; T1 t1;
      if(s >> t1 && s.eof()) { return t1; }
      else {
	throw std::runtime_error("bad conversion " + s.str());
	return T1();
      }
    }

    // retrieve input arguments
    // because constexper cannot be used on gcc 4.8, a hack is used here
    template<typename T> T ParseScalar(const int ac, const char** av, int &i, T& v) {
      const int init = i;
      if (init + 1 < ac) {
	v = lexical_cast<double, const char*>(av[++i]);
      } else {
	throw std::runtime_error("value required for " + std::string(av[init]));
      }
    }
    template<int N, typename T> T Parse(const int ac, const char** av, int &i, T& v) {
      const int init = i;
      if (init + N < ac) {
	for (int k = 0; k < N; ++k) {
	  v[k] = lexical_cast<double, const char*>(av[++i]);
	}	
      } else {
	throw std::runtime_error(std::to_string(N) + " values required for " + av[init]);
      }
    }
    template<> int Parse<1, int>(const int ac, const char** av, int &i, int& v) {
      return ParseScalar<int>(ac, av, i, v);
    }
    template<> float Parse<1, float>(const int ac, const char** av, int &i, float& v) {
      return ParseScalar<float>(ac, av, i, v);
    }
    template<> double Parse<1, double>(const int ac, const char** av, int &i, double& v) {
      return ParseScalar<double>(ac, av, i, v);
    }

  };
};
