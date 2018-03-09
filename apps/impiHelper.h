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
