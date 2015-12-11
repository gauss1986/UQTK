/* =====================================================================================
                      The UQ Toolkit (UQTk) version 2.1.1
                     Copyright (2013) Sandia Corporation
                      http://www.sandia.gov/UQToolkit/

     Copyright (2013) Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
     with Sandia Corporation, the U.S. Government retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is free software: you can redistribute it and/or modify
     it under the terms of the GNU Lesser General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
 
     You should have received a copy of the GNU Lesser General Public License
     along with UQTk.  If not, see <http://www.gnu.org/licenses/>.
 
     Questions? Contact Bert Debusschere <bjdebus@sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */

#include <stdlib.h>
#include <cstdio>
#include <string>
#include <iostream>
#include "UtilsDuffing.h"
#include "Array1D.h"

void WriteToFile(Array1D<double>& data, const char* filename)
{
  int nx=data.XSize();
  
  FILE* f_out;
  if(!(f_out = fopen(filename,"w"))){ 
    printf("WriteToFile: could not open file '%s'\n",filename); 
    exit(1); 
  }

 for(int ix = 0 ; ix < nx ; ix++){
   fprintf(f_out, "%lg\n", data(ix));
 }

 if(fclose(f_out)){ 
   printf("WriteToFile: could not close file '%s'\n",filename); 
   exit(1); 
 }

 printf("Data written to '%s'\n", filename);

 return; 
}

void WriteModesToFilePtr(const double tym, const double* u, const int n, FILE* f_dump)
{
  // Write time
  fprintf(f_dump, "%lg ", tym);

  // Write modes
  for (int ip=0; ip < n; ip++){
    fprintf(f_dump, "%lg ", u[ip]);
  }

  fprintf(f_dump, "\n");
  
  return;
}

void WriteMeanStdDevToFilePtr_lorenz(const double tym, const double x, const double y, const double z,  FILE* f_dump)
{
  // Write time
  fprintf(f_dump, "%lg ", tym);

  // Write x, y, z
  fprintf(f_dump, "%lg %lg %lg", x,y,z);

  // New line
  fprintf(f_dump, "\n");

  return;
}

void WriteMeanStdDevToStdOut_lorenz(const int step, const double tym, const double x, const double y, const double z)
{
  // write time, x/y/z to screen 
  cout << "Time Step: " << step << ", time: " << tym;
  cout << ", x: " << x << ", y: " << y << ", z:" << z;
  cout << endl;

  return;
}

