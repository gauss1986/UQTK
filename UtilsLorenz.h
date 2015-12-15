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
#include "Array1D.h"

/// \brief Write mean and std. dev. to a file
///
/// Arguments:
///   \li const double tym: current time
///   \li const double u0: means of solution components u
///   \li const double u_std: std. dev. of three solution components u
///   \li FILE* f_dump: C file pointer to write to
void WriteMeanStdDevToFilePtr_lorenz(const double tym, const double x, const double y, const double z,  FILE* f_dump);

/// \brief Write x/y/z to screen
///
/// Arguments:
///   \li const int step: current time step
///   \li const double x/y/z: value of x/y/z in the Lorenz model
void WriteMeanStdDevToStdOut_lorenz(const int step, const double tym, const double u, const double y, const double z);
