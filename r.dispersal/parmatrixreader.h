/* parmatrixreader.h -- read the parameter file

   Copyright (C) 2001 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  */

/* Written by Andreas Gros */


#include <stdlib.h>
#include <stdio.h> 

#ifndef DEFINITIONS
#include "definitions.h"
#endif

/* Function to get the parameters out of the parameter file into a temporary matrix.
 * Parameters: matrix to be filled, 
 *             fp - open file pointer to the parameter file
 */
parameters read_parameter_file_contents(FILE *fp) {
  int i, j;
  parameters parmatrix;
  int localCols, localRows;
  const int MAXLINELENGTH = 1024;
  char *dummy;

  char seps[]   = " ,\t\n";
  /*get the first line with the number of rows and cols*/
  if (fscanf(fp, "%i %i", &localRows, &localCols) != 2) {
    G_fatal_error(_("The first line of the parameter file has to consist of the number of rows and cols of the file"));
  }

  parmatrix.cols = localCols;
  parmatrix.rows = localRows;
  fprintf(stderr, "here 0");
  /*read the rest of the line (if there is any)*/
  fgets(dummy, MAXLINELENGTH,  fp);
  /*allocate enough space for everything*/
  parmatrix.entrymatrix = (double**)malloc(localRows * sizeof(double*));
  
  for (i = 0; i < localRows; i++) {
    parmatrix.entrymatrix[i] = (double*)malloc(localCols * sizeof(double));
    for (j = 0; j < localCols; j++) {
      if(fscanf(fp, "%lf",  &(parmatrix.entrymatrix[i][j]))!=1) {
        G_fatal_error(_("The conversion of the %d-th number in line %d of the parameter file failed. Bailing out!"), j, i);
      }
      else {
        fprintf(stderr, "conversion result: %lf \n ", parmatrix.entrymatrix[i][j]);
      }
    }
  }

  /*output the contents of the parameter file*/
  for (i = 0; i < localRows; i++) {
    for (j = 0; j < localCols; j++) {
      fprintf(stderr, "%lf ",parmatrix.entrymatrix[i][j]);  
    }
    fprintf(stderr, "::\n");  
  }
  return parmatrix;
}

double* fillParameterVector(const int INDEX, parameters parmatrix) {
  int i;
  if(parmatrix.rows + INDEX < parmatrix.cols) {  
    double *slv = (double*)malloc(parmatrix.rows * sizeof(double));
    for (i = 0; i < parmatrix.rows; i++) {
      slv[i] = parmatrix.entrymatrix[i][parmatrix.rows+INDEX];
    }
    return slv;
  } else {
    fprintf(stderr, "in parameterfile: ratio rows to cols does not fit: cols must be at least rows+%d: rows are %d, cols are %d ", INDEX, parmatrix.rows, parmatrix.cols);  
    return NULL;
  }
}

double** fillTransitionMatrix(parameters parmatrix) {
  const int nrows = parmatrix.rows;
  double **transmat = (double**)malloc(nrows * sizeof(double*));
  int i, j;
  for(i = 0; i < nrows; i++) {
    transmat[i] = (double*)malloc(nrows * sizeof(double));
    for(j = 0; j < nrows; j++) {
      transmat[i][j] = parmatrix.entrymatrix[i][j];
    }
  }
  return transmat;  
}

double* fillSteplengthVector(const int STEPLENGTHINDEX, parameters parmatrix) {
  return fillParameterVector(STEPLENGTHINDEX, parmatrix);
}

//double* fillFrictionVector(const int FRICTIONINDEX, parameters parmatrix)
//  {
//  return fillParameterVector(FRICTIONINDEX, parmatrix);
//  }

double* fillEnergycostVector(const int ENERGYINDEX, parameters parmatrix) {
  return fillParameterVector(ENERGYINDEX, parmatrix);
}

double* fillHabitattypeVector(const int HABITATINDEX, parameters parmatrix) {
  return fillParameterVector(HABITATINDEX, parmatrix);
}
