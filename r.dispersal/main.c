/****************************************************************************
 *
 * MODULE:       r.dispersal 
 * AUTHOR(S):    Andreas Gros - andigros72@gmail.com
 * PURPOSE:      To provide an assessment tool for single habitat patches 
 *               settled by a metapopulation. Output is the number of immigrants
 *               per patch.
 *               Required input: A map with different landscape types, a map
 *               of starting points for the migrants and a map with the replacement 
 *               habitat types in case patches are taken out for assessment.
 * COPYRIGHT:    (C) 2013 by the GRASS Development Team
 *
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with GRASS
 *               for details.
 *
 * Version       0.90 
 *                - relative number of immigrants (the number of immigrants 
 *                  from each patch in relation to all immigrants)  
 *                - number of successfull emigrants (count emigrants from each
 *                  patch)    
 *****************************************************************************/

#define MAIN

#include <stdlib.h>
#include <stdio.h> //for reading the parameter file
#include <float.h> //for DBL_MAX
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <fcntl.h>
#include <grass/gis.h>
#include <grass/site.h>
#include <grass/segment.h>
#include <grass/glocale.h>
#include "local_proto.h"
#include "definitions.h"
#include "inttostr.h"

#include "parmatrixreader.h"

struct Cell_head window;

/***************************************************************************************
** find parameter "currtype" in a given parameter vector "searchvec"
** and return the value at the same position of parameter vector
** "resvec"
 ***************************************************************************************/
double findInVec(double *resvec, double *searchvec, double currtype, int searchveclength) {
  int i;
  double result = 0;
  for (i = 0; i < searchveclength; i++) {
    if (currtype == searchvec[i]) {
      result = resvec[i];
      break;
    }
  }
  return result;
}

/***************************************************************************************
** find parameter value "tofind" in parameter vector "vec" up to a
** given parameter "length" and return its position if found, -1 if it
** was not found
***************************************************************************************/
int findInIntVec(int tofind, int *vec, int length) {
  int i = 0;
  if ((tofind < length) && (vec[tofind] == tofind)) i = tofind;
  else while((i < length) && (vec[i] != tofind)) i++;
  return (i == length) ? NONVALID : i;
}

/******************************************************************************
** Allocate space for the habitat matrix and initialize it with NONHABITAT
******************************************************************************/
habitat InitializeMatrixOfHabitat(habitat h) {
  int i, j;
  const int yentr = (h.extension.ymax - h.extension.ymin + 1);
  const int xentr = (h.extension.xmax - h.extension.xmin + 1);

  h.entrymatrix = (int**)malloc(yentr * sizeof(int*));
  for (i = 0; i < yentr; i++) {
    h.entrymatrix[i] = (int*)malloc(xentr * sizeof(int));
    for (j = 0; j < xentr; j++) {
      h.entrymatrix[i][j] = NONHABITAT;
    }
  }
  return h;
}

/******************************************************************************
** Fill the raster buffer from the landscape input file
******************************************************************************/
int fillRasterBuffer(
    int row, 
    int minrow, 
    int predictedMaximumNumberOfSteps, 
    void **buffervec, 
    int in_landscape_fd, 
    RASTER_MAP_TYPE landscape_data_type) {
  int i, focusIndex = 0;
  //fill the raster buffer
  for (i = 0; i < predictedMaximumNumberOfSteps; i++) {
    Rast_get_row(in_landscape_fd, buffervec[i], (minrow + i), landscape_data_type);
    if (minrow + i == row) focusIndex = i;    
  }
  return focusIndex;
}

/******************************************************************************
** Calculate the current buffer size
******************************************************************************/
int assessCurrentBufferMemsize(ncols, predictedMaximumNumberOfSteps) {
  return (ncols * predictedMaximumNumberOfSteps); 
}



/******************************************************************************
** Simulate the movement of individuals released from the habitat cell
** at x, y (parameters)
******************************************************************************/
void searchPath(
    int in_landscape_fd, 
    int in_replacement_fd, 
    int from_x, 
    int to_x, 
    int y, 
    int nrows, 
    int ncols, 
    RASTER_MAP_TYPE landscape_data_type, 
    RASTER_MAP_TYPE replacement_data_type, 
    double cell_ext, 
    int release_num_ind, 
    double **transitionProbabilities, 
    int numhabitattypes, 
    double *steplengths, 
    double *habitattypes, 
    double *energycosts, 
    double initial_energy, 
    double perceptual_range, 
    int *patchhash, 
    int numStartPatches, 
    int missingPatchIndex, 
    habitat missinghabitat, 
    int **resmap ) {

  int bufferWasShifted = FALSE;
  int allowedToEnlargeBuffer = TRUE;
  int tdir, direction;         // one of eight neighbouring cells
  int i, j, numi, stepcount;   
  int row = y;                 // current row
  int col = from_x;            // current col
  int tcol, trow;              // temporal indices for column and row
  int aim_index;
  void  *tbuf,                 //a temporary buffer for the current line, 
                               //the moving individual is in
        *replacement_buf;      //if it hits a "missing" habitat, the 
                               //replacement habitat is opened   
  int replacement_buf_row = 0;
  register individual ind;     // the individual
  int initial_habitat_type, tht;
  register double initial_steplength, tsl;
  register double initial_energycosts, tec;
  register int setMaxRow = 0, setMinRow = 0;
  CELL surrounding[SEARCHMASK_WIDTH]; //surroundingmask_WIDTH 
                                      //defined in definitions.h 
  double transitionprobs_surrounding[SEARCHMASK_WIDTH];

  int headed_north, headed_south, headed_west, headed_east;
  const int CELLDIST_HEADED = 10;

  //find the habitat type with the smallest energy loss:
  double lowest_energy_cost = DBL_MAX;
  for (i = 0; i < numhabitattypes; i++) {
    if (energycosts[i] < lowest_energy_cost) {
      lowest_energy_cost = energycosts[i];
    }
  }

  //calculate the maximum number of line buffers for individual movement
  //can go straight up or down (-> *2)
  int predictedMaximumNumberOfSteps = 
    (int)(ceil(initial_energy / lowest_energy_cost) * 2 + 1);

  if (ROW_ALLOC_MAX < predictedMaximumNumberOfSteps) {  
    predictedMaximumNumberOfSteps = ROW_ALLOC_MAX;
  }
  
  DEBUG fprintf(
    stderr, 
    "the predicted maximum number of steps is [%d]\n", 
    predictedMaximumNumberOfSteps);
  
  //allocate space for the line buffers
  void **buffervec = (void**)malloc(
    predictedMaximumNumberOfSteps * sizeof(void*)
  );
  void **repl_buffervec;

  for (i = 0; i < predictedMaximumNumberOfSteps; i++) {
    buffervec[i] = Rast_allocate_c_buf();//  G_allocate_c_raster_buf();
  }

  //allocate space for the cell-data to be read
  //center_buf      = G_allocate_c_raster_buf();
  //lower_buf       = G_allocate_c_raster_buf();
  //upper_buf       = G_allocate_c_raster_buf();
  //replacement_buf = G_allocate_c_raster_buf();
  replacement_buf = Rast_allocate_c_buf();

  //calculate the center element of the raster buffer
  int pivot = (int)floor(predictedMaximumNumberOfSteps/2);

  //the lowest possible resp. neccessary row
  int minrow = (row - pivot) > 0 ? row - pivot : 0;
  DEBUG fprintf(stderr, "Minrow is [%d]\n", minrow);

  //the maximum possible resp. neccessar row
  int maxrow = (minrow + predictedMaximumNumberOfSteps) < nrows ? 
    minrow + predictedMaximumNumberOfSteps : 
    nrows - 1;

  DEBUG fprintf(stderr, "Maxrow is [%d]\n", maxrow);
  
  //now we can calculate the possible maximumNumberOfSteps (straight upwards 
  //or downwards); it can be that we hit a boundary 
  predictedMaximumNumberOfSteps = (maxrow - minrow); 

  DEBUG fprintf(
    stderr, 
    "the actual predicted maximum number of steps is [%d]\n", 
    predictedMaximumNumberOfSteps
  );

  //fill the buffer and get the index in the vector of buffered rows pointing 
  //to the focus row
  int focusIndex = fillRasterBuffer(
    row, 
    minrow, 
    predictedMaximumNumberOfSteps, 
    buffervec, 
    in_landscape_fd,
    landscape_data_type
  );

  setMaxRow = maxrow;
  setMinRow = minrow;

  int currentBufferMemsize = 
    assessCurrentBufferMemsize(ncols, predictedMaximumNumberOfSteps);
  if (currentBufferMemsize > MAX_MEM_SIZE_USED_BY_BUFFER) {
    allowedToEnlargeBuffer = FALSE;
  }

  initial_habitat_type = ((CELL *) buffervec[focusIndex])[col];
  if ((Rast_is_c_null_value(&initial_habitat_type) == 1) || 
      (initial_habitat_type == 0)) {
    initial_habitat_type = NONHABITAT;
  } else {
    //find initial steplength
    initial_steplength = findInVec(
      steplengths, 
      habitattypes, 
      initial_habitat_type, 
      numhabitattypes
    );
    initial_energycosts = findInVec(
      energycosts, 
      habitattypes, 
      initial_habitat_type, 
      numhabitattypes
    );
    if ((initial_energycosts == 0) || (initial_steplength == 0)) {
      G_fatal_error(
        _("Faulty parameter habitattype [%d]: found no energycosts in parameterfile [%lf], steplength is [%lf]"), 
        initial_habitat_type, 
        initial_energycosts, 
        initial_steplength
      );
    }

    //check if there is a habitat that should be left out of the picture:
    if (missingPatchIndex >= 0) {
      missingPatchIndex = patchhash[missingPatchIndex]; 
    }

    const int orig_row = y;
    int orig_col = 0;
    COUNTMIGRANTS {
      headed_north = 0;
      headed_south = 0;
      headed_west  = 0;
      headed_east  = 0;
    }

    //for each cell in the focus row from which individuals are supposed 
    //to be released
    int release_col; 
    int icheck;
    for (release_col = from_x; release_col <= to_x; release_col++) {
      DEBUG fprintf(stderr, "operating on col [%d]\n", release_col);

      //release the individuals
      for (numi = 0; numi < release_num_ind; numi++) {
        icheck = FALSE;

        //if neccessary, but back the old raster buffer
        if ( bufferWasShifted == TRUE ) {
          minrow = setMinRow;
          maxrow = setMaxRow;
          focusIndex = fillRasterBuffer(
            y, 
            minrow, 
            predictedMaximumNumberOfSteps, 
            buffervec, 
            in_landscape_fd,
            landscape_data_type
          );
          bufferWasShifted = FALSE; 
        }

        row = y;
        col = release_col;
        COUNTMIGRANTS orig_col = release_col;
        ind.energylevel = initial_energy; //charge an individual  
        ind.x = release_col; // now resolution only at the cell-level
        ind.y = y; // to be used later
        direction = (int)floor(frand()*SURROUNDINGMASK_WIDTH);  
        //(direction \in {0;7} )
        stepcount = 0;                                                      
        tht = initial_habitat_type;
        tsl = initial_steplength;
        tec = initial_energycosts;
        register double totalTransprob = 0; 
        register int c_index = focusIndex; //the current index in buffervec
        DEBUG fprintf(stderr, "first c_index [%d]\n", c_index);
        do { //handle a single individual till its energy runs out 
          totalTransprob = 0; 
          //collect information about the surrounding
          for (i = LEFTRANGE; i <= RIGHTRANGE; i++) {
            tdir = direction + i;
            if (tdir < 0) {
              tdir += SURROUNDINGMASK_WIDTH;
            } else if (tdir >= SURROUNDINGMASK_WIDTH) {
              tdir -= SURROUNDINGMASK_WIDTH;
            }
            tcol = col + surroundingmask[tdir][XPOS];
            trow = row + surroundingmask[tdir][YPOS];
            tbuf = NULL;
            if (surroundingmask[tdir][YPOS] == -1) { //downwards
              if (c_index > 0) {
                tbuf = buffervec[c_index - 1]; 
              }
            } else if (surroundingmask[tdir][YPOS] == 1) { //or upwards
              if (c_index < predictedMaximumNumberOfSteps-1) {
                tbuf = buffervec[c_index + 1]; 
              }
            } else { 
              tbuf = buffervec[c_index]; //left or right
            }

            //check if the habitat type is in range of the parameter 
            //vector "habitattypes"
            //if not, it's classified as nonhabitat
            if ((tbuf != NULL) && (((CELL *)tbuf)[tcol] < numhabitattypes) && 
                (Rast_is_c_null_value(&((CELL *)tbuf)[tcol]) != 1)) {
              aim_index = ((CELL *)tbuf)[tcol];
              if ((missingPatchIndex >= 0) && 
                  (tcol >= missinghabitat.extension.xmin) && //check if cell in 
                  (tcol <= missinghabitat.extension.xmax) && //range of the 
                  (trow >= missinghabitat.extension.ymin) && //habitat to be 
                  (trow <= missinghabitat.extension.ymax) && //left out of the
                  (missinghabitat.entrymatrix[               //picture
                    trow - missinghabitat.extension.ymin
                  ][
                    tcol - missinghabitat.extension.xmin
                  ] == missingPatchIndex)) {            //if it is, take the 
                                                        //replacement quality
                Rast_get_row(in_replacement_fd, replacement_buf, trow, replacement_data_type); 
                replacement_buf_row = trow;  
                aim_index = 
                  (Rast_is_c_null_value(&((CELL *)replacement_buf)[tcol]) == 1) ? 
                    NONHABITAT : 
                    ((CELL *)replacement_buf)[tcol];
              }
            }
            else {
              aim_index = NONHABITAT;
            }
            //store the index of the habitat type the 
            //individual is looking into
            surrounding[i - LEFTRANGE] = aim_index; 
            
            //store the transition probabilities
            transitionprobs_surrounding[i - LEFTRANGE] = 
              transitionProbabilities[aim_index][tht]; 

            //and sum up the transition probabilities
            totalTransprob += transitionprobs_surrounding[i - LEFTRANGE]; 
          }
          const double threshold = frand()*totalTransprob; //now pick one

          i = -1;
          double tp = 0;
          do {
            i++; 
            tp += transitionprobs_surrounding[i]; 
          } while ((tp < threshold) && (i < SEARCHMASK_WIDTH-1));
          //chosen cell is direction + i;
          tht = surrounding[i]; //new habitat type

          ind.energylevel -= findInVec(
            energycosts, 
            habitattypes, 
            tht, 
            numhabitattypes
          );
          tdir = direction + i + LEFTRANGE;
          if (tdir < 0) { 
            tdir += SURROUNDINGMASK_WIDTH; 
          } else if (tdir >= SURROUNDINGMASK_WIDTH) { 
            tdir -= SURROUNDINGMASK_WIDTH; 
          }
          //now check to where it moved:
          row += surroundingmask[tdir][YPOS];
          c_index += surroundingmask[tdir][YPOS];
          col += surroundingmask[tdir][XPOS];

          COUNTMIGRANTS {
            if ((icheck == FALSE)) {
              //check in which direction it went
              if ((row - orig_row) > CELLDIST_HEADED) {
                headed_north++;
                icheck = TRUE;
              } else if ((row - orig_row) < -CELLDIST_HEADED) {
                headed_south++;
                icheck = TRUE;
              }
              if ((col - orig_col) > CELLDIST_HEADED) {
                headed_east++;
                icheck = TRUE;
              } else if ((col - orig_col) < -CELLDIST_HEADED) {
                headed_west++;
                icheck = TRUE;
              }
            }
          }
          if ((row <= 0) || (row >= nrows - 1) || 
              (col <= 0) || (col >= ncols - 1)) { 
            ind.energylevel = 0;
          }

          DEBUG fprintf(
            stderr, 
            "calculate the next step c_index [%d]\n", 
            c_index
          );

          if (ind.energylevel > 0) {
            if (surroundingmask[tdir][YPOS] < 0) {
              //check if individual runs out of predicted range:
              if ((c_index == 0) && (row > 0)) {
                DEBUG fprintf(stderr, "hitting zero, still going south\n");

                //if there's still memory space (MAX_MEM_SIZE_USED_BY_BUFFER) 
                //left, try to enlarge buffervec
                if (allowedToEnlargeBuffer == TRUE) { 
                  DEBUG fprintf(stderr, "try to enlarge\n");

                  int increment = (minrow - ROW_ALLOC_INCREMENT > 0) ? 
                    ROW_ALLOC_INCREMENT : 
                    minrow; 
                  minrow -= increment;
                  setMinRow = minrow;
                  repl_buffervec = 
                    (void**)malloc((maxrow - minrow) * sizeof(void*));
                  
                  //allocate the missing buffers
                  for (i = 0; i < increment; i++) {
                    
                    repl_buffervec[i] = Rast_allocate_c_buf();//G_allocate_c_raster_buf();
                    Rast_get_row(in_landscape_fd, repl_buffervec[i], minrow + i, landscape_data_type); 
                  }
                  for (i = 0; i < predictedMaximumNumberOfSteps; i++) {
                    repl_buffervec[i+increment] = buffervec[i];
                  }
                  predictedMaximumNumberOfSteps = maxrow - minrow;
                  DEBUG fprintf(
                    stderr, 
                    _("predictedMaximumNumberOfSteps is now [%d]\n"),
                    predictedMaximumNumberOfSteps  
                  );

                  //adapt focusIndex to the enlarged buffervec:
                  focusIndex = focusIndex + increment;
                  c_index = increment;

                  //now swap the buffer pointer
                  buffervec = repl_buffervec;
                  repl_buffervec = NULL;
                  currentBufferMemsize = assessCurrentBufferMemsize(
                    ncols, 
                    predictedMaximumNumberOfSteps
                  );
                  if (currentBufferMemsize > MAX_MEM_SIZE_USED_BY_BUFFER) {
                    //from then on only shifting is allowed
                    allowedToEnlargeBuffer = FALSE;  
                    DEBUG fprintf(
                      stderr, 
                      "Maximum buffersize reached, switch to shifting buffer\n"
                    );
                  }
                  DEBUG fprintf(
                    stderr, 
                    _("Enlarged the row buffer at the bottom! maxrow: [%d]; minrow: [%d]\n"), 
                    maxrow, 
                    minrow
                  );
                } else { //shift the buffervector
                  int increment = (minrow - ROW_ALLOC_INCREMENT > 0) ? 
                    ROW_ALLOC_INCREMENT : 
                    minrow; 
                  minrow -= increment;
                  maxrow -= increment;
                  c_index += increment;
                  focusIndex = fillRasterBuffer(
                    y, 
                    minrow, 
                    predictedMaximumNumberOfSteps, 
                    buffervec, 
                    in_landscape_fd,
                    landscape_data_type
                  );
                  //attention: focusIndex could be 0 afterwards
                  bufferWasShifted = TRUE; 
                }
              } //end if ((c_index == 0) && (row > 0)) 
            } //end if (surroundingmask[tdir][YPOS] == -1)
            else if (surroundingmask[tdir][YPOS] > 0 ) {
              //check if individual runs out of predicted range:
              if ((c_index == predictedMaximumNumberOfSteps) && 
                  (maxrow < nrows - 1) && 
                  (row < nrows - 1)) { //no higher buffer possible
                //try to enlarge buffervec
                if (allowedToEnlargeBuffer == TRUE) { 
                  DEBUG fprintf(stderr, "going north, try to enlarge\n");

                  int increment = (maxrow + ROW_ALLOC_INCREMENT < nrows) ? 
                    ROW_ALLOC_INCREMENT : 
                    nrows - 1 - maxrow; 
                  maxrow += increment;
                  setMaxRow = maxrow;
                  repl_buffervec = 
                    (void**)malloc((maxrow - minrow) * sizeof(void*));
                  for (i = 0; i < predictedMaximumNumberOfSteps; i++) {
                    repl_buffervec[i] = buffervec[i];
                  }

                  //allocate the missing buffers
                  for (i = predictedMaximumNumberOfSteps; 
                      i < predictedMaximumNumberOfSteps + increment; 
                      i++) {
                    repl_buffervec[i] = Rast_allocate_c_buf(); //G_allocate_c_raster_buf();
                    Rast_get_row(in_landscape_fd, repl_buffervec[i], minrow + i, landscape_data_type);
                  }

                  predictedMaximumNumberOfSteps = maxrow - minrow;
                  
                  DEBUG fprintf(
                    stderr, 
                    _("predictedMaximumNumberOfSteps is now [%d]\n"),
                    predictedMaximumNumberOfSteps  
                  );

                  //now swap the buffer pointer
                  buffervec = repl_buffervec;
                  repl_buffervec = NULL;
                  currentBufferMemsize = assessCurrentBufferMemsize(
                    ncols, predictedMaximumNumberOfSteps
                  );
                  if (currentBufferMemsize > MAX_MEM_SIZE_USED_BY_BUFFER) {
                    //from then on only shifting is allowed
                    allowedToEnlargeBuffer = FALSE; 
                    DEBUG fprintf(
                      stderr, 
                      _("Maximum buffersize reached, switch to shifting buffer\n")
                    );
                  }
                  DEBUG fprintf(
                    stderr, 
                    _("Enlarged the row buffer at the top! maxrow: [%d]; minrow: [%d]\n"), 
                    maxrow, 
                    minrow
                  );
                }
                else {//resort to shifting the buffer
                  int increment = (maxrow + ROW_ALLOC_INCREMENT < nrows) ? 
                    ROW_ALLOC_INCREMENT : 
                    nrows - 1 - maxrow; 
                  maxrow += increment;
                  minrow += increment;
                  c_index -= increment;
                  focusIndex = fillRasterBuffer(
                    y, 
                    minrow, 
                    predictedMaximumNumberOfSteps, 
                    buffervec, 
                    in_landscape_fd,
                    landscape_data_type
                  );
                  //attention: focusIndex could be 0 afterwards
                  bufferWasShifted = TRUE; 
                }
              }//end if ((c_index == predictedMaximumNumberOfSteps-1) && 
               //(row < nrows-1)) 
            }
          }
          stepcount++; 
          DEBUG fprintf(
            stderr, 
            _("stepcount [%d], row [%d], col [%d]\n"),
            stepcount, 
            row, 
            col 
          );
        } while ((ind.energylevel > 0) && (stepcount < MAXNUMSTEPS));
        //now count it:
        if ((row >= 0) && (row < nrows) && (col > 0) && (col < ncols) ) {
          resmap[row][col]++;
        }
      }//end for (i = 0; i < release_num_ind; i++)
    }//end for (release_col = from_x; release_col <= to_x; release_col++)
  }

  for (i=0; i < predictedMaximumNumberOfSteps; i++) {
    if (buffervec[i]) { 
      G_free(buffervec[i]);
    }
  }
  if (replacement_buf) { 
    G_free(replacement_buf);
  }
  COUNTMIGRANTS fprintf(
    stderr, 
    _("headed north [%d], south [%d], west [%d], east [%d] \n"),
    headed_north, 
    headed_south, 
    headed_west, 
    headed_east 
  );

}

/***************************************************************************************
** Main action: Go through the habitat map, find and measure all
** habitats, call searchPath for each habitat cell and register, where
** the released individuals died. Then do the same repeatedly for 
** sequentially taking out one patch and measure the difference to the
** scenario with no patch missing. The "missing" patches are replaced
** by the corresponding entry in "replacement_map".
***************************************************************************************/

void doTheModeling(
    char landscape_layer[64], 
    char replacement_layer[64], 
    char startpatches_layer[64], 
    parameters parmatrix, 
    double *steplengthvec, 
    double *frictionvec, 
    double *energycostsvec,
    double *habitattypesvec,
    int numberOfIndividualsReleasedPerCell,
    double initial_energy,
    double perceptual_range,
    char* outputpath,
    int flag_assessmentmodePAMM,
    int flag_assessmentmodePAMI) {

  char *landscape_fname,
       *replacement_fname,
       *startpatches_fname,
       *startpatches_cat_fname,
       *search_mapset,
       *landscape_mapset,
       *replacement_mapset,
       *startpatches_mapset,
       *startpatches_cat_mapset;
  int i, j;
  int numStartPatches, startPatchCount;
  const int BUFSIZE = 1024;
  int in_startpatches_fd,  
      in_landscape_fd, 
      in_replacement_fd, 
      out_result_fd;
  RASTER_MAP_TYPE startpatches_data_type, 
                  landscape_data_type, 
                  replacement_data_type;
  void *startpatches_cell_buf, 
       *startpatches_res_cell_buf,
       *landscape_cell_buf, 
       *replacement_cell_buf, 
       *output_cell_buf;
  struct Cell_head landscape_cellhd, 
                   replacement_cellhd,
                   startpatches_cellhd;
  int landscape_head_ok, 
      replacement_head_ok, 
      startpatches_head_ok, 
      startpatches_cat_head_ok;

  double **transitionMatrix;

  int **finalresult = NULL, **cellresult = NULL;
  int **patchresult = NULL, **resmap = NULL;

  int nrows, ncols, row, col;
  register int source, destination;
  int missingPatchIndex;

  fprintf(stderr, "starting\n");
  //nrows = G_window_rows();
  //ncols = G_window_cols();
  nrows = Rast_window_rows();
  ncols = Rast_window_cols();
  fprintf(stderr, "number of rows: %d\n", parmatrix.rows);
  
  if ((steplengthvec = 
      fillSteplengthVector(STEPLENGTHINDEX, parmatrix)) == NULL) {
    G_fatal_error(_("cannot fill steplength vector!")); 
  }
  if ((energycostsvec = 
      fillEnergycostVector(ENERGYINDEX, parmatrix)) == NULL) {
    G_fatal_error(_("cannot fill energycost vector!")); 
  }
  if ((habitattypesvec = 
      fillHabitattypeVector(HABITATINDEX, parmatrix)) == NULL) {
    G_fatal_error(_("cannot fill habitattype vector!")); 
  }
  if ((transitionMatrix = fillTransitionMatrix(parmatrix)) == NULL) {
    G_fatal_error(_("cannot fill transitionMatrix!")); 
  }

  //print out steplengthvec
  fprintf(stderr, "\n\nsteplengthvec: \n");
  for (i = 0; i < parmatrix.rows; i++) {
    fprintf(stderr, "%lf \n", steplengthvec[i]);
  }  
  fprintf(stderr, "\n\nenergycostsvec: \n");
  for (i = 0; i < parmatrix.rows; i++) {
    fprintf(stderr, "%lf\n", energycostsvec[i]);
  }
  fprintf(stderr, "\n\nhabitattypesvec: \n");
  for (i = 0; i < parmatrix.rows; i++) {
    fprintf(stderr, "%lf\n", habitattypesvec[i]);
  }
 
  //number_of_rows = parameters_cellhd.rows;
  //number_of_cols = parameters_cellhd.cols;
  //nrows = G_window_rows();
  //ncols = G_window_cols();

  landscape_mapset = G_find_raster(landscape_layer, "");
  startpatches_mapset = G_find_raster(startpatches_layer, "");
  replacement_mapset = G_find_raster(replacement_layer, "");

  if (landscape_mapset == NULL) {
    G_fatal_error(_("%s - not found"), landscape_layer);
  }
  if (startpatches_mapset == NULL) {
    G_fatal_error(_("%s - not found"), startpatches_layer);
  }
  if (replacement_mapset == NULL) {
    G_fatal_error(_("%s - not found"), replacement_layer);
  }
 
  startpatches_fname  = startpatches_layer;
  landscape_fname = landscape_layer;
  replacement_fname = replacement_layer;
 
  //count the number of different start patches
  struct Range range;
  int zmin, zmax;
  
  //if (G_read_range(startpatches_fname, startpatches_mapset, &range) < 0) {
  if (Rast_read_range(startpatches_fname, startpatches_mapset, &range) < 0) {
    G_fatal_error(_("could not read range file"));
  }

  //G_get_range_min_max(&range, &zmin, &zmax);
  //G_get_range_min_max(&range, &zmin, &zmax);

  Rast_get_range_min_max(&range, &zmin, &zmax);
  Rast_get_range_min_max(&range, &zmin, &zmax);

  numStartPatches = zmax - zmin + 1;
  if (numStartPatches <= 0) {
    numStartPatches = INITIAL_NUMBER_OF_STARTPATCHES;
  }
  fprintf(
    stderr, 
    _("the range of startpatches is from %d to %d, which gives us a number of %d startpatches."), 
    zmin, 
    zmax, 
    numStartPatches
  );

  if ((in_startpatches_fd = 
      //G_open_cell_old(startpatches_fname, startpatches_mapset)) < 0) {
      Rast_open_old(startpatches_fname, startpatches_mapset)) < 0) {

    G_fatal_error(
      _("Cannot open cell file [%s], in mapset [%s]"), 
      startpatches_fname, 
      startpatches_mapset
    );
  }
  if ((in_landscape_fd = 
      //G_open_cell_old(landscape_fname, landscape_mapset)) < 0) {
      Rast_open_old(landscape_fname, landscape_mapset)) < 0) {
    G_fatal_error(
      _("Cannot open cell file [%s], in mapset [%s]"), 
      landscape_fname, 
      landscape_mapset
    );
  }
  if ((in_replacement_fd = 
      //G_open_cell_old(replacement_fname, replacement_mapset)) < 0) {
      Rast_open_old(replacement_fname, replacement_mapset)) < 0) {
    G_fatal_error(
      _("Cannot open cell file [%s], in mapset [%s]"), 
      replacement_fname, 
      replacement_mapset
    );
  }
  Rast_get_cellhd(landscape_layer, landscape_mapset, &landscape_cellhd);
  Rast_get_cellhd(startpatches_layer, startpatches_mapset, &startpatches_cellhd);
  Rast_get_cellhd(replacement_layer, replacement_mapset, &replacement_cellhd);

  /* determine the inputmap type (CELL/FCELL/DCELL) */
  //landscape_data_type = G_raster_map_type(landscape_layer, "");
  //startpatches_data_type = G_raster_map_type(startpatches_layer, "");
  //replacement_data_type = G_raster_map_type(replacement_layer, "");
  
  landscape_data_type = Rast_get_map_type(in_landscape_fd);
  startpatches_data_type = Rast_get_map_type(in_startpatches_fd);
  replacement_data_type = Rast_get_map_type(in_replacement_fd);

  fprintf(stderr, "landscape_mapset: [%s]\n", landscape_mapset);
  fprintf(stderr, "startpatches_mapset: [%s]\n", startpatches_mapset);

  //startpatches_cell_buf = G_allocate_raster_buf(startpatches_data_type);
  //startpatches_res_cell_buf = G_allocate_raster_buf(startpatches_data_type);
  //replacement_cell_buf = G_allocate_raster_buf(replacement_data_type);

  startpatches_cell_buf = Rast_allocate_buf(startpatches_data_type);
  startpatches_res_cell_buf = Rast_allocate_buf(startpatches_data_type);
  replacement_cell_buf = Rast_allocate_buf(replacement_data_type);
  //make space for habitats:
  habitat *starthabitatvec = (habitat*)malloc(numStartPatches*sizeof(habitat));

  //initialize them
  for (i = 0; i < numStartPatches; i++) {
    starthabitatvec[i].extension.xmin = ncols;
    starthabitatvec[i].extension.xmax = 0;
    starthabitatvec[i].extension.ymin = nrows;
    starthabitatvec[i].extension.ymax = 0;
  }

  //the patchhash is neccessary for cases where patches are not numbered
  //sequentially or in case some patch numbers are missing
  //it stores the category number of the startpatches
  int *patchhash = (int*)malloc(numStartPatches * sizeof(int)); 
  startPatchCount = 0; //count the actual number of start patches
  for (row = 0; row < nrows; row++) {
    CELL c;

    /* read input map and copy it into the tmp_matrix*/
    //if (G_get_c_raster_row(in_startpatches_fd, startpatches_cell_buf, row) < 0) {
    Rast_get_row(in_startpatches_fd, startpatches_cell_buf, row, startpatches_data_type);

    for (col = 0; col < ncols; col++) {
      switch (startpatches_data_type) {
        case CELL_TYPE:
          c = ((CELL *) startpatches_cell_buf)[col];
          if ((c != NONHABITAT) && (Rast_is_c_null_value(&c) != 1)) {
            if (findInIntVec(c, patchhash, startPatchCount) < 0) {
              patchhash[startPatchCount] = c;
              startPatchCount++; //a new startpatch was found
            }
            i = findInIntVec(c, patchhash, startPatchCount);
            if (i >= 0) {
              if (starthabitatvec[i].extension.xmin > col) {
                starthabitatvec[i].extension.xmin = col;
              }  
              if (starthabitatvec[i].extension.xmax < col) { 
                starthabitatvec[i].extension.xmax = col;
              }
              if (starthabitatvec[i].extension.ymin > row) { 
                starthabitatvec[i].extension.ymin = row;
              }
              if (starthabitatvec[i].extension.ymax < row) { 
                starthabitatvec[i].extension.ymax = row;
              }
              starthabitatvec[i].extension.habitattype = c;
            }
          }
          break;
        default: G_fatal_error(
          _("Habitat patch numbers in <%s> are to be integers!"), 
          startpatches_fname
        );
      }
    }  
  }

  numStartPatches = startPatchCount;
  fprintf(stderr, "the number of startpatches is: %d \n", numStartPatches);

  //organize Space and fill them:
  for (i = 0; i < numStartPatches; i++) {
    starthabitatvec[i] = InitializeMatrixOfHabitat(starthabitatvec[i]);
    fprintf(stderr, "habitattype %d, xmin %d xmax %d, ymin %d ymax %d\n", 
      starthabitatvec[i].extension.habitattype,
      starthabitatvec[i].extension.xmin,
      starthabitatvec[i].extension.xmax,
      starthabitatvec[i].extension.ymin,
      starthabitatvec[i].extension.ymax
    );
  }

  //the filling:
  for (row = 0; row < nrows; row++) {
    CELL c;
    Rast_get_row(in_startpatches_fd, startpatches_cell_buf, row, startpatches_data_type);

    for (col = 0; col < ncols; col++) {
      c = ((CELL *) startpatches_cell_buf)[col];
      if ((c != NONHABITAT) && (Rast_is_c_null_value(&c) != 1)) {
        i = findInIntVec(c, patchhash, numStartPatches);
        if (i >= 0) {
          starthabitatvec[i].entrymatrix[
            row - starthabitatvec[i].extension.ymin
          ][
            col - starthabitatvec[i].extension.xmin
          ] = c;
          starthabitatvec[i].numcells++;
        }
      }
    }  
  }

  int *numReleasedPerPatch = (int*)malloc(numStartPatches * sizeof(int));
  for (i = 0; i < numStartPatches; i++) {
    numReleasedPerPatch[i] = 
      starthabitatvec[i].numcells * numberOfIndividualsReleasedPerCell;
  }
  //the search function for the primitive hash is:
  //findInIntVec(int tofind, int *vec, int length)


  //create first assessment matrix
  //it consists of n-startpatch columns and (n+1)-startpatch rows
  int **assmatrix = (int**)malloc((numStartPatches+1)*sizeof(int*));
  for (i = 0; i <= numStartPatches; i++) {
    assmatrix[i] = (int*)malloc(numStartPatches*sizeof(int));
    //initialize it with zeros
    for (j = 0; j < numStartPatches; j++ ) {
      assmatrix[i][j] = 0; 
    }
  }  
  //not the matrix of who arrived where (how many individuals from patch x
  //arrive in patch y....

  int p;
  int ***totArrivalmatrix = (int***)malloc((numStartPatches + 1) * sizeof(int**));
  for (p = 0; p <= numStartPatches; p++) {
    totArrivalmatrix[p] = (int**)malloc((numStartPatches) * sizeof(int*));
    for (i = 0; i < numStartPatches; i++) {
      totArrivalmatrix[p][i] = (int*)malloc(numStartPatches * sizeof(int));
      //initialize it with zeros
      for (j = 0; j < numStartPatches; j++ ) { 
        totArrivalmatrix[p][i][j] = 0; 
      }
    }  
  }

  //get space for the resulting map
  resmap = (int**)malloc(nrows * sizeof(int*));
  for (i = 0; i < nrows; i++) {
    resmap[i] = (int*)malloc(ncols * sizeof(int)); //there's a good chap
    //now initialize it
    for (j = 0; j < ncols; j++) {
      resmap[i][j] = 0;
    }
  }
  finalresult = (int**)malloc(nrows * sizeof(int*));
  for (i = 0; i < nrows; i++) {
    finalresult[i] = (int*)malloc(ncols * sizeof(int)); //there's a good chap
    //now initialize it
    for (j = 0; j < ncols; j++) {
      finalresult[i][j] = 0;
    }
  }

  int **totalEmigrants = (int**)malloc((numStartPatches + 1) * sizeof(int*));
  for (i = 0; i <= numStartPatches; i++) {
    totalEmigrants[i] = (int*)malloc(numStartPatches * sizeof(int));
  }

  for (missingPatchIndex = 0; 
      missingPatchIndex <= numStartPatches; 
      missingPatchIndex++) {
    for (i = 0; i < numStartPatches; i++) { 
      totalEmigrants[missingPatchIndex][i] = 0;
    }

    /* for each row */
    for (row = 1; row < nrows - 1; row++) {
      CELL c;
      Rast_get_row(in_startpatches_fd, startpatches_cell_buf, row, startpatches_data_type);
      /* process the data */
      col = 0;
      while (col < ncols-1) {
        col++;
        /* use different function for each data type */
        c = ((CELL *) startpatches_cell_buf)[col];
        if ((c != NONHABITAT) && (Rast_is_c_null_value(&c) != 1)) {
          //search ahead to check how far the patch row reaches to start search
          //path for several cells at once
          int tocol = col + 1;
          CELL cc;
          while ( (tocol < ncols - 1) && 
              (c == (cc = ((CELL *) startpatches_cell_buf)[tocol]))) {
            tocol++; 
          }
          tocol--; //it was one too far

          if (missingPatchIndex > 0) {
            if (c != patchhash[missingPatchIndex - 1]) {
              searchPath(in_landscape_fd, in_replacement_fd, col, tocol, row, 
                  nrows, ncols, landscape_data_type, replacement_data_type, 1, 
                  numberOfIndividualsReleasedPerCell,
                  transitionMatrix, parmatrix.rows, steplengthvec,
                  habitattypesvec, energycostsvec, initial_energy, 
                  perceptual_range, patchhash, numStartPatches, 
                  missingPatchIndex-1, starthabitatvec[missingPatchIndex-1], 
                  resmap);
            }
          } else {  //the default case without leaving anything out of 
                    //the picture
            searchPath(in_landscape_fd, in_replacement_fd, col, tocol, row, 
                nrows, ncols, landscape_data_type, replacement_data_type, 1, 
                numberOfIndividualsReleasedPerCell,
                transitionMatrix, parmatrix.rows, steplengthvec,
                habitattypesvec, energycostsvec, initial_energy, 
                perceptual_range, patchhash, numStartPatches, 
                missingPatchIndex-1, starthabitatvec[missingPatchIndex], 
                resmap);
          }
          col = tocol;  
          for (i = 0; i < nrows; i++) {
            CELL cr;
            Rast_get_row(in_startpatches_fd, startpatches_res_cell_buf, i, startpatches_data_type);
            for (j = 0; j < ncols; j++) {
              if (resmap[i][j] > 0) {
                finalresult[i][j] += resmap[i][j]; 
                cr = ((CELL *) startpatches_res_cell_buf)[j];
                source = findInIntVec(c, patchhash, numStartPatches);
                destination = NONVALID;
                if ((cr != NONHABITAT) && (Rast_is_c_null_value(&cr) != 1)) {
                  destination = findInIntVec(cr, patchhash, numStartPatches);
                  if ((destination != NONVALID) && 
                      (source != NONVALID) && 
                      (source != destination)) {  //don't count an individual 
                    totArrivalmatrix[missingPatchIndex][source][destination] += 
                      resmap[i][j];               //that stayed in it's own patch
                    assmatrix[missingPatchIndex][destination] += resmap[i][j];
                  }
                }

                if ((source != NONVALID) && (source != destination)) {
                  totalEmigrants[missingPatchIndex][source]+=resmap[i][j];
                }
                resmap[i][j] = 0;
              }//if (resmap[i][j] > 0)
            }//for (j = 0; j < ncols; j++)
          }//for (i = 0; i < nrows; i++)
        }//if ((c != NONHABITAT) && ( G_is_c_null_value(&c) != 1))
      }//for (col = 1; col < ncols-1; col++)
    } // for (row = 1; row < nrows-1; row++)
    if (finalresult != NULL) {
      //leaves space for 10^9 patches
      char *bf = (char*)malloc((strlen(RESULTMAPFILENAME) + 10) * 
        sizeof(char)); 
      if (missingPatchIndex == 0) { 
        sprintf(bf, "%s%d", RESULTMAPFILENAME, 0);
      } else { 
        sprintf(
          bf, 
          "%s%d", 
          RESULTMAPFILENAME, 
          patchhash[missingPatchIndex - 1]
        );
      }

      fprintf(stderr, "start writing dispersalmap <%s>\n",bf);
      out_result_fd = Rast_open_new(bf, landscape_data_type);
      //out_result_fd = G_open_raster_new(bf,landscape_data_type);
      //output_cell_buf = G_allocate_raster_buf(landscape_data_type);
      output_cell_buf = Rast_allocate_buf(landscape_data_type);

      for (i = 0; i < nrows; i++) {
        output_cell_buf = (CELL*)finalresult[i];
        //G_put_map_row(out_result_fd, output_cell_buf);
        Rast_put_row(out_result_fd, output_cell_buf, landscape_data_type);
        for (j = 0; j < ncols; j++) {
          finalresult[i][j] = 0;
        }
      }
      //Rast_close(out_result_fd);
      Rast_close(out_result_fd);
      free(bf);
    }
  }
  
  FILE *fp_out;
  char outname[BUFSIZE];
  int **totalImmigrantsMissing = 
    (int**)malloc((numStartPatches + 1) * sizeof(int*));
  for (i = 0; i <= numStartPatches; i++) {
    totalImmigrantsMissing[i] = (int*)malloc(numStartPatches*sizeof(int));
  }
  //write the rest of the results:
  for (p = 0; p <= numStartPatches; p++) {
    double *percentSuccessfulEmigrants = 
      (double*)malloc(numStartPatches * sizeof(double));
    int **arrivalmatrix = totArrivalmatrix[p];
    //write the total number of immigrants into the patches
    if (p == 0) {
      sprintf( 
        outname, 
        "%s/%s-total.csv", 
        outputpath,
        RESULT_IMMIGRANTS_PER_EMIGRANTS_FILENAME
      );
    } else { 
      sprintf(
        outname, 
        "%s/%s-without-%d.csv",
        outputpath, 
        RESULT_IMMIGRANTS_PER_EMIGRANTS_FILENAME, 
        patchhash[p - 1]
      );
    }
    fprintf(
      stderr, 
      "Writing total ratio of immigrants to emigrants to\n %s\n", 
      outname
    );
    if ((fp_out = fopen(outname, "w")) == NULL) {
      G_fatal_error(_("Cant write to file <%s>"), outname);
    } else {
      //write the header:
      fprintf(fp_out, "fromhabitat ");
      for (i = 0; i < numStartPatches; i++) {
        fprintf(fp_out, "tohabitat-%d ", patchhash[i]);
      }
      fprintf(fp_out, "totalEmigrants\n");

      double tval = 0;
      for (i = 0; i < numStartPatches; i++) {
        percentSuccessfulEmigrants[i] = 0;
        fprintf(fp_out, "%d ", patchhash[i]);
        for (j = 0; j < numStartPatches; j++) {
          if (totalEmigrants[p][i] > 0) {
            tval = (1.0 * arrivalmatrix[i][j]) / totalEmigrants[p][i];
            percentSuccessfulEmigrants[i] += tval; 
          }
          else { 
            tval = 0;
          }
          fprintf(fp_out, "%7.5lf ", tval);
        }
        fprintf(fp_out, " %d\n", totalEmigrants[p][i]);  
      }
      fclose(fp_out);
    }

    //write the maps for totalEmigrants:
    for (row = 0; row < nrows; row++) {
      CELL c;
      /*if (G_get_c_raster_row(
          in_startpatches_fd, 
          startpatches_cell_buf, 
          row) < 0) {*/
      Rast_get_row(in_startpatches_fd, startpatches_cell_buf, row, startpatches_data_type);
      /* process the data */
      for (col = 0; col < ncols; col++) {
        finalresult[row][col] = 0;
        /* use different function for each data type */
        c = ((CELL *) startpatches_cell_buf)[col];
        if ((c != NONHABITAT) && (Rast_is_c_null_value(&c) != 1)) {
          i = findInIntVec(c, patchhash, numStartPatches);
          if (i >= 0) {
            finalresult[row][col] = totalEmigrants[p][i];
          }
        }
      } 
    } 
    //write the emigrant map:
    if (p == 0) { 
      sprintf(outname, "%s0", TOTALEMIGRANTMAP);
    }
    else {
      sprintf(outname, "%s%d", TOTALEMIGRANTMAP, patchhash[p - 1]);
    }

    fprintf(stderr, "start writing totalEmigrantsmap <%s>\n", outname);

    out_result_fd = Rast_open_new(outname, landscape_data_type);
    output_cell_buf = Rast_allocate_buf(landscape_data_type);
    //out_result_fd = G_open_raster_new(outname,landscape_data_type);
    //output_cell_buf = G_allocate_raster_buf(landscape_data_type);
    for (i = 0; i < nrows; i++) {
      output_cell_buf = (CELL*)finalresult[i];
      Rast_put_row(out_result_fd, output_cell_buf, landscape_data_type);
      //G_put_map_row(out_result_fd, output_cell_buf);
      for (j = 0; j < ncols; j++) {
        finalresult[i][j] = 0;
      }
    }
    Rast_close(out_result_fd);

    //prepare map with the successfull emigrants
    for (row = 0; row < nrows; row++) {
      CELL c;
      Rast_get_row(
        in_startpatches_fd, 
        startpatches_cell_buf, 
        row,
        startpatches_data_type
      );
      /* process the data */
      for (col = 0; col < ncols; col++) {
        finalresult[row][col] = 0;
        /* use different function for each data type */
        c = ((CELL *) startpatches_cell_buf)[col];
        if ((c != NONHABITAT) && (Rast_is_c_null_value(&c) != 1)) {
          i = findInIntVec(c, patchhash, numStartPatches);
          if (i >= 0) {
            finalresult[row][col] = (int)floor(
              percentSuccessfulEmigrants[i] * 100
            );
          }
        }
      }
    } 

    //write map with the successfull emigrants
    if (p == 0) { 
      sprintf(outname, "%s0", PERCENTSUCCESSFULEMIGRANTSMAP );
    }
    else {
      sprintf(outname, "%s%d",PERCENTSUCCESSFULEMIGRANTSMAP, patchhash[p - 1]);
    }

    fprintf( 
      stderr, 
      "start writing percent_successful_emigrantsmap <%s>\n",
      outname
    );

    out_result_fd = Rast_open_new(outname, landscape_data_type);
    //out_result_fd = G_open_raster_new(outname,landscape_data_type);
    output_cell_buf = Rast_allocate_buf(landscape_data_type);
    //output_cell_buf = G_allocate_raster_buf(landscape_data_type);
    for (i = 0; i < nrows; i++) {
      output_cell_buf = (CELL*)finalresult[i];
      //G_put_map_row(out_result_fd, output_cell_buf);
      Rast_put_row(out_result_fd, output_cell_buf, landscape_data_type);
      for (j = 0; j < ncols; j++) { 
        finalresult[i][j] = 0;
      }
    }
    Rast_close(out_result_fd);

    //assess the relative number of immigrants
    if (p == 0) {
      sprintf(
        outname, 
        "%s/%s-total.csv", 
        outputpath, 
        RESULTRELATIVENUMIMMIGRANTS
      );
    } else { 
      sprintf(
        outname, 
        "%s/%s-without-%d.csv",
        outputpath, 
        RESULTRELATIVENUMIMMIGRANTS, 
        patchhash[p - 1]
      ); 
    }
    fprintf(
      stderr, 
      "Writing relative number of immigrants to\n %s\n", 
      outname
    );
    if ((fp_out = fopen(outname, "w")) == NULL) {
      G_fatal_error(_("Cant write to file <%s>"), outname);
    } else {
      int totalImmigrants;
      fprintf(fp_out, "tohabitat ");
      for (i = 0; i < numStartPatches; i++) {
        fprintf(fp_out, "fromhabitat-%d ", patchhash[i]);
      }
      fprintf(fp_out, "totalImmigrants\n");
      for (i = 0; i < numStartPatches; i++) {
        totalImmigrants = 0;
        //count all immigrants that died in patch i
        for (j = 0; j < numStartPatches; j++) {
          totalImmigrants += arrivalmatrix[j][i]; //sum over column (from j to i)
        }
        totalImmigrantsMissing[p][i] = totalImmigrants; //all immigrants that died in patch i
      }

      for (i = 0; i < numStartPatches; i++) {
        //start writing relative number of immigrants  
        fprintf(fp_out, "%d ", patchhash[i]);
        for (j = 0; j < numStartPatches; j++) {
          if (totalImmigrantsMissing[p][i] > 0) {
            fprintf(
              fp_out, 
              "%7.5lf ", 
              ((1.0 * arrivalmatrix[j][i]) / totalImmigrantsMissing[p][i]));
          } else {
            fprintf(fp_out, "%7.5lf ", 0.0);
          }
        }
        fprintf(fp_out, " %d\n", totalImmigrantsMissing[p][i]);  
      }//end for (i = 0; i < numStartPatches; i++)
      fclose(fp_out);
    }//end if ((fp_out = fopen(outname, "w")) == NULL)

    //write the maps for totalImmigrants:
    for (row = 0; row < nrows; row++) {
      CELL c;
      Rast_get_row(
        in_startpatches_fd, 
        startpatches_cell_buf, 
        row,
        startpatches_data_type
      );
      /* process the data */
      for (col = 0; col < ncols; col++) {
        finalresult[row][col] = 0;
        /* use different function for each data type */
        c = ((CELL *) startpatches_cell_buf)[col];
        if ((c != NONHABITAT) && (Rast_is_c_null_value(&c) != 1)) {
          i = findInIntVec(c, patchhash, numStartPatches);
          if (i >= 0) {
            finalresult[row][col] = totalImmigrantsMissing[p][i];
          }
        }
      }//end for (col = 0; col < ncols; col++)
    }

    //write the immigrant map:
    if (p == 0) {
      sprintf(outname, "%s0", TOTALIMMIGRANTMAP); //makeing up a name
    }  
    else {
      sprintf(outname, "%s%d", TOTALIMMIGRANTMAP, patchhash[p - 1]);  //add the missing patch number
    }

    fprintf(stderr, "start writing totalimmigrantsmap <%s>\n",outname);

    out_result_fd = Rast_open_new(outname, landscape_data_type);
    //out_result_fd = G_open_raster_new(outname,landscape_data_type);
    //output_cell_buf = G_allocate_raster_buf(landscape_data_type);
    output_cell_buf = Rast_allocate_buf(landscape_data_type);
    for (i = 0; i < nrows; i++) {
      output_cell_buf = (CELL*)finalresult[i];
      //G_put_map_row(out_result_fd, output_cell_buf);
      Rast_put_row(out_result_fd, output_cell_buf, landscape_data_type);
      for (j = 0; j < ncols; j++) { 
        finalresult[i][j] = 0;
      }
    }
    Rast_close(out_result_fd);
  }

  int **totalArrivalMatrix =  totArrivalmatrix[0];
  //now compare that with the total result:
  //and summarize the difference matrix
  int flagsum = 0;
  flagsum += (flag_assessmentmodePAMI == FLAGSWITCH_ON) ? PAMI : 0;
  flagsum += (flag_assessmentmodePAMM == FLAGSWITCH_ON) ? PAMM : 0;

  for (p = 1; p <= numStartPatches; p++) {
    double patchValue = 0;
    if (flag_assessmentmodePAMI == FLAGSWITCH_ON)  {
      for (i = 0; i < numStartPatches; i++) {
        if (totalImmigrantsMissing[0][i]>0) {                            
          patchValue += (1.0*totalArrivalMatrix[p-1][i])/totalImmigrantsMissing[0][i];                          
        }
      }
    }  
    if (flag_assessmentmodePAMM == FLAGSWITCH_ON) {
      int **arrivalmatrix = totArrivalmatrix[p]; //now check what happens if p is missing
      for (i = 0; i < numStartPatches; i++) {
        for (j = 0; j < numStartPatches; j++) {
          if ((i != (p-1)) && (j != (p-1)) && (totalEmigrants[0][i]>0)) { 
            patchValue -= ((1.0*arrivalmatrix[i][j] - totalArrivalMatrix[i][j])/totalEmigrants[0][i]); 
          }
        }
      }
    }
    starthabitatvec[p - 1].value = patchValue;
    fprintf(
      stderr, 
      "patchvalue of patch %d is %lf \n", 
      patchhash[p - 1],  
      patchValue
    );
  }

  //and write it into a mapfile, so that each patch gets a numerical value
  for (row = 0; row < nrows; row++) {
    CELL c;
    Rast_get_row(
      in_startpatches_fd, 
      startpatches_cell_buf, 
      row,
      startpatches_data_type
    );
    /* process the data */
    for (col = 0; col < ncols; col++) {
      finalresult[row][col] = 0;
      /* use different function for each data type */
      c = ((CELL *) startpatches_cell_buf)[col];
      if ((c != NONHABITAT) && (Rast_is_c_null_value(&c) != 1)) {
        i = findInIntVec(c, patchhash, numStartPatches);
        if (i >= 0) {
          if (starthabitatvec[i].value > 0) {
            finalresult[row][col] = 
              (int)floor(starthabitatvec[i].value*PATCHASSESSMENTMULTIPLIER);
          } else { 
            finalresult[row][col] = NULLREPLACEMENT; 
          }
          //fprintf(stderr, "row %d, col %d, value %6.5lf\n", row, col, starthabitatvec[i].value);
        }
      }
    }
  }//end for (row = 0; row < nrows; row++) 

  //the actual writing starts
  const char *pname;
  switch(flagsum) {
    case PAMM:
      pname = PATCHASSESSMENTMAP_PAMM;
      break;
    case PAMI:
      pname = PATCHASSESSMENTMAP_PAMI;
      break;
    case PAMMI:
      pname = PATCHASSESSMENTMAP_PAMMI;
  }

  fprintf(stderr, "start writing patchassessmentmap <%s>\n",pname);
  out_result_fd = Rast_open_new(pname, landscape_data_type);
  output_cell_buf = Rast_allocate_buf(landscape_data_type);
  //out_result_fd = G_open_raster_new(pname,landscape_data_type);
  //output_cell_buf = G_allocate_raster_buf(landscape_data_type);
  for (i = 0; i < nrows; i++) {
    output_cell_buf = (CELL*)finalresult[i];
    Rast_put_row(out_result_fd, output_cell_buf, landscape_data_type);
    //G_put_map_row(out_result_fd, output_cell_buf);
    for (j = 0; j < ncols; j++) {
      finalresult[i][j] = 0;
    }
  }
  Rast_close(out_result_fd);

  for (i = 0; i <= numStartPatches; i++) { 
    free(totalImmigrantsMissing[i]);
  }
  free(totalImmigrantsMissing); //release memory of the immigrant count

  //now write the assmatrix into a file
  sprintf(outname, "%s/%s", outputpath, RESULTDISPERSALFILENAME);

  if ((fp_out = fopen(outname, "w")) == NULL) {
    G_fatal_error(_("Cant write to file <%s>"), outname);
  } else {
    //write the header:
    fprintf(fp_out, "withouthabitat ");
    for (i = 0; i < numStartPatches; i++) { 
      fprintf(fp_out, "ended_in_habitat-%d ", patchhash[i]);
    }
    fprintf(fp_out, "\n");

    for (i = 0; i <= numStartPatches; i++) {
      if (i == 0) { 
        fprintf(fp_out, "%d ", 0);
      } else { 
        fprintf(fp_out, "%d ", patchhash[i - 1]);  
      }
      for (j = 0; j < numStartPatches; j++) {
        fprintf(fp_out, "%d ", assmatrix[i][j] );
      }
      fprintf(fp_out, "\n");  
    }
    fclose(fp_out);
  }//end else if ((fp_out = fopen(outname, "w")) == NULL)

  if (finalresult != NULL) {
    fprintf(stderr, "releasing memory\n");
    for (i = 0; i < nrows; i++) { 
      free(finalresult[i]); 
    }
    free(finalresult);
  }

  for (i=0; i<=numStartPatches; i++) {
    for (j=0; j<numStartPatches; j++) { 
      free(totArrivalmatrix[i][j]);
    }
    free(totArrivalmatrix[i]);
  }
  free(totArrivalmatrix); 

  //free memory
  for (i = 0; i < nrows; i++) {
    free(resmap[i]);
  }
  free(resmap);
  for (i=0; i<numStartPatches; i++) {
    free(assmatrix[i]);
  }
  free(assmatrix);
  for (i=0; i<numStartPatches; i++) {
    free(totalEmigrants[i]);
  }
  free(totalEmigrants);
  free(numReleasedPerPatch);
  free(patchhash);
  G_free(startpatches_cell_buf);
  G_free(startpatches_res_cell_buf);
  G_free(replacement_cell_buf);
  fprintf(stderr, "done!\n\n");  
}

int main(int argc, char *argv[]) {

  double initial_energy_level;

  int number_of_individual_released_by_cell;

  int i,
      j;

  char *landscape_mapset, 
       *startpatches_mapset,
       *replacement_mapset,
       *landscape_in_file, 
       *startpatches_in_file,  
       *replacement_in_file,  
       *parameters_in_file, 
       *out_file;

  char startpatches_layer[64], 
       replacement_layer[64], 
       landscape_layer[64], 
       parameters_file[64], 
       outputpath[OUTPUTPATHLENGTH];

  int landscape_fd, 
      startpatches_fd, 
      replacement_fd, 
      parameters_fd;

  double NS_fac, 
         EW_fac, 
         DIAG_fac, 
         H_DIAG_fac, 
         V_DIAG_fac;

  int col, 
      row, 
      nrows, 
      ncols, 
      nr_of_habitat_types; 

  double *habitattypesvec, 
         *energycostsvec, 
         *frictionvec, 
         *steplengthvec;       

  int *ncols_parameterfile, 
      *nrows_parameterfile;

  struct Option *opt1, 
                *opt2, 
                *opt3, 
                *opt4, 
                *opt5, 
                *opt6;

  struct Flag *flag1, 
              *flag2;

  struct GModule *module; 

  int number_of_rows,  
      number_of_cols;  /*number of cols and rows in the landscape*/

  FILE *param_filepointer;     /*file pointer to the parameter file*/

  double **tmp_parameter_file_entry_matrix; /* a matrix to be filled with all entries of the
                                             * parameter file */
  double **transition_matrix;  /*the transition matrix for the habitat types (how)
                                *likely is it to go from one type to another type 
                                *matrix properties: square (n=m), t[i,i]=1*/    

  double *habitat_types_vec;  /* list of habitat types present in the landscape
                               * which have to correspond to the number of rows
                               * in the transition matrix*/

  double *friction_habitat_types_vec; /* vector of the amount of friction the habitat
                                       * tpyes pose to the individual passing through 
                                       * a cell of that type: 
                                       * friction[1] = friction in habitat type 1 ...*/

  double *steplength_habitat_types_vec; /* vector of the distance an individual covers 
                                         * by doing one step when passing through habitat
                                         * of the corresponding type:
                                         * steplength[2] = steplength in habitat type 2 ... */

  double energy_level_individuals;      /* the level of energy that an individual has at
                                         * its disposal. it's being reduced by friction of
                                         * the type of habitat it is currently in */

  double perceptual_range = 1.0;  
  parameters parmatrix;

  int flag_assessmentmodePAMM, flag_assessmentmodePAMI;   
  /* flag whether or not to include relative 
   * streams of migrants into the calculation
   * of the patchassessmentvalue*/

  G_gisinit(argv[0]);

  module = G_define_module();
  module->description =
    _("Outputs a raster map layer showing the "
        "the number of immigrants of a given landscape cell. ");

  opt1 = G_define_option();
  opt1->key = "landscape_types";
  opt1->type = TYPE_STRING;
  opt1->required = YES;
  opt1->gisprompt = "old,cell,raster";
  opt1->description = _("Name of landscape type input raster map");

  opt2 = G_define_option();
  opt2->key = "start_patches";
  opt2->type = TYPE_STRING;
  opt2->required = YES;
  opt2->gisprompt = "old,cell,raster";
  opt2->description = _("Name of input raster map of habitats");

  opt6 = G_define_option();
  opt6->key = "replacement_types";
  opt6->type = TYPE_STRING;
  opt6->required = YES;
  opt6->gisprompt = "old,cell,raster";
  opt6->description = _("Name of input raster map of replacement habitats");

  opt3 = G_define_option();
  opt3->key = "parameter_file";
  opt3->type = TYPE_STRING;
  opt3->required = YES;
  opt3->gisprompt = "old_file,file,input";
  opt3->description = _("Name of parameter file");

  opt4 = G_define_option();
  opt4->key = "energylevel";
  opt4->type = TYPE_DOUBLE;
  opt4->required = YES;
  opt4->description = _("Initial energy level of individuals");

  opt5 = G_define_option();
  opt5->key = "nrofindividualsreleasedpercell";
  opt5->type = TYPE_INTEGER;
  opt5->required = YES;
  opt5->description = _("Number of individuals to be released per start cell");

  flag1 = G_define_flag();
  flag1->key = 'm';
  flag1->description = _("PAMM: Map includes changes in migrant streams of other patches");

  flag2 = G_define_flag();
  flag2->key = 'i';
  flag2->description = _("PAMI: Map includes changes in nr. of immigrants into other patches");

  if (G_parser(argc, argv))
    exit(EXIT_FAILURE);

  G_get_window(&window);

  /*  Find north-south, east_west and diagonal factors */
  EW_fac = window.ew_res; /*East-West; Must be the physical distance */
  NS_fac = window.ns_res; /*North-South*/
  /*diagonal factor*/
  DIAG_fac    = (double)sqrt((double)(NS_fac * NS_fac + EW_fac * EW_fac)); 
  V_DIAG_fac  = (double)sqrt((double)(4 * NS_fac * NS_fac + EW_fac * EW_fac)); 
  H_DIAG_fac  = (double)sqrt((double)(NS_fac * NS_fac + 4 * EW_fac * EW_fac));

  strcpy(landscape_layer, opt1->answer);
  strcpy(startpatches_layer, opt2->answer);
  strcpy(replacement_layer, opt6->answer);
  strcpy(parameters_file, opt3->answer);

  getcwd(outputpath, OUTPUTPATHLENGTH );
  fprintf(stderr, "outputpath: %s\n", outputpath);

  initial_energy_level = atof(opt4->answer);
  number_of_individual_released_by_cell = atoi(opt5->answer);

  flag_assessmentmodePAMM = flag1->answer ? FLAGSWITCH_ON  : FLAGSWITCH_OFF; 
  flag_assessmentmodePAMI = flag2->answer ? FLAGSWITCH_ON  : FLAGSWITCH_OFF; 

  if ((flag_assessmentmodePAMM == FLAGSWITCH_OFF) && 
      (flag_assessmentmodePAMI == FLAGSWITCH_OFF)) {
    G_fatal_error(_("At least one Flag (PAMI and/or PAMM) must be set!")); 
  }

  fprintf(
    stderr, 
    "flags: %d %d \n", 
    flag_assessmentmodePAMI, 
    flag_assessmentmodePAMM 
  );

  param_filepointer = fopen(parameters_file, "r");

  if ( param_filepointer) {
    parmatrix = read_parameter_file_contents(param_filepointer);
    fprintf(stderr, "back\n"); 
    for (i = 0; i < parmatrix.rows; i++) {
      for (j = 0; j < parmatrix.cols; j++) {
        fprintf(stderr, "%lf ",parmatrix.entrymatrix[i][j]);  
      }
      fprintf(stderr, "::\n");  
    }
  }
  else {
    G_fatal_error(_("%s - cannot be openend!"), parameters_file); 
  }

  landscape_in_file = G_tempfile();
  startpatches_in_file = G_tempfile();
  parameters_in_file = G_tempfile();
  out_file = G_tempfile();

  doTheModeling( landscape_layer, replacement_layer, startpatches_layer, parmatrix, steplengthvec, frictionvec, energycostsvec, habitattypesvec,number_of_individual_released_by_cell, initial_energy_level, perceptual_range, outputpath, flag_assessmentmodePAMM, flag_assessmentmodePAMI);

  exit (EXIT_SUCCESS);
}
