/*
 * This is an example module call:
 * r.dispersal -m -i landscape_types=matrix start_patches=start_patches replacement_types=repl_patches parameter_file=/home/andi/dev/grass7_trunk/raster/r.dispersal/param.in energylevel=4 nrofindividualsreleasedpercell=1000 
*/
/* the following searchmask is being applied:
 * S S S
 * S X S
 * 0 0 0
 */      
#define DEFINITIONS 1
#define XPOS 0
#define YPOS 1
#define INITIAL_NUMBER_OF_STARTPATCHES 200
#define PATCHASSESSMENTMULTIPLIER 1000
#define SEARCHMASK_WIDTH 5
#define SEARCHMASK_HEIGHT 2
#define SURROUNDINGMASK_WIDTH 8
#define SURROUNDINGMASK_HEIGHT 2        
#define OUTPUTPATHLENGTH 256
#define NONVALID -1
#define ROW_ALLOC_INCREMENT 5
#define ROW_ALLOC_MAX 200
#define MAX_MEM_SIZE_USED_BY_BUFFER 1024*1024*80/sizeof(CELL)
#define TRUE 1
#define FALSE 0
#define FLAGSWITCH_ON 17
#define FLAGSWITCH_OFF 42

#define PAMM 1
#define PAMI 3
#define PAMMI 4

#define PRINT 0
#define DEBUG if(PRINT==1) 

#define COUNT 0
#define COUNTMIGRANTS if(COUNT == 1)

#define NULLREPLACEMENT 0.00001

const   int surroundingmask[SURROUNDINGMASK_WIDTH][SURROUNDINGMASK_HEIGHT] = 
  {{-1, 0},
   {-1, 1},
   { 0, 1}, 
   { 1, 1}, 
   { 1, 0}, 
   { 1,-1}, 
   { 0,-1}, 
   {-1,-1}};

const int LEFTRANGE = -2;
const int RIGHTRANGE = 2;

const char *RESULTMAPFILENAME = "dispersalmap_missing_";
const char *PATCHASSESSMENTMAP_PAMI  = "PAMI";
const char *PATCHASSESSMENTMAP_PAMM  = "PAMM";
const char *PATCHASSESSMENTMAP_PAMMI = "PAMMI";

const char *TOTALEMIGRANTMAP = "total_emigrants_missing_";
const char *TOTALIMMIGRANTMAP = "total_immigrants_missing_";
const char *RESULT_IMMIGRANTS_PER_EMIGRANTS_FILENAME = "immigrants_per_emigrants";
const char *RESULTRELATIVENUMIMMIGRANTS = "immigrants_per_totalimmigrants";
const char *RESULTDISPERSALFILENAME = "dispersalnumbers.csv";
const char *PERCENTSUCCESSFULEMIGRANTSMAP = "percent_successful_emigrants_missing_";

const int HABITATINDEX    = 0;
const int STEPLENGTHINDEX = 1;
const int ENERGYINDEX     = 2; 

const int NONHABITAT     = 0;
const int MAXNUMSTEPS    = 500;

struct point {
  int x, y;
};


struct individual {
  double energylevel;
  int sourcehabitat;
  double x, y;
  double dx, dy;
};

struct parameters {
  int rows, cols;  
  double **entrymatrix;
};


struct habitatextension {
  int habitattype;
  int xmin, xmax, ymin, ymax;
};

typedef struct habitatextension habitatextension;

struct habitat {
  habitatextension extension;
  int numcells;
  double value;
  int **entrymatrix;
};

typedef struct parameters parameters;
typedef struct habitat habitat;
typedef struct individual individual;
typedef struct point point;

#define frand() ((double) rand() / (RAND_MAX+1.0))
