# grass7_r_dispersal

## DESCRIPTION

<code>r.dispersal</code> is a module for Grass7 to estimate 
the importance of habitat patches as stepping stones for animal dispersal. 
For a series of habitat patches, each is taken in turn as a start patch for 
a mechanistic animal movement/survival model and each other patch is 
substitued in turn by a replacement patch. The number of
arriving animals is recorded for every patch.


The program will be run non-interactively if the user
specifies program arguments (see OPTIONS) on the command
line.  

## Command line options
<code>r.dispersal -m -i --verbose landscape_types=habitat_base@march07 start_patches=habitat_start@march07 replacement_types=habitat_post_mortem@march07 parameter_file=~/grassdata/parameterfile.in energylevel=2 nrofindividualsreleasedpercell=100</code>

## OPTIONS
The user must specify the names of the raster map layers to
be used for *base habitat types*, *start_patches* and *replacement patches*, 
the *method* used (i.e., Proximity Index), the *min*, the *max* and the *keyval* 
of the class of interest of the input raster map.


### Patchassessment Options

the *-m* and *-i* flags define that changes in number of migrants and changes in number of immigrants of the focus patch are calculated separetely. If you select both flags than it will add up these two changes.
Generally these values are multiplied by 1000 for display reasons.

### Parameter file:
The *parameter file* holds all information concerning extinction probability, edge-traversability-probability, matrix values. 

If you have four landscape cell types, your parameterfile looks like this:

<pre>
4 7
0.01 0.01 0.4 0.48 0 1 0.001
0.1 0.01 0.1 0.79 1 1 0.01
0.01 0.01 0.38 0.6 2 1 0.001
0.1 0.01 0.7 0.19 3 1 0.0001
</pre>

The IDs of the four cell types are in then in the fifth column and the probabilities to cross from one cell type to another are in the first four columns. 
The transition probabilities (dispersal costs) are calculated from horizontal category to vertical category, e.g.:

<pre>
0.01 0.01 0.4 0.48 0 1 0.001
0.1 0.01 0.1 0.79 1 1 0.01
0.01 0.01 0.38 0.6 2 1 0.001
0.1 0.01 0.7 0.19 3 1 0.0001
</pre>
means, that the probability to cross from category 0 to 3 is 48%, and the probability to cross from 3 to 2 is 70%.

The 5th row corresponds to the matrix categories. The 6th row it the step length (not yet used) and the 7th row holds information concerning enery costs per matrix type and step (0 = reduction of energy level per step is zero, 1= reduction of energy level by 1 in one step. The higher the parameter setting "energylevel", the more steps individual animals can take.

## NOTES

If you have problems with memory allocation you have to change 
the following line in the file definitions.h:

<code>#define MAX_MEM_SIZE_USED_BY_BUFFER 1024*1024*20/sizeof(CELL)</code>

the value 20 (approx. MB) can be changed e.g. 80, 100 up to the amount of your memory.

### Authors

Dr. Andreas Gros, Dr. Martin Wegmann

Last changed: 20130804
