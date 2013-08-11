g.move
g.copy
r.support habitat
g.remove rast=habitat,habitat_cat
g.copy
g.remove rast=habitat_base
r.mapcalc "habitat_base = habitat"
r.mapcalc "habitat_start_cat = habitat_cat"
r.mapcalc "habitat_start = habitat"
g.remove rast=habitat_base
r.mapcalc "habitat_base = habitat"
r.mapcalc "habitat_base_cat = habitat_cat"
g.remove
r.mapcalc "habitat_start_1 = if(habitat_base_cat == 579
r.extract
r.mapcalc "habitat_start_1 = if(habitat_base_cat == 579,602,1,null())"
r.mapcalc "habitat_start_1 = if(habitat_base_cat == 579;602,1,null())"
r.to.vect rast=habitat_base_cat vect=habitat_base_cat --h
r.to.vect output=habitat_base_cat input=habitat_base_cat feature=area
v.extract
v.to.rast
r.to.vect
v.to.rast
r.mapcalc "habitat_start_1_post_mortem_A = if(habitat_start_1==6,2,3)"
r.mapcalc "habitat_start_1_post_mortem_B = if(habitat_start_1==6,1,3)"
r.mapcalc "habitat_start_cat_1 = habitat_start_1"
r.mapcalc "habitat_start_1 = if(habitat_start_1>0,3,nul())"
r.mapcalc "habitat_start_1 = if(habitat_start_1>0,3,null())"
r.mapcalc "habitat_start_1_post_mortem_C = if(habitat_start_1==1,1,3)"
r.mapcalc "habitat_start_1_post_mortem_C = if(habitat_start_cat_1==1,1,3)"
r.mapcalc "habitat_start_1_post_mortem_D = if(habitat_start_cat_1==3,1,3)"
r.mapcalc "habitat_start_1_post_mortem_D = if(habitat_start_cat_1==3&&habitat_start_cat_1==4&&habitat_start_cat_1==7,1,3)"
r.mapcalc "habitat_start_1_post_mortem_D = if(habitat_start_cat_1==3||habitat_start_cat_1==4||habitat_start_cat_1==7,1,3)"
r.mapcalc "habitat_start_1_post_mortem_E = if(habitat_start_cat_1==5||habitat_start_cat_1==6,1,3)"
r.mapcalc "habitat_start_1_post_mortem_F = if(habitat_start_cat_1==8,1,3)"
g.list rast
r.mapcalc "habitat_start_1_post_mortem_A-F = if(habitat_start_cat_1>0,1,3)"
r.mapcalc "habitat_start_1_post_mortem_AF = if(habitat_start_cat_1>0,1,3)"
r.mapcalc "habitat_start_1_post_mortem_AFb = if(habitat_start_1_post_mortem_AF>2||habitat_start_cat_1==1,1,3)"
r.mapcalc "habitat_start_1_post_mortem_AF = if(habitat_start_cat_1==5,2,habitat_start_1_post_mortem_AF)"
r.mapcalc "habitat_start_1_post_mortem_AFb = if(habitat_start_cat_1==5,2,habitat_start_1_post_mortem_AF)"
r.mapcalc "habitat_post_mortem = habitat_start_1_post_mortem_AFb"
r.mapcalc "habitat_start_1_post_mortem_AFc = if(habitat_start_cat_1==5&&habitat_start_cat_1==6&&habitat_start_cat_1==8,2,habitat_start_1_post_mortem_AF)"
r.mapcalc "habitat_start_1_post_mortem_AFc = if(habitat_start_cat_1==5||habitat_start_cat_1==6||habitat_start_cat_1==8,2,habitat_start_1_post_mortem_AF)"
r.mapcalc "habitat_post_mortem_2 = habitat_start_1_post_mortem_AFc"
g.remove
g.remove
r.mapcalc "habitat_start = if(habitat_start_1>0,1,null())"
g.remove
r.mapcalc "habitat_start_cat = habitat_start_cat_1"
g.remove
g.remove
exit
exit
exit
r.dispersal landscape_types=habitat_base@march07 start_patches=habitat_start@march07 replacement_types=habitat_post_mortem@march07 parameter_file=/Users/andi/grass_metapop2/parameterfile.in energylevel=1 nrofindividualsreleasedpercell=10
r.dispersal -m -i landscape_types=habitat_base@march07 start_patches=habitat_start@march07 replacement_types=habitat_post_mortem@march07 parameter_file=/Users/andi/grass_metapop2/parameterfile.in energylevel=1 nrofindividualsreleasedpercell=10
r.dispersal -m -i --verbose landscape_types=habitat_base@march07 start_patches=habitat_start@march07 replacement_types=habitat_post_mortem@march07 parameter_file=/Users/andi/grass_metapop2/parameterfile.in energylevel=2 nrofindividualsreleasedpercell=100
