g.region -p
d.mon x0
d.rast glc00
d.what.rast
r.mapcalc "glc00.1 = if(glc00==0,null(),glc00)"
r.colors map=glc00.1 rast=glc00
d.rast glc00.1
d.zoom
exit
r.mapcalc "glc00.2 = if(glc00.1==26,null(),glc00.1)"
r.colors map=glc00.2 raster=glc00.1
r.colors map=glc00.2 rast=glc00.1
r.info glc00.2
g.region -p
r.le.setup
r.le.pixel --h
r.le.pixel map=glc00.2 sam=m att=b1,b2,b3,b4 div:d1,d2,d3,d4 te1=m7 te2=t1,t2,t3,t4,t5 edg=e1
d.rast africa_mask
g.copy rast=africa_mask,MASK
r.le.pixel map=glc00.2 sam=m att=b1,b2,b3,b4 div=d1,d2,d3,d4 te1=m7 te2=t1,t2,t3,t4,t5 edg=e1
d.rast b1
d.rast d1
d.rast d2
d.rast d3
d.rast t5
d.what.rast
d.rast e1
d.m.
r.in.gdal
exit
r.support --h
r.support map=b1
r.support map=b1
r.info b1
r.info b2
r.info t5
ecit
exit
r.info glc00.2
r.support glc00.2
r.info glc00.2
r.info b1
r.out.gdal
r.proj
exit
r.info b1
g.region -p
r.in.gdal
v.proj
v.proj input=eilidat_61_90.1 location=srtm_gtopo30 mapset=topo_analysis.50m_step output=eilidat_resample_del
r.proj input=eilidat_61_90.1 location=srtm_gtopo30 mapset=topo_analysis.50m_step output=eilidat_resample_del
r.proj input=eilidat_61_90.1 location=srtm_gtopo30 mapset=topo_analysis.50m_steps output=eilidat_resample_del
g.region align=eilidat_61_90.1
g.region align=eilidat_resample_del
g.region -p
g.region align=eilidat_resample_del
g.copy rast=MASK,x_MASK
g.remove rast=MASK
g.region align=eilidat_resample_del
g.region -p
r.proj
g.remove rast=eilidat_resample_del
g.region -p
g.region nsres=0:02:29.99999 ewres=0:02:29.99999
g.region -p
r.resample input=b1 output=glc_mean10_resample
g.copy rast=xMASK,MASK
g.copy rast=x_MASK,MASK
r.resample input=b1 output=glc_mean10_resample
g.remove rast=glc_mean10_resample
r.resample input=b1 output=glc_mean10_resample
r.resample input=d1 output=glc_div10_resample
r.resample input=e1 output=glc_edge10_resample
r.resample input=t5 output=glc_contrast10_resample
r.out.tiff --h
r.out.tiff inpu=glc_mean10_resample output=glc_mean10_resample 
pwd
cd jakob/10km_radius/
r.out.tiff inpu=glc_mean10_resample output=glc_mean10_resample 
r.out.tiff inpu=glc_mean10_resample output=glc_mean10_resample -t
r.out.tiff inpu=glc_div10_resample output=glc_div10_resample -t
r.out.tiff inpu=glc_contrast10_resample output=glc_contrast10_resample -t
r.out.tiff inpu=glc_edge10_resample output=glc_edge10_resample -t
exit
g.region -p
g.region align=glc00
g.region -p
r.out.tiff input=b1 output=glc_div10 -t
r.out.tiff input=b1 output=glc_mean10 -t
r.out.tiff input=d1 output=glc_div10 -t
r.out.tiff input=e1 output=glc_edge10 -t
r.out.tiff input=e1 output=glc_edge10 -t
exit
exit
exit
g.remove st
exit
r.out.gdal2
r.out.
exit
g.remvoe
exit
g.remove
g.mremove rast=b*
g.mremove rast=d*
g.list rast
g.mremove rast=t*
g.mremove rast=e*
g.list rast
g.list vect
g.remove rast=glc00.1,glc_div10_resample,x_MASK
g.remove rast=africa_mask,glc_edge10_resample,glc_contrast10_resample,glc_mean10_resample
g.list vect
g.list rast
g.remove rast=glc00
g.rename rast=glc00.2,GLC2000
r.support GLC2000
r.info GLC2000
g.remove
exit
g.remove
exit
exit
