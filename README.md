# SomSeg
 Image analysis codes used for segmentation of somite boundaries

This repository contains image analysis codes for segmenting somite and notochord boundaries in images of zebrafish embryo development.
The reference preprint for this work can be found here: https://doi.org/10.1101/2020.08.14.251645
Map-projected data set for this work can be found here: http://doi.org/10.5281/zenodo.4146919

The repository contains the following image analysis codes:

1. Get_embryo_surface_SN.m
This code detects nuclei in fused images obtained from a multiview microscope and performs a sphere fit of the point cloud

2. Map_projection_SN.m
This code peforms map projection (equidistant cylindrical or Mercator) of fused spherical data

3. Map_projection_pixel_map_SN.m
This code determines the pixel size in the projected map

4. Segment_somites_notochord_SN.m
This code segments notochord and somite boundaries in map projected images.
This code should perform equally well in non-map-projected images as well. You will have to comment/uncomment a few lines as mentioned within this code

5. Determine_median_position_boundary_SN.m
This code determines median positions of segmented somite boundaries and notochord


