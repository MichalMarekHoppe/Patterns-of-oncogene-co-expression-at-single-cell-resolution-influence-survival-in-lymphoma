Spatial Analysis of DLBCL Tissue Microarrays: Pair Correlation Function

Gayatri Kumar
06/02/2023


# Contact: gkumar@mdanderson.org or gayatri@iisc.ac.in

# Scripts used  for the analysis of the two datasets.

# 1. MaskWithCellBoundaries_.ipynb - Python notebook overlays the X/Y coordinates of the cells on the  RGB images with segmentation masks.

# 2. The two Geojson annotations are exported from QuPath for - 
#       - RGB+ segmentation  mask 
#       - RGB+ segmentation  mask with overlaid points

# 3. Read_geojson_pointpattern.Rmd - 
#    - Point patterns with convex windows are generated for each image
#    - The geojson annotations replace the convex window for each image.

# 4. PairCorrelationFunction_Clustering.Rmd - 
#    Calculates clustering using the pair correlation function for a range of radius 
#    values

Sample input files included, corresponding to two patients, 03_08 and 10_02:
 
-.jpg files of 'RGB+ segmenation mask' images 
-.tiff files of 'RGB+ segmenation mask with overlaid points'images
-.csv files with the X/Y coordinates and subpopulation labels for each cell per image
