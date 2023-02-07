https://github.com/MichalMarekHoppe/Spatial-analysis---calculating-delta-between-cellular-phenotypes.git
author: Michal Marek Hoppe (mmlhoppe@gmail.com)
source: https://www.medrxiv.org/content/10.1101/2020.10.20.20216101

Spatial analysis - calculating delta % between cellular phenotypes among k nearest neighbours

MYC-BCL2-BCL6_spatial_sample.txt -  sample data from three high-power microscopic images of DLBCL cases. input data consists of
a) id - sample id
b) label - cellular pheotype. every cell has to belong to a group. In this work, CD20+ cells have been phenotyped into subpopulations depending on their MYC, BCL2, BCL6 oncogene expression status. CD20 negatiove cells constistute all non-B cells in the analyzed image.
d) x and y - spatial coordiantes for each cell.

Number of k neighbours need to can be defined at the beginning of the process.

Summary results sheet provides results for interactions among all phenotypes for each sample. 