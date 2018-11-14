

-- Overview --------

The program extracts a collection of image statistics called minPS from an image. minPS is a simplified version of statistical parameters used in Portilla & Simoncelli algorithm (PS statistics). 

-- Usage -----------

You should install and add path to texture synthesis codes provided by Portilla & Simoncelli.

im2minPS.m - convert images into minPS statistics.
minPS2V4cell.m - convert minPS statistics into hypothetical V4 neurons' activities. Note that activities are normalized and do not correspond to the actual firing rates.

show_ac_pattern.m - show principal components of 'Position' statistics. minPS includes coefficients for those principal components (see minPS_list.pdf).




G Okazawa (2014)





