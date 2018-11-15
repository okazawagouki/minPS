# Code to generate minPS stats in Okazawa et al. (2015) PNAS

---

## Overview

The program extracts a collection of image statistics called minPS from an image. minPS is a simplified version of statistical parameters used in Portilla & Simoncelli (2000) used to fit texture selectivity of V4 neurons in Okazawa et al. (2015).

## Requirement
The code is tested under Matlab 2018a.It requires [matlab pyramid tool](http://www.cns.nyu.edu/~eero/STEERPYR/) and [texture synthesis tool](http://www.cns.nyu.edu/~lcv/texture/) developped in Simoncelli lab. Those codes are already included in this package.
## Package

[im2minPS.m](./im2minPS.m) - main code to convert images into minPS statistics.
[tutorial_im2minPS.m](./tutorial_im2minPS.m) - tutorial code of im2minPS.m including some details necessary to correctly run the program.


## Distribution of this code

Because this package contains codes developped in Simoncelli lab, distribution of this package is prohibited.


## Reference

Okazawa G, Tajima S, Komatsu H (2015) Image statistics underlying natural texture selectivity of neurons in macaque V4. *PNAS* 112(4):E351-60

Portilla J, Simoncelli EP (2000) A Parametric Texture Model based on Joint Statistics of Complex Wavelet Coefficients. *International Journal of Computer Vision* 40(1):49-71


---
Gouki Okazawa (2018/11/15)





