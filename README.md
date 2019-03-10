# Fractal Image Compression Using Upper Bound on Scaling Parameter

This repository contains an implementation of the paper "Fractal image compression using upper bound on scaling parameter" - ([link to the ScienceDirect page](https://www.sciencedirect.com/science/article/pii/S0960077917304691))
```
Roy, S.K., Kumar, S., Chanda, B., Chaudhuri, B.B. and Banerjee, S., 2018. Fractal image compression 
using upper bound on scaling parameter. Chaos, Solitons & Fractals, 106, pp.16-22.
```            
            
## Description
This paper presents a novel approach to calculate the affine parameters of fractal encoding, in order to reduce its computational complexity. A simple but efficient approximation of the scaling parameter is derived which satisfies all properties necessary to achieve convergence. It allows us to substitute to the costly process of matrix multiplication with a simple division of two numbers. We have also proposed a modified horizontal-vertical (HV) block partitioning scheme, and some new ways to improve the en- coding time and decoded quality, over their conventional counterparts. Experiments on standard images show that our approach yields performance similar to the state-of-the-art fractal based image compres- sion methods, in much less time


## Experimental Data
 
 The effectiveness of the proposed method has been tested on bechmark images from [USC-SIPI Image Database](http://sipi.usc.edu/database/).

## Results

Comparison of proposed scheme with other conventional fractal coders, w.r.t. PSNR, BPP, and Time (in sec), for four standard images.


## Citation

If you use this code or a derivative thereof in your research, we would appreciate a citation to the original paper:

```
@article{roy2018fractal,
        title={Fractal image compression using upper bound on scaling parameter},
        author={Roy, Swalpa Kumar and Kumar, Siddharth and Chanda, Bhabatosh and Chaudhuri, Bidyut B and Banerjee, Soumitro},
        journal={Chaos, Solitons \& Fractals},
        volume={106},
        pages={16--22},
        year={2018},
        publisher={Elsevier}}
```

