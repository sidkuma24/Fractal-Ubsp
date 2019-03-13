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
 
 
 | Image   | Quality | Fisher | Polvere | Saupe | Hurtgen | Saupe-MC | Proposed-Quadtree | Propsoed-HV |
|---------|---------|--------|---------|-------|---------|----------|-------------------|-------------|
| Baboon  | PSNR    | 24.28  | 24.43   | 25.54 | 24.41   | 24.45    | 25.21             | 25.80       |
|         | BPP     | 0.60   | 0.61    | 0.60  | 0.60    | 0.60     | 0.60              | 0.60        |
|         | Time    | 1.60   | 2.80    | 25.40 | 2.8     | 3.37     | 1.57              | 8.90        |
| Boat    | PSNR    | 33.95  | 34.82   | 35.61 | 34.62   | 34.63    | 34.21             | 35.83       |
|         | BPP     | 0.60   | 0.60    | 0.60  | 0.60    | 0.61     | 0.60              | 0.60        |
|         | Time    | 2.40   | 4.70    | 12.40 | 4.07    | 4.60     | 2.19              | 11.20       |
| Lenna   | PSNR    | 35.58  | 36.10   | 36.20 | 35.90   | 36.90    | 36.12             | 37.67       |
|         | BPP     | 0.60   | 0.60    | 0.60  | 0.59    | 0.60     | 0.60              | 0.60        |
|         | Time    | 2.20   | 3.60    | 13.20 | 3.50    | 4.30     | 1.90              | 12.19       |
| Man     | PSNR    | 30.19  | 30.50   | 31.80 | 30.45   | 30.44    | 31.25             | 32.60       |
|         | BPP     | 0.60   | 0.60    | 0.60  | 0.61    | 0.60     | 0.60              | 0.71        |
|         | Time    | 1.8    | 2.9     | 15.4  | 3.30    | 3.90     | 1.15              | 11.60       |
| Average | PSNR    | 31     | 31.46   | 31.25 | 31.35   | 31.35    | 31.69             | 37.72       |
|         | BPP     | 0.60   | 0.60    | 0.60  | 0.60    | 0.60     | 0.60              | 0.60        |
|         | Time    | 2.00   | 3.50    | 14.10 | 3.41    | 4.04     | 1.70              | 10.97       |

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

