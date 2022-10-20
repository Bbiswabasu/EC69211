
# Digital Image And Video Processing Lab
## Experiment 6

## By
- _Bbiswabasu Roy , 19EC39007_
- _Jothi Prakash , 19EC39023_

# Problem Statement 

Write Python functions to perform the following operations on the given test image,
ricegrains_mono.bmp. All functions must support binary images.

Q1. Make separate functions for erosion, dilation, opening, and closing of binary images

# Tech

- Python 2+
- Open CV 

# Imports
- math
- numpy
- cv2

# How to run

- Open ipynb file and run all cells
- Enter the choice of kernel to apply the operations
- Enter the choise of Morphological operation you want to perform on the image

# Functions

# ErodeBinary(x,b)
- Takes as input the image - x, and the kernel - b.
- Performs erosion on the image using the defined algorithm.
- Returns the modified image.

# DilateBinary(x,b)
- Takes as input the image - x, and the kernel -b.
- Performs dilation on the image using the defined algorithm.
- Returns the modified image.

# kernel_gen(dim)
- Takes as input the dimension of the kernel.
- Generates a square kernel with the given dimensions.
- Returns the generated kernel.
