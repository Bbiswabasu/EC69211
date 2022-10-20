
# Digital Image And Video Processing Lab
## Experiment 4

## By
- _Bbiswabasu Roy , 19EC39007_
- _Jothi Prakash , 19EC39023_

# Problem Statement 

Q1. Write python modular functions/subroutines to design spatial filters - Mean, Median, Prewitt, Laplacian, Sobel kernels (horizontal, vertical, diagonal), Gaussian Blur, Laplacian of Gaussian on a stack of grayscale images (say, 15 images per stack).

Q2. Create a filter called Gaussian_Unblur to undo the effects of blurring.

# Tech

- Python 2+
- Open CV 

# Imports
- math
- numpy
- cv2
- os

# How to run

- Open terminal and run `python3 experiment4.py`
- Enter 1 if you want to apply filter and 2 if you want to unblur image
- Enter the path on which to operate
- Then choose appropriate options to execute

# Functions

# mean_filter_kernel(n)
- Generate the mean filter kernel of dimension nxn
- Returns the generated kernel

# median_filter(img,filter_size)
- Takes as input the img and the dimension of the filter
- Applies the required filtering and returns the new_img

# prewitt_filter_kernel(orientation)
- returns the prewitt filter for the desired orientation

# laplacian_filter_kernel()
- returns the laplacian filter kernel

# sobel_filter_kernel(orientation)
- returns the sobel filter for the desired orientation

# gaussian_filter_kernel(n,sigma)
- returns the gaussian kernel with the given n and sigma

# convolve(img,kernel)
- Takes the given image and type of kernel as input
- Flips the kernel horizontally and vertically and performs the convolution
- Returns the convoled resultant image 

# gaussian_unblur(img)
- Takes the img on which gaussian filter has been apllied
- Performs the unblur operation to sharpen the image 
- Returns the gaussian unblurred image



