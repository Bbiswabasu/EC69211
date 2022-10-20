
# Digital Image And Video Processing Lab
## Experiment 3

## By
- _Bbiswabasu Roy , 19EC39007_
- _Jothi Prakash , 19EC39023_

# Problem Statement 

Q1. Write your own function for computing Fast Fourier Transform (FFT) and inverse FFT of an Image with any random spatial dimension

Q2. Read the image dip.tiff, and perform the required operation on the image.

# Tech

- C++ 11.0 or above
- Open CV 

# Header Files
- iostream
- cmath
- opencv2/core.hpp
- opencv2/imgcodecs.hpp
- opencv2/highgui.hpp

# How to run

- Open terminal and run `g++ Experiment3.cpp -o  Experiment && ./Experiment`
- Enter the file path on which to operate
- Choose one of the options (1 or 2) based on the function to be tested

# Functions

# struct my_complex
- Contains the custom defined complex data type to handle the computations
- real and imag variables are used to denote the real and imaginary part of complex number
- Major functions such as add, sub, multiply and div are overloaded for complex operations

# fft_1D(a,n,inverse)
- Takes input the 1D array as a and the length of it as n
- Performs the Fast Fourier Transform algorithm to compute coefficients in O(NlogN)
- inverse -> 0 denotes fft and inverse -> 1 denotes inverse fft
- Performs inplace computation in the a array itself

# fft_image(pixel_array,row,col,inverse)
- Takes the pixel_array and the row and column dimensions in row and col variables
- inverse -> 0 denotes fft and inverse -> 1 denotes inverse fft
- Performs inplace calculation and updates the pixel_array

# fft_shift_image(pixel_array,int height,int width)
- Takes the pixel_array and the row and column as heigth and width respectively
- Shifts the fft result to make the result centered around the DC component
- Preforms inplace calculation and updates the pixel_array

# q2(pixel_array,height,width)
- Takes as input the pixel_array and the heigth and width of the array
- Performs the operations mentioned in Question 2 of assignment
- Performs inplace computation and updates the result in pixel_array





