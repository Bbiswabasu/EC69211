
# Digital Image And Video Processing Lab
## Experiment 2

## By
- _Bbiswabasu Roy , 19EC39007_
- _Jothi Prakash , 19EC39023_

# Problem Statement 

Write C++ modular functions to read, perform operations, and write BMP image files. All functions must support 24-bit RGB and 8-bit grayscale image formats.

# Tech

- C++ 11.0 or above

# Header Files
- iostream
- fstream
- cmath
- cstring

# How to run

- Open terminal and run `g++ Experiment2.cpp - o Experiment2 && ./Experiment2`
- Choose one of the options (1 or 2) based on the function to be tested
- Enter the name of the image file as <filename>.bmp (file should exist in same directory) 

# Generated Input and Output files
- The files taken for input are placed int input/ directory
- The files obtained after applying various operations are stored in output/ directory

# Functions

# Struct BMPdata
- Stores all the required header information and the pixel array of the bmp file
- Includes all the necessary functions to read, write and display the data

# Struct BMPdata::readBMP(FilePointer fptr)
- Takes as input the file pointer of the opened BMP Image File
- Reads the necessary header information and the pixel array of the Image

# Struct BMPdata::writeBMP(FilePointer fptr)
- Writes all the header and pixel array information to the output BMP file

# Struct BMPdata::printHeader()
- Prints all the header information in a structured manner for easy reading

# scaleImage(BMPdata &bmp_data, double scaleX, double scaleY)
- Takes input the bmp_data and the x and y scaling factor
- Performs the scaling operation and interpolation using Nearest Neighbor Interpolation
- Updates the required changes in the bmp_data header and pixel array

# rotateImage(BMPdata &bmp_data,double angle,int type)
- Takes input the bmp_data and the angle of rotation as well as the Interpolation type
- Type 1 -> Bilinear , Type 2 -> Nearest Neighbour
- Performs the requires rotation and the interpolation based on the parameters
- Updates the required changes in the bmp_data header and pixel array





