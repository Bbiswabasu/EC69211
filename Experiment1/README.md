
# Digital Image And Video Processing Lab
## Experiment 1

## By
- _Bbiswabasu Roy , 19EC39007_
- _Jothi Prakash , 19EC39023_

# Problem Statement 

Write python modular functions to perform the following task. All functions must support 24-bit RGB, 8-bit grayscale and 8- bit color indexed image formats.

# Tech

- Python 3
- Matplotlib Library (to display the images)

# How to run

- Open terminal and run `python3 experiment.py`
- Choose one of the options (1,2 or 3) based on the function to be tested
- Enter the input file path. All input images are stored under `input/` directory, so enter file path as `input/<filename>.bmp`
- If you are asked to enter path of output file, enter `output/<filename>.bmp`

# Functions

## read_bmp(filename)
- Takes the filename with .bmp extension and reads the entire file using file pointers.
- Stores all the header files related data in the bmp_header variable
- Stores the entire image in a 3D list in RGB format in the pixel_array variable
- return bmp_header, pixel_array

## write_bmp(filename,bmp_header,pixel_array)
- Creates a new .bmp extension file with the name filename
- Stores the bmp_header and the pixel_array in binary format to the created file

## display_bmp(bmp_header,pixel_array)
- Displays the stored bmp_header dictionary
- Displays the image formed by pixel_array using matplotlib

## manipulate_channel(filename,remove_color)
- Takes the the name of the .bmp file to be read with the filename parameter
- Takes the remove_color paramter of the color that will be masked in the RGB channel
- remove_color - 'r' = Red , 'g' = Green , 'b' = Blue

## bytes_to_int(b)
- Converts the given bytearray to integer

## int_to_bytes
- Converts the given integer to bytearray




