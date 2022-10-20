
# Digital Image And Video Processing Lab
## Experiment 5

## By
- _Bbiswabasu Roy , 19EC39007_
- _Jothi Prakash , 19EC39023_

# Problem Statement 

Q1. Perform the following frequency domain filtering by writing your own Python functions (LPF= Low Pass Filter, HPF = High Pass Filter).
a. Ideal LPF, Ideal HPF
b. Gaussian LPF, Gaussian HPF
c. Butterworth LPF, Butterworth HPF

Q2. Generate a hybrid image using "einstein.png" and "marilyn.png".

Q3. Remove the noise from the input image and store the denoised image.

# Tech

- Python 2+
- Open CV 

# Imports
- math
- numpy
- cv2
- matplotlib

# How to run

- Open ipynb file and run all cells
- Enter 1 if you want to apply Ideal, Gaussian or Butterworth Filter
- Enter 2 if you want to generate the hybrid image
- Enter 3 if you want to remove noise from the noised image

# Functions

# q1()
- Take as input the file path of the image, filter name, cut off type and cut off frequency.
- Generates the output after applying selected filter on the given image

# q2()
- Passes "einstein.png" through a high pass gaussian filter
- Passes "marlyne.png" through a low pass gaussian filter
- Combines the result to generate the hybrid image

# q3()
- Take as input the file path the image
- Generates the de noised image from the noised image by using the median filter

# ideal_filter(f0,type,rows,cols)
- f0 - Cut Off Frequency, type - High / Low
- rows, cols are the dimensions
- Generates the ideal filter in frequency domain

# gaussian_filter(f0,type,rows,cols)
- f0 - Cut Off Frequency, type - High / Low
- rows, cols are the dimensions
- Generates the gaussian filter in frequency domain

# butterworth_filter(f0,type,rows,cols,order)
- f0 - Cut Off Frequency, type - High / Low, order - order of filter
- rows, cols are the dimensions
- Generates the butterworth filter in frequency domain

# fft_image(img,inverse)
- Performs the Fast Fourier Transform operation on the given image img
- inverse - 0 is FFT and inverse - 1 is inverse FFT

# fftshift_image(img) 
- Performs fft shift on the image by centering the low frequency components

# pad_image(img)
- pads extra rows and columns to img to make the rows and cols a power of 2
