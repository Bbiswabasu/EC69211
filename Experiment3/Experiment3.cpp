#include <iostream>
#include <cmath>
#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
using namespace std;
using namespace cv;

#define PI (4*atan(1))
struct my_complex
{
    double real,imag;
    my_complex()
    {
        real = imag = 0;
    }
    my_complex(double t1,double t2)
    {
        real=t1;
        imag=t2;
    }

    my_complex operator+(my_complex b)
    {
        return my_complex(real+b.real,imag+b.imag);
    }

    my_complex operator-(my_complex b)
    {
        return my_complex(real-b.real,imag-b.imag);
    }

    my_complex operator*(my_complex b)
    {
        return my_complex(real*b.real-imag*b.imag,real*b.imag+imag*b.real);
    }

    my_complex operator/(double t)
    {
        return my_complex(real/t,imag/t);
    }

    double getMag()
    {
        double t=real*real+imag*imag;
        t=sqrt(t);
        return t;
    }

    double getPhase()
    {
        return atan(imag/real);
    }

};

ostream &operator<<(ostream &os, my_complex t)
{
    os<<t.real << " + j" << t.imag;
    return os;
}

void fft_1d(my_complex *a, int n, bool inverse)
{
    if (n == 1)
        return;

    my_complex a_even[n/2], a_odd[n/2];
    for (int i = 0; 2*i < n; i++)
    {
        a_even[i] = a[2*i];
        a_odd[i] = a[2*i+1];
    }

    // recursively compute FFT for each half
    fft_1d(a_even, n/2, inverse);
    fft_1d(a_odd, n/2, inverse);

    // merge the FFT of two halves as per the equations
    double ang = 2*PI/n * (inverse ? -1 : 1);
    my_complex w(1,0), wn(cos(ang), sin(ang));
    for (int i = 0; 2*i < n; i++)
    {
        a[i] = a_even[i] + w*a_odd[i];
        a[i + n/2] = a_even[i] - w*a_odd[i];
        if (inverse)
        {
            // scale down by N for inverse FFT
            // since this operation is applied at all stages, it will be scaled down by N at the end
            a[i] = a[i]/2;
            a[i+n/2] = a[i+n/2]/2;
        }
        w = w*wn;
    }
}

void fftshift_image(my_complex **pixel_array, int height, int width)
{
    for(int col=0;col<width/2;col++)
    {
        for(int row=0;row<height;row++)
        {
            swap(pixel_array[row][col],pixel_array[row][col+width/2]);
        }
    }

    for(int row=0;row<height/2;row++)
    {
        for(int col=0;col<width;col++)
        {
            swap(pixel_array[row][col],pixel_array[row+height/2][col]);
        }
    }
}

void fft_image(my_complex **pixel_array, int height, int width, bool inverse)
{
    if(inverse)
       fftshift_image(pixel_array,height,width);

    // apply FFT along each column
    for(int col=0;col<width;col++)
    {
        my_complex current_col[height];
        for(int row=0;row<height;row++)
        {
            current_col[row]=pixel_array[row][col];
        }
        fft_1d(current_col, height, inverse);
        for(int row=0;row<height;row++)
        {
            pixel_array[row][col]=current_col[row];
        }
    }

    // apply FFT along each row
    for(int row=0;row<height;row++)
    {
        my_complex current_row[width];
        for(int col=0;col<width;col++)
        {
            current_row[col]=pixel_array[row][col];
        }
        fft_1d(current_row, height, inverse);
        for(int col=0;col<width;col++)
        {
            pixel_array[row][col]=current_row[col];
        }
    }
    if(!inverse)
        fftshift_image(pixel_array,height,width);
}

void q2(my_complex **pixel_array,int height,int width)
{
    // apply (-1)^(x+y) operation on pixel_array
    for(int row=0;row<height;row++)
    {
        for(int col=0;col<width;col++)
        {
            if((row+col)%2)
                pixel_array[row][col]=pixel_array[row][col]*my_complex(-1,0);
        }
    }
    
    // take FFT and then complex conjugate
    fft_image(pixel_array,height,width,0);
    for(int row=0;row<height;row++)
    {
        for(int col=0;col<width;col++)
        {
            pixel_array[row][col].imag*=-1;
        }
    }
    
    // take inverse FFT
    fft_image(pixel_array,height,width,1);
    
    // apply (-1)^(x+y) on real part of pixel_array elements
    for(int row=0;row<height;row++)
    {
        for(int col=0;col<width;col++)
        {
            if((row+col)%2)
                pixel_array[row][col].real*=-1;
        }
    }
}

int main()
{
    string image_path;
    cout<<"Enter image path: ";
    cin>>image_path;
    Mat img = imread(image_path, IMREAD_GRAYSCALE);
    if(img.empty())
    {
        std::cout << "Could not read the image: " << image_path << std::endl;
        return 1;
    }
    
    // Read the image pixel data in an array of type unsigned char
    unsigned char **array = new unsigned char*[img.rows];
    for (int i=0; i<img.rows; ++i)
        array[i] = new unsigned char[img.cols*img.channels()];

    for (int i=0; i<img.rows; ++i)
        array[i] = img.ptr<unsigned char>(i);
    
    // Load the same pixel data in an array which stores complex numbers as its elements
    my_complex **pixel_array=new my_complex*[img.rows];
    for(int i=0;i<img.rows;i++)
    {
        pixel_array[i]=new my_complex[img.cols];
        for(int j=0;j<img.cols;j++)
        {
            pixel_array[i][j]=my_complex(array[i][j],0.0);
        }
    }

    cout<<"Which part to execute (1 or 2): ";
    int x;
    cin>>x;
    if(x==1)
    {
    	// Apply fft and result gets stored in pixel_array
        fft_image(pixel_array,img.rows,img.cols,0);

	// Load the magnitude and phase of coefficients in Mat object so that it can be written to new jpg file
        Mat mag_plot(img.rows,img.cols,CV_16U);
        Mat phase_plot(img.rows,img.cols,CV_16U);
        for(int i=0;i<img.rows;i++)
        {
            for(int j=0;j<img.cols;j++)
            {
                mag_plot.at<short>(i,j)=25*log10(pixel_array[i][j].getMag());
                phase_plot.at<short>(i,j)=255/PI*pixel_array[i][j].getPhase();
            }
        }
        image_path.pop_back();
        image_path.pop_back();
        image_path.pop_back();
        image_path.pop_back();
        string opFile=image_path+"_mag_out.jpg";
        imwrite(opFile, mag_plot);
        opFile=image_path+"_phase_out.jpg";
        imwrite(opFile,phase_plot);
    }
    else
    {
    	// Call q2() function to apply the required sequence of operations and then store new pixel_array in Mat object
        Mat dip_op(img.rows,img.cols,CV_16U);
        q2(pixel_array,img.rows,img.cols);
        for(int i=0;i<img.rows;i++)
        {
            for(int j=0;j<img.cols;j++)
            {
                dip_op.at<short>(i,j)=pixel_array[i][j].real;
            }
        }
        image_path.pop_back();
        image_path.pop_back();
        image_path.pop_back();
        image_path.pop_back();
        string opFile=image_path+"_out.jpg";
        imwrite(opFile, dip_op);
    }  
}
