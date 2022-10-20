#include <iostream>
#include <fstream>
#include<cmath>
#include <cstring>

using namespace std;

#define PI (4*atan(1))

struct BMPdata
{
    unsigned char type[3];
    unsigned int size_bytes;
    unsigned int reserved_1;
    unsigned int reserved_2;
    unsigned int offset;
    unsigned int dib_header_size;
    unsigned int width_pixels;
    unsigned int height_pixels;
    unsigned int colour_planes;
    unsigned int bits_per_pixel;
    unsigned int compression_method;
    unsigned int raw_image_size_bytes;
    unsigned int horizontal_resolution;
    unsigned int vertical_resolution;
    unsigned int num_colors;
    unsigned int important_colors;
    unsigned int *color_table;
    unsigned int ***pixel_array;

    unsigned int bytes_to_int(unsigned char *b,int len)
    {
        unsigned int n=0,shift = 0;
        for(int i=0;i<len;i++)
        {
            unsigned int t = b[i];
            n = n | (t<<shift);
            shift += 8;
        }
        return n;
    }

    void readBMP(FILE* fptr)
    {
        unsigned char buf[4];
        
        fread(type, sizeof(unsigned char), 2, fptr);
        type[2]='\0';
        char bm[3]="BM";
        if(strcmp((char*)type,"BM")!=0)
        {
            cout<<"Not BMP file\n";
            return;
        }
        
        fread(buf, sizeof(unsigned char), 4, fptr);
        size_bytes = bytes_to_int(buf, 4);
        
        fread(buf, sizeof(unsigned char), 2, fptr);
        reserved_1 = bytes_to_int(buf, 2);
        
        fread(buf, sizeof(unsigned char), 2, fptr);
        reserved_2 = bytes_to_int(buf, 2);
        
        fread(buf, sizeof(unsigned char), 4, fptr);
        offset = bytes_to_int(buf, 4);
        
        fread(buf, sizeof(unsigned char), 4, fptr);
        dib_header_size = bytes_to_int(buf, 4);
        
        fread(buf, sizeof(unsigned char), 4, fptr);
        width_pixels = bytes_to_int(buf, 4);
        
        fread(buf, sizeof(unsigned char), 4, fptr);
        height_pixels = bytes_to_int(buf, 4);
        
        fread(buf, sizeof(unsigned char), 2, fptr);
        colour_planes = bytes_to_int(buf, 2);
        
        fread(buf, sizeof(unsigned char), 2, fptr);
        bits_per_pixel = bytes_to_int(buf, 2);
        
        fread(buf, sizeof(unsigned char), 4, fptr);
        compression_method = bytes_to_int(buf, 4);
        
        fread(buf, sizeof(unsigned char), 4, fptr);
        raw_image_size_bytes = bytes_to_int(buf, 4);
        
        fread(buf, sizeof(unsigned char), 4, fptr);
        horizontal_resolution = bytes_to_int(buf, 4);
        
        fread(buf, sizeof(unsigned char), 4, fptr);
        vertical_resolution = bytes_to_int(buf, 4);
        
        fread(buf, sizeof(unsigned char), 4, fptr);
        num_colors = bytes_to_int(buf, 4);
        
        fread(buf, sizeof(unsigned char), 4, fptr);
        important_colors = bytes_to_int(buf, 4);
        
        color_table = new unsigned int[num_colors];
        for(int i=0;i<num_colors;i++)
        {
            fread(buf, sizeof(unsigned char), 4, fptr);
            color_table[i] = bytes_to_int(buf, 4);
        }

        pixel_array = new unsigned int**[height_pixels];
        for(int i=0;i<height_pixels;i++)
        {
            pixel_array[height_pixels-1-i] = new unsigned int*[width_pixels];
            for(int j=0;j<width_pixels;j++)
            {
                pixel_array[height_pixels-1-i][j] = new unsigned int[3];
                fread(buf, sizeof(unsigned char), bits_per_pixel/8, fptr);
                unsigned int current_color = bytes_to_int(buf, bits_per_pixel/8);
                for(int k=0;k<3;k++)
                {
                    pixel_array[height_pixels-1-i][j][2-k] = current_color&(0xff);
                    current_color>>=8;
                }
            }
        }
    }

    void int_to_bytes(unsigned int n, unsigned char *b, int len)
    {
        for(int i=0;i<len;i++)
        {
            b[i] = n&0x000000ff;
            n>>=8;
        }
    }
    void writeBMP(FILE* fptr)
    {
        unsigned char buf[4];
        fwrite(type, sizeof(unsigned char), 2, fptr);
        
        int_to_bytes(size_bytes, buf, 4);
        fwrite(buf, sizeof(unsigned char), 4, fptr);

        int_to_bytes(reserved_1, buf, 2);
        fwrite(buf, sizeof(unsigned char), 2, fptr);
        
        int_to_bytes(reserved_2, buf, 2);
        fwrite(buf, sizeof(unsigned char), 2, fptr);
        
        int_to_bytes(offset, buf, 4);
        fwrite(buf, sizeof(unsigned char), 4, fptr);
        
        int_to_bytes(dib_header_size, buf, 4);
        fwrite(buf, sizeof(unsigned char), 4, fptr);
        
        int_to_bytes(width_pixels, buf, 4);
        fwrite(buf, sizeof(unsigned char), 4, fptr);
        
        int_to_bytes(height_pixels, buf, 4);
        fwrite(buf, sizeof(unsigned char), 4, fptr);
        
        int_to_bytes(colour_planes, buf, 2);
        fwrite(buf, sizeof(unsigned char), 2, fptr);
        
        int_to_bytes(bits_per_pixel, buf, 2);
        fwrite(buf, sizeof(unsigned char), 2, fptr);
        
        int_to_bytes(compression_method, buf, 4);
        fwrite(buf, sizeof(unsigned char), 4, fptr);
        
        int_to_bytes(raw_image_size_bytes, buf, 4);
        fwrite(buf, sizeof(unsigned char), 4, fptr);
        
        int_to_bytes(horizontal_resolution, buf, 4);
        fwrite(buf, sizeof(unsigned char), 4, fptr);
        
        int_to_bytes(vertical_resolution, buf, 4);
        fwrite(buf, sizeof(unsigned char), 4, fptr);
        
        int_to_bytes(num_colors, buf, 4);
        fwrite(buf, sizeof(unsigned char), 4, fptr);
        
        int_to_bytes(important_colors, buf, 4);
        fwrite(buf, sizeof(unsigned char), 4, fptr);
        
        for(int i=0;i<num_colors;i++)
        {
            int_to_bytes(color_table[i], buf, 4);
            fwrite(buf, sizeof(unsigned char), 4, fptr);
        }

        for(int i=0;i<height_pixels;i++)
        {
            for(int j=0;j<width_pixels;j++)
            {
                unsigned int current_color = pixel_array[height_pixels-1-i][j][2] | (pixel_array[height_pixels-1-i][j][1]<<8) | (pixel_array[height_pixels-1-i][j][0]<<16);
                int_to_bytes(current_color, buf, 3);
                for(int k=0;k<bits_per_pixel/8;k++)
                {
                    fwrite(&(buf[k]), sizeof(unsigned char), 1, fptr);
                }
            }
        }
    }

    void printHeader()
    {
        cout<<"type: "<<type[0]<<type[1]<<"\n"
            <<"size_bytes: "<<size_bytes<<"\n"
            <<"reserved_1: "<<reserved_1<<"\n"
            <<"reserved_2: "<<reserved_2<<"\n"
            <<"offset: "<<offset<<"\n"
            <<"dib_header_size: "<<dib_header_size<<"\n"
            <<"width_pixels: "<<width_pixels<<"\n"
            <<"height_pixels: "<<height_pixels<<"\n"
            <<"colour_planes: "<<colour_planes<<"\n"
            <<"bits_per_pixel: "<<bits_per_pixel<<"\n"
            <<"compression_method: "<<compression_method<<"\n"
            <<"raw_image_size_bytes: "<<raw_image_size_bytes<<"\n"
            <<"horizontal_resolution: "<<horizontal_resolution<<"\n"
            <<"vertical_resolution: "<<vertical_resolution<<"\n"
            <<"num_colors: "<<num_colors<<"\n"
            <<"important_colors: "<<important_colors<<"\n";
    }
};

void updateHeader(BMPdata &bmp_data, unsigned int new_height_pixels,unsigned int new_width_pixels,unsigned int ***new_pixel_array)
{
    // updating heaeder fields
    unsigned int header_size = bmp_data.size_bytes-bmp_data.raw_image_size_bytes;

    bmp_data.height_pixels = new_height_pixels;
    bmp_data.width_pixels = new_width_pixels;
    bmp_data.raw_image_size_bytes = new_width_pixels*new_height_pixels*(bmp_data.bits_per_pixel/8);
    bmp_data.size_bytes = bmp_data.raw_image_size_bytes + header_size;
    bmp_data.pixel_array = new_pixel_array;
}

void scaleImage(BMPdata &bmp_data, double fx, double fy)
{
    unsigned new_height_pixels = bmp_data.height_pixels*fx;
    unsigned new_width_pixels = bmp_data.width_pixels*fy;
    new_width_pixels = (new_width_pixels+3)/4*4;

    unsigned int ***new_pixel_array = new unsigned int**[new_height_pixels];
    for(int i=0;i<new_height_pixels;i++)
    {
        new_pixel_array[i] = new unsigned int*[new_width_pixels];
        for(int j=0;j<new_width_pixels;j++)
        {
            new_pixel_array[i][j] = new unsigned int[3];
            int nearest_i = min((int)round(i*1.0/fx), (int)bmp_data.height_pixels-1);
            int nearest_j = min((int)round(j*1.0/fy), (int)bmp_data.width_pixels-1);
            for(int k=0;k<3;k++)
            {
                new_pixel_array[i][j][k] = bmp_data.pixel_array[nearest_i][nearest_j][k];
            }
        }
    }

    updateHeader(bmp_data,new_height_pixels,new_width_pixels,new_pixel_array);
}

bool withinRange(int i1,int j1,int i2,int j2,int n,int m)
{
  if(i1>=0 && i1<n && j1>=0 && j1<m && i2>=0 && i2<n && j2>=0 && j2<m) return true;
  return false;
}

void rotateImage(BMPdata &bmp_data,double angle,int type)
{
    angle=angle*PI/180;
    unsigned new_height_pixels = bmp_data.height_pixels*cos(abs(angle))+bmp_data.width_pixels*sin(abs(angle));
    unsigned new_width_pixels = bmp_data.height_pixels*sin(abs(angle))+bmp_data.width_pixels*cos(abs(angle));
    new_width_pixels = (new_width_pixels+3)/4*4;

    unsigned int ***int_pixel_array = new unsigned int**[new_height_pixels];
    for(int i=0;i<new_height_pixels;i++)
    {
        int_pixel_array[i] = new unsigned int*[new_width_pixels];
        for(int j=0;j<new_width_pixels;j++)
        {
            int_pixel_array[i][j] = new unsigned int[3];
            for(int k=0;k<3;k++)
            {
                if(i>=(new_height_pixels-bmp_data.height_pixels)/2 && i<bmp_data.height_pixels+(new_height_pixels-bmp_data.height_pixels)/2 && j>=(new_width_pixels-bmp_data.width_pixels)/2 && j<bmp_data.width_pixels+(new_width_pixels-bmp_data.width_pixels)/2)
                {
                    int_pixel_array[i][j][k] = bmp_data.pixel_array[i-(new_height_pixels-bmp_data.height_pixels)/2][j-(new_width_pixels-bmp_data.width_pixels)/2][k];
                }
                else
                {
                    int_pixel_array[i][j][k]=255;
                }
            }
        }
    }


    for(int i=0;i<new_height_pixels;i++)
    {
        for(int j=0;j<new_width_pixels;j++)
            cout<<int_pixel_array[i][j][2]<<" ";
        cout<<"\n";
    }

    angle = -angle; // since we do backward transform
    
    // Rotation of Image about the Image Center
    double origin_x=new_width_pixels/2;
    double origin_y=new_height_pixels/2;
    origin_y=-origin_y;


    unsigned int ***new_pixel_array = new unsigned int**[new_height_pixels];
    for(int i=0;i<new_height_pixels;i++)
    {
        new_pixel_array[i] = new unsigned int*[new_width_pixels];
        for(int j=0;j<new_width_pixels;j++)
        {
            new_pixel_array[i][j] = new unsigned int[3];

            // Inverse Rotation to get the coordinates of final to initial Image
            double cartesian_x=j;
            double cartesian_y=-i;
            double cartesian_newx = (cartesian_x-origin_x)*cos(angle)-(cartesian_y-origin_y)*sin(angle);
            double cartesian_newy = (cartesian_y-origin_y)*cos(angle)+(cartesian_x-origin_x)*sin(angle);

            double exact_i=-(cartesian_newy+origin_y);
            double exact_j=cartesian_newx+origin_x;
            
            for(int k=0;k<3;k++)
            {
                if(type==1)
                {
                    int i_=exact_i;
                    int j_=exact_j;
                    if(withinRange(i_,j_,i_+1,j_+1,new_height_pixels,new_width_pixels)==false)
                    {
                        new_pixel_array[i][j][k]=255;
                        continue;
                    }

                    // Bilinear Interpolation
                    double t1=int_pixel_array[i_][j_][k];
                    double t2=int_pixel_array[i_+1][j_][k];
                    double v1=(t1*(i_+1-exact_i)+t2*(exact_i-i_));
                    t1=int_pixel_array[i_][j_+1][k];
                    t2=int_pixel_array[i_+1][j_+1][k];
                    double v2=t1*(i_+1-exact_i)+t2*(exact_i-i_);

                    double res=v1*(j_+1-exact_j)+v2*(exact_j-j_);
                    new_pixel_array[i][j][k]=int(res);
                }
                else if(type==2)
                {
                    // Nearest Neighbour Interpolation
                    int i_=round(exact_i);
                    int j_=round(exact_j);
                    if(i_>=0 && i_<new_height_pixels && j_>=0 && j_<new_width_pixels)
                    {
                        new_pixel_array[i][j][k] = int_pixel_array[i_][j_][k];
                    }
                    else
                    {
                        new_pixel_array[i][j][k]=255;
                    }
                }
            }

        }
    }
    
    updateHeader(bmp_data,new_height_pixels,new_width_pixels,new_pixel_array);
}

int main()
{
    char filename[100];
    cout<<"Enter name of file to be read (*.bmp): ";
    cin>>filename;
    FILE* fptr = fopen(filename, "rb");
    strncpy(filename,filename,strlen(filename)-5);
    filename[strlen(filename)-4]='\0';
    
    BMPdata bmp_data;
    bmp_data.readBMP(fptr);
    bmp_data.printHeader();
    fclose(fptr);

    int op;
    cout<<"1. Scale\n2. Rotate\nEnter operation to be done: ";
    cin>>op;
    if(op==1)
    {
        double fx,fy;
        cout<<"Enter scale factor along height: ";
        cin>>fx;
        cout<<"Enter scale factor along width: ";
        cin>>fy;
        scaleImage(bmp_data, fx, fy);
    }
    else
    {
        double angle;
        int type;
        cout<<"Enter angle (degree) by which to rotate: ";
        cin>>angle;
        cout<<"1. Bilinear\n2. Nearest neighbour\nEnter type of interpolation: ";
        cin>>type;
        rotateImage(bmp_data,angle,type);
    }

    fptr = fopen(strcat(filename,"_out.bmp"), "w");
    bmp_data.writeBMP(fptr);
    fclose(fptr);
    return 0;
}