#include <iostream>
#include <cmath>
using namespace std;

const int N = 8; // Block size

void removeDC(double **a)
{
	for(int row=0;row<N;row++)
	{
		for(int col=0;col<N;col++)
			a[row][col] -= 128;
	}
}

void dct(double *a, double *temp, int len)
{
	if (len == 1)
		return;
	for (int i = 0; i < len/2; i++)
	{
		double x = a[i];
		double y = a[len - 1 - i];
		temp[i] = x + y;
		temp[i + len/2] = (x - y) / (cos((i + 0.5) * M_PI / len) * 2);
	}
	dct(temp, a, len/2);
	dct(&temp[len/2], a, len/2);
	for (int i = 0; i < len/2 - 1; i++)
	{
		a[i * 2] = temp[i];
		a[i * 2 + 1] = temp[i + len/2] + temp[i + len/2 + 1];
	}
	a[len - 2] = temp[len/2 - 1];
	a[len - 1] = temp[len - 1];
}

void idct(double *a, double *temp, int len)
{
	if (len == 1)
		return;
	temp[0] = a[0];
	temp[len/2] = a[1];
	for (size_t i = 1; i < len/2; i++) {
		temp[i] = a[i * 2];
		temp[i + len/2] = a[i * 2 - 1] + a[i * 2 + 1];
	}
	idct(temp, a, len/2);
	idct(&temp[len/2], a, len/2);
	for (size_t i = 0; i < len/2; i++) {
		double x = temp[i];
		double y = temp[i + len/2] / (std::cos((i + 0.5) * M_PI / len) * 2);
		a[i] = x + y;
		a[len - 1 - i] = x - y;
	}
}

void dct2(double **a,bool inverse)
{
	double *temp = (double*)malloc(N*sizeof(double));
	double scaling = sqrt(2.0/N);
	for(int col=0;col<N;col++)
	{
		double cur[N];
		for(int row=0;row<N;row++)
			cur[row]=a[row][col];
		if(!inverse)
			dct(cur,temp,N);
		else
		{
			cur[0]/=2;
			idct(cur,temp,N);
		}
		for(int row=0;row<N;row++)
			a[row][col]=cur[row]*scaling;
	}

	for(int row=0;row<N;row++)
	{
		double cur[N];
		for(int col=0;col<N;col++)
			cur[col]=a[row][col];
		if(!inverse)
			dct(cur,temp,N);
		else
		{
			cur[0]/=2;
			idct(cur,temp,N);
		}
		for(int col=0;col<N;col++)
			a[row][col]=cur[col]*scaling;
	}
}

void quantize(double **a,double **q,int **b)
{
	for(int row=0;row<N;row++)
	{
		for(int col=0;col<N;col++)
			b[row][col] = round(a[row][col]/q[row][col]);
	}
}

void zigzag_order(int **b,int *z)
{
	int ind=0;
	bool down=0;
	for(int sum=0;sum<=2*N-2;sum++)
	{
		int row,col;
		if(down)
		{
			row=0,col=sum-row;
			if(col>=N)
				col=N-1,row=sum-col;
		}
		else
		{
			col=0,row=sum-col;
			if(row>=N)
				row=N-1,col=sum-row;
		}
		while(row>=0 && row<N && col>=0 && col<N)
		{
			z[ind]=b[row][col];
			ind++;
			if(down)
				row++,col--;
			else
				row--,col++;
		}
		down^=1;
	}
}

int num_bits(int x)
{
	if(x<0)
		x=-x;
	int count = 0;
	while(x != 1)
		count++, x/=2;
	return count+1;
}
void run_length_encode(int *z, unsigned char *rle)
{
	unsigned int l=0,r=1;
	unsigned char *extend;
	int len=0;
	int lst=N*N;
	while(lst>=0 && z[lst]==0)
		lst--;
	while(r<lst)
	{
		if(z[r]==0)
		{
			if(r-l+1==16)
			{
				len+=2;
				extend = (unsigned char*)realloc(rle,len*sizeof(unsigned char));
				if(extend)
					rle = extend;
				else
				{
					cout<<"Could not allocate memory!!";
					exit(0);
				}
				rle[len-2] = (15<<4), rle[len-1] = 0;
				l=r;
			}
		}
		else
		{
			len+=2;
			extend = (unsigned char*)realloc(rle,len*sizeof(unsigned char));
			if(extend)
				rle = extend;
			else
			{
				cout<<"Could not allocate memory!!";
				exit(0);
			}
			rle[len-2] = ((r-l-1)<<4)+num_bits(z[r]), rle[len-1] = z[r];
			l=r;
		}
		r++;
	}
	len+=2;
	extend = (unsigned char*)realloc(rle,len*sizeof(unsigned char));
	if(extend)
		rle = extend;
	else
	{
		cout<<"Could not allocate memory!!";
		exit(0);
	}
	rle[len-2] = rle[len-1] = 0;

	for(int i=0;i<len;i+=2)
	{
		cout<<"("<<((unsigned int)rle[i]>>4)<<","<<(unsigned int)(rle[i]&0x000000ff)<<") "<<(unsigned int)rle[i+1]<<"\n";
	}
}
void display(double **a)
{
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			cout<<a[i][j]<<" ";			
		}
		cout<<"\n";
	}
	cout<<"\n";
}
void encode(double **a)
{
	// removeDC(a);
	// display(a);
	// dct2(a,0);
	// display(a);

	int **b = (int**)malloc(N*sizeof(int*));
	double **q = (double**)malloc(N*sizeof(double*));
	for(int i=0;i<N;i++)
	{
		b[i]=(int*)malloc(N*sizeof(int));
		q[i]=(double*)malloc(N*sizeof(double));
		for(int j=0;j<N;j++)
			q[i][j]=1;
	}
	quantize(a,q,b);
	int *z = (int*)malloc(N*N*sizeof(int));
	zigzag_order(b,z);
	for(int i=0;i<N*N;i++)
	{
		cout<<z[i]<<" ";
	}
	cout<<"\n";

	unsigned char* rle=(unsigned char *)malloc(sizeof(unsigned char));
	run_length_encode(z,rle);
}
int main()
{
	int n;
	cin>>n;
	double **a = (double**)malloc(n*sizeof(double*));
	for(int i=0;i<n;i++)
	{
		a[i]=(double*)malloc(n*sizeof(double));
		for(int j=0;j<n;j++)
			cin>>a[i][j];
	}
	encode(a);
}