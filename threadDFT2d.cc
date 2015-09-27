// Threaded two-dimensional Discrete FFT transform
// Yushan Cai
// ECE8893 Project 2


#include <iostream>
#include <string>
#include <cstring>
#include <math.h>

#include "Complex.h"
#include "InputImage.h"

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.

using namespace std;

int N = 1024; // the number of points in the 1D transform.
// Function to reverse bits in an unsigned integer
unsigned ReverseBits(unsigned v)
{
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

//transpose function for 1d array
void Transpose(Complex* src, int& height, int& width)
{
  Complex* temp = new Complex[height*width];
  memcpy(temp, src, sizeof(Complex) * height * width);
  int r = 0, c = 0, index = 0;
  for(int j = 0; j < height * width; j++)
  {
    r = j / width;
    c = j % width;
    index = c * height + r;
    src[index] = temp[j];
  }
  delete [] temp;
  return;
}

void testReverse()
{
  unsigned index[N];
  for (int i = 0; i < N; ++i)
  {
    index[i] = ReverseBits(i);
    cout << index[i] << ' ';
  }
  cout << endl;
}

// Danielson-Lanczos DFT, where "h" is an input/output parameter
// "N" is the size of the array (assume even power of 2)
void Transform1D(Complex* h, int N)
{
  if (N == 1)
  {
    return;
  }
  Complex * tmp = new Complex[N];
  int index = 0;
  //reindexing array into tmp
  //cout << "\nReformed order of array" << endl;
  for (int i = 0; i < N; ++i)
  {
    index = ReverseBits(i);
    memcpy(&tmp[i], &h[index], sizeof(Complex));
    //cout << index << ' ';
  }
  //cout << endl;
  //Precompute W_n, since W_n[n+N/2] = -W_n[n], only half is computed
  int n = N / 2;
  Complex * W_n = new Complex[n];
  for (int i = 0; i < n; ++i)
  {
    W_n[i].real = cos(double(2 * M_PI * i / N));
    W_n[i].imag = - sin(double(2 * M_PI * i / N));
  }
  int step = 2;
  while (step <= N)
  {
    for (int i = 0; i < N/step; i += step)
    {
      for (int j = 0; j < step/2; ++j)
      {
        Complex buf = tmp[i*step+j];
        Complex W = W_n[j*N/step];
        tmp[i*step+j] = tmp[i*step+j] + W * tmp[i*step+j+step/2];
        tmp[i*step+j+step/2] = buf - W * tmp[i*step+j+step/2];
      }
    }
    //cout << step << "finished" << endl;
    step *= 2;
  }
  memcpy(h, tmp, N*sizeof(Complex));
}

void test1D(Complex* h, int N)
{
  cout << "\nAfter transform" << endl;
  Transform1D(h, N);
  for (int i = 0; i < N; ++i)
  {
    cout << h[i].Mag() << ' ';
  }
  cout << endl;
}

void* Transform2DTHread(void* v)
{ // This is the thread starting point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete
  return 0;
}

void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.
  InputImage image(inputFN);  // Create the helper object for reading the image
  int width = image.GetWidth();
  int height = image.GetHeight();
  Complex * data = image.GetImageData();
  Complex * test = new Complex[width];
  // for (int i = 0; i < height; ++i)
  // {
  // cout<< "before transform" << endl;
  //     for (int j = 0; j < width; ++j)
  //   {
  //     memcpy(&test[j], &data[j], sizeof(Complex));
  //     cout << test[j].Mag() << " ";
  //   }
  //   cout << endl;
  //   test1D(test , width);
  // }
  // Create the global pointer to the image array data
  // Create 16 threads
  pthread_t threads[16];
  // Wait for all threads complete
  // Write the transformed data
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  Transform2D(fn.c_str()); // Perform the transform.
  //testReverse();
}