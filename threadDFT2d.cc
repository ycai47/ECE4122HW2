// Threaded two-dimensional Discrete FFT transform
// Yushan Cai
// ECE8893 Project 2

#include <iostream>
#include <string>
#include <cstring>
#include <math.h>
#include <pthread.h>
#include <cstdlib>
#include <stdint.h>

#include "Complex.h"
#include "InputImage.h"

using namespace std;

//Gobal variables
pthread_mutex_t imageMutex;
pthread_mutex_t exitMutex;
pthread_mutex_t startCountMutex;
pthread_cond_t exitCond;
pthread_barrier_t bar;
int startCount;
Complex* ImageData;
int ImageWidth;
int ImageHeight;
int nThreads;


//transpose function for 1d array
void Transpose(Complex * src, int & width, int & height)
{
  Complex * temp = new Complex[height * width];
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

// Function to reverse bits in an unsigned integer
unsigned ReverseBits(unsigned v, int N)
{
  bool check = (N > 0) && !(N & (N-1));
  if (check == 0)
  {
    cout << "Unacceptable file dimension. Width should be an int which is power of 2." << endl;
    return 0;
  }
  else
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
}

void testReverse()
{
  int length = 1024;
  unsigned index[length];
  for (int i = 0; i < length; ++i)
  {
    index[i] = ReverseBits(i, length);
    cout << index[i] << ' ';
  }
  cout << endl;
}

// Danielson-Lanczos DFT, where "h" is an input/output parameter
// "N" is the size of the array (assume even power of 2)
void Transform1D(Complex* h, int N)
{
  //Make sure size is reasonable to do the transform
  if (N == 1)
  {
    return;
  }

  //reindexing array into tmp
  Complex * tmp = new Complex[N];
  int index = 0;
  for (int i = 0; i < N; ++i)
  {
    index = ReverseBits(i, ImageWidth);
    memcpy(&tmp[i], &h[index], sizeof(Complex));
  }

  //Precompute W_n, since W_n[n+N/2] = -W_n[n], only half is computed
  int n = N / 2;
  Complex * W_n = new Complex[n];
  for (int i = 0; i < n; ++i)
  {
    W_n[i].real = cos(double(2 * M_PI * i / N));
    W_n[i].imag = - sin(double(2 * M_PI * i / N));
    if (fabs(W_n[i].real) < 1E-10) W_n[i].real = 0;
    if (fabs(W_n[i].imag) < 1E-10) W_n[i].imag = 0; //adjust double precision
  }

  //Actual butterfly computation
  int step = 2;
  while (step <= N)
  {
    for (int i = 0; i < N / step; ++i)
    {
      for (int j = 0; j < step / 2; ++j)
      {
        Complex buf = tmp[i * step + j];
        Complex W = W_n[j * N / step];
        tmp[i * step + j] = tmp[i * step + j] + W * tmp[i * step + j + step / 2];
        tmp[i * step + j + step / 2] = buf - W * tmp[i * step + j + step / 2];
      }
    }
    step *= 2;
  }
  memcpy(h, tmp, N * sizeof(Complex));
  delete [] tmp;
}

void test1D(Complex * h, int N)
{
  for (int i = 0; i < ImageHeight; ++i)
  {
    Transform1D(&h[i * ImageWidth], N);
  }
  InputImage image("Tower.txt");
  image.SaveImageData("Tower-DFT1D.txt ",ImageData, ImageWidth, ImageHeight); 
}

// This is the thread starting point.  "v" is the thread number
void* Transform2DThread(void* v)
{
  unsigned long myId = (unsigned long) v;
  // Calculate 1d DFT for assigned rows
  int rowsPerThread = ImageHeight / nThreads;
  int startingRow = myId * rowsPerThread;
  for (int r = 0; r < rowsPerThread; ++r)
  {
    pthread_mutex_lock(&imageMutex);
    Transform1D(&ImageData[(startingRow + r) * ImageWidth], ImageWidth);
    // cout << "Transform on thread" << myId << "row" <<r << "done.\n";
    pthread_mutex_unlock(&imageMutex);
  }
  // wait for all to complete
  pthread_barrier_wait(&bar);
  //Let thread 0 to transpose the overall array to prepare for column transformation
  pthread_mutex_lock(&imageMutex);
  if(myId == 0)
  {
    Transpose(ImageData, ImageWidth, ImageHeight);
  }
  pthread_mutex_unlock(&imageMutex);
  pthread_barrier_wait(&bar);
  // Calculate 1d DFT for assigned columns
  int colsPerThread = ImageWidth / nThreads;
  int startingCol = myId * colsPerThread;
  for (int c = 0; c < colsPerThread; ++c)
  {
    pthread_mutex_lock(&imageMutex);
    Transform1D(&ImageData[(startingCol + c) * ImageHeight], ImageHeight);
    pthread_mutex_unlock(&imageMutex);
  }
  //Decrement active count and signal TRansform2D if all complete
  pthread_mutex_lock(&startCountMutex);
  startCount--;
  if (startCount == 0)
  {
    Transpose(ImageData, ImageWidth, ImageHeight);
    pthread_mutex_unlock(&startCountMutex);
    pthread_mutex_lock(&exitMutex);
    pthread_cond_signal(&exitCond);
    pthread_mutex_unlock(&exitMutex);
  }
  else
  {
    pthread_mutex_unlock(&startCountMutex);
  }
  return 0;
}

void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.
  InputImage image(inputFN);  // Create the helper object for reading the image
  ImageWidth = image.GetWidth();
  ImageHeight = image.GetHeight();
  ImageData = image.GetImageData();
  // test1D(ImageData , ImageWidth);

  //Initializations
  pthread_mutex_init(&startCountMutex, 0);
  pthread_mutex_init(&imageMutex, 0);
  pthread_mutex_init(&exitMutex, 0);
  pthread_cond_init(&exitCond, 0);
  pthread_barrier_init(&bar, NULL, nThreads);
  startCount = nThreads;

  // Create 16 threads
  int rc;
  for (int i = 0; i < nThreads; ++i)
  {
    pthread_t pt;
    rc = pthread_create(&pt, 0, Transform2DThread, (void*)(uintptr_t)i);
    if(rc)
    {
      cout<<"Fail to create thread. rc: " << rc<< endl;
      exit(1);
    }
  }
  // Wait for all threads complete
  pthread_cond_wait(&exitCond, &exitMutex);
  // Write the transformed data
  image.SaveImageData("Tower-DFT2D.txt ", ImageData, ImageWidth, ImageHeight); 
}

int main(int argc, char** argv)
{
  string fn("Tower.txt"); // default file name
  nThreads = 16; // default thread numbers
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  if (argc > 2) nThreads = atol(argv[2]); // if thread number specified on cmd line
  Transform2D(fn.c_str()); // Perform the transform.
  //testReverse();
}