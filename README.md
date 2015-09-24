Objective: Create 16 threads to work in parallel to compute the two–dimensional DFT.

1. Prepare W_n.
2. Prepare transpose function.
3. Write Transform1D function using the Danielson–Lanczos approach and test it.
4. Use a barrier to insure all threads have completed the rows.
5. Transpose
6. Implement Transform2D function and test it. (column calculations with pthreads)
7. Wait (using a condition variable) until all threads are done.
8. Transpose
9. Save to file
10. Make sure it is able to take in customized inputs.