## **CacheMulti: Optimize matrix multiplication in C**

### **Team members: Khang Vu, Yichen Jiang**

### **Project Goals**

The primary objectives of this project are to learn and implement various matrix multiplication techniques. Our MVP is to finish the General Matrix Multiply (GEMM) optimization tutorials, understand how to optimize matrix multiplication in code, and write matrix multiplication code that is more optimized than conventional methods. We want to graph performance metrics to show that the new method is faster than the conventional one. Our medium range goal is to apply a divide and conquer algorithm with parallelization such as the Strassen’s algorithm. The stretch goal involves understanding and utilizing a mesh network to optimize matrix multiplication.

### **Personal Learning Goals**

- Khang: I am interested in learning about different matrix multiplication algorithms and their tradeoffs in terms of time and space complexity. Also, I would like to be better at programming in C and learn some low-level specific topics such as memory parallelism and distributed computing for optimization.
- Yichen: My primary goal in this project is to become a better C programmer. I think because C is such a low-level language, learning C can help me understand more about other programming languages and how I can optimize my code for time and space efficiency. I am also very interested in linear algebra and it sounds fascinating to me to understand how the most fundamental operation in linear algebra, matrix multiplication, is implemented and optimized in C. I want to learn more about the intersection of math and programming through this project.

### **Resources**

The following resources are referenced in the initial stage of the project. As the project continues, we will expand this list with more links to tutorials and papers.

- [Github tutorial: Efficient matrix multiplication](https://gist.github.com/nadavrot/5b35d44e8ba3dd718e595e40184d03f0)
- [Github tutorial: GEMM optimization tutorial](https://github.com/flame/how-to-optimize-gemm)
- [High Performance Scientific Computing with C on Udemy](https://www.udemy.com/high-performance-scientific-computing-with-c/)
- [Wikipedia: Strassen algorithm](https://en.wikipedia.org/wiki/Strassen_algorithm)
- [Matrix multiplication: Strassen's algorithm](https://stanford.edu/~rezab/classes/cme323/S16/notes/Lecture03/cme323_lec3.pdf)
- [Wikipedia: Coppersmith–Winograd algorithm](https://en.wikipedia.org/wiki/Coppersmith%E2%80%93Winograd_algorithm)
- [Anatomy of high-performance matrix multiplication](https://dl.acm.org/doi/10.1145/1356052.1356053)
- Head First C Chapter 12: Threads

### **What we did**

#### **Overview**

To familiarize ourselves with scientific computing, we first finished the High Performance Scientific Computing with C course on Udemy and two Github tutorials. The Udemy course didn’t seem super useful as it only explains general concepts of time complexity analysis. Follow the Github tutorials, we were able to understand different optimization methods for matrix multiplication. Register blocking and matrix tiling are two optimization techniques that yield the biggest performance gain. Matrix tiling means to multiply an m x n block of the matrix each time instead of multiplying each element one by one, and register blocking allows the program to save part of the matrix values directly in the register and minimizes fetching time from the memory.

After understanding the concepts, we implemented a test pipeline in C ([compare_matrix_multi.c](https://github.com/yjiang0929/SoftSysCacheMulti/blob/master/compare_matrix_multi.c)) to test the performance of different matrix multiplication methods explained in the [General Matrix Multiply (GEMM) optimization tutorial](https://github.com/flame/how-to-optimize-gemm). We also created a function to generate random matrices for testing ([utils.c](https://github.com/yjiang0929/SoftSysCacheMulti/blob/master/utils.c#L31)). The test pipeline will generate random matrices of different sizes, load the given matrix multiplication algorithm and record the time it takes to compute matrix multiplication for different sizes. The result will be output into a .csv file to make plotting easier.

We then implemented the Strassen’s algorithm described in the [Matrix multiplication: Strassen's algorithm](https://stanford.edu/~rezab/classes/cme323/S16/notes/Lecture03/cme323_lec3.pdf) paper. The performance of the Strassen’s algorithm is much slower than those of the GEMM optimization methods, even though it seems to improve more slightly as the matrix size gets bigger. The advantage of Strassen’s algorithm is that it allows us to divide a task into smaller sub-tasks, which can be handled in parallel. Therefore, we decided to learn about multi-threading by reading Chapter 12 of Head First C and implemented another Strassen version using 7 threads. The final result with parallelism is astounding and is better than the GEMM approaches when the matrix size gets bigger.

#### **GEMM methods explained**

GEMM stands for General Matrix Multiplication. For the first half of the project, we looked at GEMM implementations and compared the performance of 5 different methods of general matrix optimization following [the tutorial](https://github.com/flame/how-to-optimize-gemm). We created a summary document for all the optimization methods presented in the tutorial [here](https://docs.google.com/document/d/1EKFd6r8zo3Vo5-Weva1ecjFY3Y4hRef2UQacKBtBrSI/edit?usp=sharing). From all the different optimization methods, we chose 4 that generated a significant performance boost and their performance curves are shown in the results section below.

The most basic optimization in `MMult_1x4_reg.c` involves using registers for elements of matrices and computes 4 elements at the same time in a single subroutine. In code, this is achieved by creating `register double` to hold values to be multiplied and updating 4 elements of C in a single function call:
```C
for ( j=0; j<n; j+=4 ) { /* Loop over the columns of C, unrolled by 4 */
  for ( i=0; i<m; i+=1 ) { /* Loop over the rows of C */
    AddDot1x4( k, &A( i,0 ), lda, &B( 0,j ), ldb, &C( i,j), ldc );
  }
}
```
`MMult_4x4_vecreg.c` takes the idea one step further by computing 4x4 blocks of C at a time. This allows us to put data into 16 vector registers and use special SSE3 instructions to speed up multiply and accumulate operations. This greatly increases the speed for lower rank matrices but the performance drops significantly when the sizes of matrices become large.

To solve this problem, we separate the resulting matrix C into smaller sub-blocks and compute each block with a separate subroutine call in `MMult_4x4_vecreg_subblock.c`. This brings up the performance of bigger matrices because it makes sure we have enough cache space to hold all the information about the multiplied matrices and don’t need to go to memory to fetch additional information.

As a final step, in `MMult_4x4_vecreg_subblock_cache.c`, we cache kx4 blocks of A and B between the operations so we can avoid duplicate computation and data fetching. This further improves the performance by a large margin.

#### **Strassen’s algorithm**

The main idea behind Strassen’s algorithm is that we can treat matrix multiplication as a recursive problem. Strassen’s algorithm takes advantage of divide-and-conquer to recursively break large matrices into smaller sub-blocks that can be used to compute the products of two matrices. To simplify the implementation, we assumed input matrices are square whose dimensions are the power of two. Below are the algorithms of the conventional approach and the Strassen approach:
```
A               B               A*B (conventional approach)
 +-----------+   +-----------+   +-------------------+-------------------+
 | A11 | A12 |   | B11 | B12 |   | A11*B11 + A12*B21 | A11*B12 + A12*B22 |
 +-----+-----+ * +-----+-----+ = +-------------------+-------------------+
 | A21 | A22 |   | B21 | B22 |   | A21*B11 + A22*B21 | A21*B12 + A22*B22 |
 +-----+-----+   +-----+-----+   +---------------------------------------+


 Strassen approach: Seven products:
 P1 = (A11+A22)(B11+B22)
 P2 = (A21+A22)B11
 P3 = A11(B12−B22)
 P4 = A22(B21−B11)
 P5 = (A11+A12)B22
 P6 = (A21−A11)(B11+B12)
 P7 = (A12−A22)(B21+B22)

         +-------------+-------------+
         | P1+P4-P5+P7 |    P3+P5    |
 A * B = +-------------+-------------+
         |    P2+P4    | P1+P2-P3+P6 |
         +-------------+-------------+
```
As we can see, the conventional approach needs 8 multiplications on sub-blocks to compute the product, while the Strassen approach only uses 7 products. It’s worth noting that the Strassen’s algorithm requires a few more summations. However, the number of summations is a constant independent of the size of our matrices ([source](https://stanford.edu/~rezab/classes/cme323/S16/notes/Lecture03/cme323_lec3.pdf)). As a result, the time complexity of the Strassen’s algorithm is O(n) = nlog2, which is better than O(n) = n3 of the conventional approach.

#### **Design decision: Strassen with multi-threading**

After we completed the implementation for the Strassen’s algorithm, we found that the algorithm performed barely better than the basic implementation. In order to improve the performance of Strassen’s algorithm, we wanted to experiment using multi-threading on the seven recursive calls to parallelize the operation. Initially, our idea was to put each of seven recursive calls into a separate thread but this quickly ran into problems. As the program recursively ran the Strassen’s algorithm, it kept generating threads and ran out of memory quickly even for small-sized matrices. To deal with this problem, we decided to only break the Strassen’s algorithm into seven parallel threads in the first level of recursion and use the non-threaded algorithm for further recursive calls. This resolved the segmentation fault that we got in the previous step and brought the Strassen’s algorithm’s performance on the same level with the final step of GEMM optimization.

### **Results**

After completing implementations for GEMM optimization and Strassen’s algorithm in this project, we used the test bench to create a graph to compare the performance between different optimization methods. We tested matrix sizes ranging from 4 to 4096 to see how well each method worked for smaller or bigger matrices. The performance curves for 6 different optimizations along with the basic matrix multiplication method are shown in the graph below.

![Performance comparison](https://github.com/yjiang0929/SoftSysCacheMulti/blob/master/comparison.png)

Most of the GEMM methods performed as expected from the tutorial. The vector-specific instructions used in MMult_4x4_vecreg showed a slight improvement over the previous method, but combined with sub-blocking and caching, it achieved a consistent performance at almost 10 GFlops / sec. The Strassen’s algorithm by itself only performed marginally better than the basic matrix multiplication. We think this is because our implementation wasn’t optimized for caching and register uses, and in this case, hardware optimized code brings more performance benefits than an optimized algorithm. After we implemented multi-threading with the Strassen’s algorithm, its performance became on par with the final step of GEMM optimization.

### **Reflection**

- Khang Vu: In this project, I learned about how to write better C programs that can utilize caching. I also had a chance to learn about the Strassen’s algorithm to take advantage of recursion and divide-and-conquer. As a new programmer in C, I feel like I gained a lot of C programming experience that could be very useful for me to become a better programmer.
- Yichen Jiang: I really enjoyed working on this project, because it gave me a chance to review concepts from linear algebra and combine them with new programming knowledge about C. I was able to learn more about how caching and threading works in C, which will be useful if I need to write high performance computation code in the future.

### **Links**

- Trello board: [https://trello. com/b/kcS0taS8/softsyscachemulti](https://trello.com/b/kcS0taS8/softsyscachemulti) 
- Github repo: [https://github.com/yjiang0929/SoftSysCacheMulti](https://github.com/yjiang0929/SoftSysCacheMulti)
