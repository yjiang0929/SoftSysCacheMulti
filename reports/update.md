## CacheMulti: Optimize matrix multiplication in C

### **Team members: Khang Vu, Yichen Jiang**

### **Project Goals:**

The MVP of this project is to finish the General Matrix Multiply (GEMM) optimization tutorials, understand how to optimize matrix multiplication in code, and write matrix multiplication code that is more optimized than conventional methods. We want to graph performance metrics to show that the new method is faster than the conventional one. Our medium range goal is to apply a divide and conquer algorithm with parallelization such as Strassen algorithm. The stretch goal involves understanding and utilizing a mesh network to optimize matrix multiplication.

### **Personal Learning Goals:**

- Khang: I am interested in learning about different matrix multiplication algorithms and their tradeoffs in terms of time and space complexity. Also, I would like to be better at programming in C and learn some low-level specific topics such as memory parallelism and distributed computing for optimization.
- Yichen: My primary goal in this project is to become a better C programmer. I think because C is such a low-level language, learning C can help me understand more about other programming languages and how I can optimize my code for time and space efficiency. I am also very interested in linear algebra and it sounds fascinating to me to understand how the most fundamental operation in linear algebra, matrix multiplication, is implemented and optimized in C. I want to learn more about the intersection of math and programming through this project.

### **Resources:**

The following resources are referenced in the initial stage of the project. As the project continues, we will expand this list with more links to tutorials and papers.

- [Github tutorial: Efficient matrix multiplication](https://gist.github.com/nadavrot/5b35d44e8ba3dd718e595e40184d03f0)
- [Github tutorial: GEMM optimization tutorial](https://github.com/flame/how-to-optimize-gemm)
- [High Performance Scientific Computing with C on Udemy](https://www.udemy.com/high-performance-scientific-computing-with-c/)
- [Wikipedia: Strassen algorithm](https://en.wikipedia.org/wiki/Strassen_algorithm)
- [Matrix multiplication: Strassen's algorithm](https://stanford.edu/~rezab/classes/cme323/S16/notes/Lecture03/cme323_lec3.pdf)
- [Wikipedia: Coppersmith–Winograd algorithm](https://en.wikipedia.org/wiki/Coppersmith%E2%80%93Winograd_algorithm)
- [Anatomy of high-performance matrix multiplication](https://dl.acm.org/doi/10.1145/1356052.1356053)

### **Updates:**

**What we have done:**

We have finished the High Performance Scientific Computing with C course on Udemy and the two Github tutorials. The Udemy course didn’t seem super useful as it only explains general concepts of time complexity analysis. We were able to follow the Github tutorials and understand different optimization methods. Register blocking and matrix tiling are two optimization techniques that yield the biggest performance gain. Matrix tiling means to multiply an m x n block of the matrix each time instead of multiplying each element one by one, and register blocking means saving part of the matrix values directly in the register so that we can minimize fetching time from the memory.

After understanding the concepts, we have implemented a test pipeline in C ([compare_matrix_multi.c](https://github.com/yjiang0929/SoftSysCacheMulti/blob/master/compare_matrix_multi.c)) to test the performance of the optimized matrix multiplication code and a function to generate random matrices for testing ([utils.c](https://github.com/yjiang0929/SoftSysCacheMulti/blob/master/utils.c#L31)). The test pipeline will generate random matrices of different sizes, load the given matrix multiplication algorithm and record the time it takes to compute matrix multiplication for different sizes. The result will be output into a .csv file to make plotting easier in the future.

**What we will do next:**

We are close to achieving our MVP goals. For the next few weeks of the project, we will mostly be working on the following items:

- Create comparison plot for basic optimization methods [Yichen]
    - Definition of done: Generate the plot to compare results of different GEMM methods.
- Understand and implement Strassen algorithms [Khang]
    - Definition of done: Read [Matrix multiplication: Strassen’s algorithm](https://stanford.edu/~rezab/classes/cme323/S16/notes/Lecture03/cme323_lec3.pdf); Implement the algorithm and run it using the test pipeline.
- Implement Strassen algorithms using multithreading [Khang, Yichen]
    - Definition of done: Read Head First C chapter 12 on multithreading; Attempt to implement Strassen with multithreading
- Compare Strassen with GEMM methods [Yichen]
    - Definition of done: Create comparison plots of GEMM and Strassen methods and analyze the difference.

### **Links:**

- Trello board: [https://trello. com/b/kcS0taS8/softsyscachemulti](https://trello.com/b/kcS0taS8/softsyscachemulti) 
- Github repo: [https://github.com/yjiang0929/SoftSysCacheMulti](https://github.com/yjiang0929/SoftSysCacheMulti)
