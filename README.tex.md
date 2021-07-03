# Mandelbrot-Set

The Mandelbrot set is defined as the set of numbers in the $\mathcal{R}^2$ that satisfies the following quadratic recurrence equation:

$$
Z_{n+1} = Z_n + C
$$

with $z_0=C$, where C are points in the complex plane for which the values of $Z_n$ does not tend to infinity. Depending on the divergence rate of each point, that is to say, the number of iterations it takes to stabilize, the final plot colors illustrates the points inside the Mandelbrot set as well as its divergence rate.

The aim of this project is to implement a model for computing the Mandelbrot set. The model at issue must be coded in C language and executed in parallel.

Two different versions of the code have been implemented. In the first case (\textit{Version 1}), the code does not use a global map since it only returns the coordinates of the points that belong to the Mandelbrot set. Moreover, this code is implemented in a single dimension.

In the other case (\textit{Version 2}), it is necessary to built a global map, where each position of the map is associated to a coordinate X and a coordinate Y, and it contains a value from 0 to 1 which indicates whether the point at issue belongs to the set or not (and other information that will be described in the below sections). 

## Version 1
The code computes the Mandelbrot set for a range of $x \in [-2, 1.5]$ and $y \in [-1.5 1.5]$. The below results are obtained by fixing a maximum number of iterations of 250 iterations and a number of divisions of $1E4 \times 1E4$. 
<img src="img/Mandelbrot_1D.png" width="750">


## Version 2
The below results are obtained by fixing a maximum number of iterations of 100 iterations and a number of divisions of $1E4 \times 1E4$. 
<img src="img/Mandelbrot_BW.png" width="750">
<img src="img/Mandelbrot_JET.png" width="750">
