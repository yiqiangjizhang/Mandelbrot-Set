# Mandelbrot-Set

The Mandelbrot set is defined as the set of numbers in the <!-- $\mathcal{R}^2$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathcal%7BR%7D%5E2"> that satisfies the following quadratic recurrence equation:

<!-- $$
Z_{ n+1 } = Z_n + C
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=Z_%7B%20n%2B1%20%7D%20%3D%20Z_n%20%2B%20C"></div>

with <!-- $z_0=C$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=z_0%3DC">, where <!-- $C$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=C"> are points in the complex plane for which the values of <!-- $Z_n$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=Z_n"> does not tend to infinity. Depending on the divergence rate of each point, that is to say, the number of iterations it takes to stabilize, the final plot colors illustrates the points inside the Mandelbrot set as well as its divergence rate.

The aim of this project is to implement a model for computing the Mandelbrot set. The model at issue must be coded in C language and executed in parallel.

Two different versions of the code have been implemented. In the first case (Version 1), the code does not use a global map since it only returns the coordinates of the points that belong to the Mandelbrot set. Moreover, this code is implemented in a single dimension.

In the other case (Version 2), it is necessary to built a global map, where each position of the map is associated to a coordinate <!-- $X$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=X"> and a coordinate <!-- $Y$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=Y">, and it contains a value from 0 to 1 which indicates whether the point at issue belongs to the set or not (and other information that will be described in the below sections).


## Version 1
The code computes the Mandelbrot set for a range of <!-- $x \in [-2, 1.5]$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=x%20%5Cin%20%5B-2%2C%201.5%5D"> and <!-- $y \in [-1.5 1.5]$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=y%20%5Cin%20%5B-1.5%201.5%5D">. The below results are obtained by fixing a maximum number of iterations of 250 iterations and a number of divisions of <!-- $1E4 \times 1E4$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=1E4%20%5Ctimes%201E4">.

<img src="img/Mandelbrot_1D.png" width="750">


## Version 2

The below results are obtained by fixing a maximum number of iterations of 100 iterations and a number of divisions of <!-- $1E4 \times 1E4$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=1E4%20%5Ctimes%201E4">.

<img src="img/Mandelbrot_BW.png" width="750">
<img src="img/Mandelbrot_JET.png" width="750">
