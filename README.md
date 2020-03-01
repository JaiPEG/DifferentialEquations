This project includes routines to solve differential equations, namely the
scalar wave equation in 1 + 1 dimensions.

It is submitted as the second assignment for the 2020 Computational Physics
course at Perimeter Institute.

## Animating a wave function

An example program included in this project that generates plots of an evolving wave function and saves them as SVG files under output/animWave/.

In order to run this program, ensure that the Julia package Gadfly is installed on your system, set your working directory to the project directory, and then type at the commandline,

    julia src/animWave.jl

If you have ImageMagick and FFmpeg installed on your system, you can use the following script to generate a video file from these plots, saved to output/animWave/vid.mp4.

    ./vid.sh

An example plot and video file is included in the project repository.

![alt text](https://github.com/JaiPEG/DifferentialEquations/raw/master/output/animWave/frame-0240.png "Frame 240 of wave animation")

## Plotting the energy over time at various resolutions

An example program included in this project generates a plot of the energy over time of evolving wave functions at several (spacial) resolutions; the time-step is also decreased in order to have good conditioning. The plot is saved as an SVG file under output/plotEnergy/.

In order to run this program, ensure that the Julia package Gadfly is installed on your system, set your working directory to the project directory, and then type at the commandline,

    julia src/plotEnergy.jl

An example plot is included in the project repository. Notice that as the number of samples increase (the resolution decreases), the energy changes less and less over the time interval, demonstrating convergence.

![alt text](https://github.com/JaiPEG/DifferentialEquations/raw/master/output/plotEnergy/energy.png "Energy change over time at various resolutions")
