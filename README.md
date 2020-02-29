This project includes routines to solve differential equations, namely the
scalar wave equation in 1 + 1 dimensions.

It is submitted as the second assignment for the 2020 Computational Physics
course at Perimeter Institute.

## Animating a wave function

An example program included in this project generates plots of an evolving wave function and saves them as SVG files under output/animWave/. In order to run this program, ensure that the Julia package Gadfly is installed on your system, set your working directory to the project directory, and then type at the commandline,

    julia --project=@. src/animWave.jl

If you have ImageMagick and FFmpeg installed on your system, you can use the following script to generate a video file from these plots, saved to output/animWave/vid.mp4.

    ./vid.sh

An example video file is included in the project repository.
