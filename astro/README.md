# Important Notes for the Astro Project

There are 2 data files that correlate to this project: astro64.vtk and astro512.vtk. The latter file is too large to push to GitHub so to access it download [this link](https://ix.cs.uoregon.edu/~hank/410/volren/astro512.vtk) and put it in the same directory as the source code for this project.

Even with optimizations, my code still takes at least 30 minutes to generate a high resolution image; this means the code uses 1000 x 1000 image dimensions, 1024 samples per ray, an early ray termination of .999, and the astro512.vtk data file. If you just want a sample of what this program does, feel free to set HIGH_RES to 0 and modify the dimensions, number of samples per ray, and early ray termination criteria at the top of the file.
