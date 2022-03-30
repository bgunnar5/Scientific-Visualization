# Scientific Visualization

The work in this repository is from a class I took in college. To respect the professor's wishes and to discourage cheating in the professor's future classes, I will not be saying where this class was taken or what the course number was. YOU DO NOT HAVE PERMISSION TO COPY THE CODE IN THIS REPOSITORY.

This class was dedicated to learning about methods of scientific visualization using VTK.

## How to Run the Projects

Step 5 is operating system dependant and is labeled as such. Note: I don't currently have a Windows machine so these instructions may be vague/incorrect. Apologies in advance.

1. Clone this repository on your system. In each project folder there will be the source code, a data file, and a CMakeLists.txt file.

2. Visit the [VTK website](https://vtk.org/download/) and download the latest version of VTK for your operating system.

3. Visit the [CMake website](https://cmake.org/download/) and download the latest version of CMake for your operating system.

4. Change directories to the project you want to run.

5.(Mac) Open the CMakeLists.txt file. Find where it says "ADD PATH TO VTK BUILD HERE" and replace it with the absolute path to your VTK build.

5.(Windows) Open the CMakeLists.txt file and modify it to work on your machine. This may already work as is or you may have to make changes to the entire file.

6. Run _cmake CMakeLists.txt_ at the command line. This will create several cmake related files, a Makefile, and an executable that can't be executed yet.

7. Run _make_ at the command line to compile the program.

8. Execute the program by running _./\<project name>.app/Contents/MacOS/\<project name>_ or just type _./\<project name>_ and tab until completion. Some projects show the data visualization instantly but others might generate images you can look at.
