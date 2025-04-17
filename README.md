# Smooth Adjustment

This repository provides the official C++ implementation of our paper "Structured Analytic Mappings for Point Set Registration".

It consists of three projects:

 * Analytic-ICP DLL code
 * SmoothAdjustment DLL code
 * Test application code

The SmoothAdjustment dynamic link library constitutes the primary framework of this project. It is designed to perform specific operations on the registered surfaces, including visualization, measurements, and more. Currently, it relies on Analytic-ICP, but it can also utilize (rely on) other registration algorithms.

# Prerequisites
 * Visual Studio 2013
 * Eigen Library
 * Boost Library
 * Opencv Library

This solution can only run on the Windows visual stdio 2013 platform for the time being. At present, we have only implemented 2D and 3D nth-order Analytic-ICP. The EXE generated in the Test project can output the visualization results of the registration and analytic mapping.

In Analytic-ICP, we utilize the point-to-point ICP code authored by Professor Andreas Geiger, necessitating the installation of the Boost library to support the KD-tree functionality. We have modified the original pure C program by integrating matrix operations through the Eigen library. This enhancement enables the program to handle the registration of large-scale point cloud models more efficiently. Additionally, within the SmoothAdjustment library and Test application code, we employ OpenCV for image processing, thereby requiring the inclusion of the OpenCV library in these two projects.

# Datasets
 * SHREC’19
 * 3Dscanrep
 * Others

All the data utilized in our experiment are contained within this package, encompassing both 2D and 3D point clouds.

# Quick Example
 * Download the Release directory
 * Download the opencv_world300.dll and place it within the Release folder
 * Edit the experiment.ini
 * Run the Test.exe
