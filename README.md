# Smooth Adjustment

**Analytic Iterive Closest Point (Analytic-ICP)** is a point-set registration algorithm, originally developed by Wei Feng et al.

This is a C++ solution for surface registration (adjustment), which uses Analytic-ICP to register moving feature point set in one surface to fixed feature point set in another surface. It consists of three projects:

1. Analytic-ICP DLL code
2. SmoothAdjustment DLL code
3. Test application code

The SmoothAdjustment dynamic link library constitutes the primary framework of this project. It is designed to perform specific operations on the registered surfaces, including adjustments, measurements, and more. Currently, it relies on Analytic-ICP, but it can also utilize (rely on) other registration algorithms.

# Prerequisites
 * Visual Studio 2013
 * Eigen Library
 * Boost Library
 * Opencv Library

This solution can only run on the Windows visual stdio 2013 platform for the time being. At present, we have only implemented 2D and 3D nth-order Analytic-ICP. The EXE generated in the Test project can output the visualization results of the registration and analytic mapping.

In Analytic-ICP, we utilize the point-to-point ICP code authored by Professor Andreas Geiger, necessitating the installation of the Boost library to support the KD-tree functionality. We have migrated it to the Eigen version to accommodate large point cloud models. Additionally, within the SmoothAdjustment library and Test application code, we employ OpenCV for image processing, thereby requiring the inclusion of the OpenCV library in these two projects.
