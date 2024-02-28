# Analytic-ICP

**Analytic Iterive Closet Point (Analytic-ICP)** is a point-set registration algorithm, originally developed by Wei Feng et al.

This is a C++ solution for visual measurement, which uses Analytic-ICP to register a template to actual measured data for subsequent measurement. It consists of three projects:

1. Analytic-ICP DLL code
2. Surface measure DLL code
3. Test application code

This solution can only run on the Windows visual stdio 2013 platform for the time being. At present, we have only implemented 2D and 3D second-order Analytic-ICP. The EXE generated in the Test project can output the visualization results of the registration and analytic mapping.

In the Analytic-ICP library, we used the point-to-point ICP code written by Professor Andreas Geiger, which requires installing the boost library to support kd tree, and no other third-party libraries are needed. In addition, in the SurfaceMeasure library and Test application code, we used OpenCV for image processing, so the OpenCV library needs to be imported into these two projects.
