# Analytic-ICP

**Analytic Iterive Closet Point (Analytic-ICP)** is a point-set registration algorithm, originally developed by Wei Feng et al.

This is a C++ solution for visual measurement, which uses Analytic-ICP to register a template to actual measured data for subsequent measurement. It consists of three projects:

1. Analytic-ICP DLL code
2. Surface measure DLL code
3. Test application code

This solution can only run on the Windows visual stdio 2013 platform for the time being. At present, we have only implemented 2D and 3D second-order Analytic-ICP. The EXE generated in the Text project can output the visualization results of the registration and analytic mapping.
