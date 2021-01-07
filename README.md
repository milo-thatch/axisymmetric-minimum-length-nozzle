# axisymmetric-minimum-length-nozzle
Computes the contour of the divergent section of an axisymmetric minimum length nozzle for a perfect gas.

This function computes the divergent section of an axisymmetric minimum length nozzle for a perfect gas. The kernel and transition regions are computed using the axisymmetric method of the characteristics. A normal sonic line is assumed at the throat.
The function requires as input:
- the exit Mach number (>1, supersonic);
- the heat capacity ratio (gamma=cp/cv);
- the number of characteristics (integer >0) in the linear kernel region;
- the number of characteristics (integer >0) in the compressed region downstream of the throat (x=0). When equal to one, the compression region is not resolved;
- the exponent of the power-law that defines the characteristics in the compressed region;
- maximum tolerance (residuals)
- aspect ratio of the mesh in the transition region (usually close to 1)
The function outputs a non-dimensional plot. Lengths are normalized to the throat radius. The code stores the coordinates of the contour points in a .txt file.

The code is based on the following references:
[1] J.D. Anderson. Modern Compressible Flow: With Historical Perspective. Aero-
nautical and Aerospace Engineering Series. McGraw-Hill Education, 2003 (Library of Congress CN: 2002067852)
[2] B. M. Argrow and G. Emanuel. Comparison of Minimum Length Nozzles. Journal of Fluids Engineering, 110(3):283-288, 09 1988.
[3] Kuno Foelsch. The analytical design of an axially symmetric Laval nozzle for a
parallel and uniform jet. Journal of the Aeronautical Sciences, 16(3):161-166,
1949.
[4] Ying-Nien Yu. A summary of design techniques for axisymmetric hypersonic
wind tunnels. Technical report, NATO - Science and Technology Organization,
1958.
