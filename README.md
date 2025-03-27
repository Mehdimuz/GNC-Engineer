# IGRF:
This directory contains all of the source code for running the IGRF.  The directory coeffs/ contains all of the coefficients available.  files titled 'dgrf' indicate that they contain definite coefficients.  the igrf2020s does not contain coefficients, but rather contains secular coefficients.

The getigrfcoefs script needs to be updated to produce a valid igrfcoefs.mat file to extend past the year 2020, using the secular coefficients for extrapolation (rather than the standard linear interpolation between two epochs)
