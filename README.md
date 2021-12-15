OpenFOAM functionObject for calculating the viscous stress tensor field using the postProcess utility.

To calculate stress vector in ParaView:

1) Calculate surface normals using "Generate surface normals" filter
2) Using the Calculator, calculate the stress vector as:

(Normals_X*viscousStressTensor_XX + Normals_Y*viscousStressTensor_XY + Normals_Z*viscousStressTensor_XZ)*iHat + (Normals_X*viscousStressTensor_XY + Normals_Y*viscousStressTensor_YY + Normals_Z*viscousStressTensor_YZ)*jHat + (Normals_X*viscousStressTensor_XZ + Normals_Y*viscousStressTensor_YZ + Normals_Z*viscousStressTensor_ZZ)*kHat