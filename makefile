
flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none

hist_3D: ~/programs/functions/precision.f08 ~/programs/functions/functions.f08 ~/programs/functions/stringlib.f90 compute_avgForce.3D.hist.f08
	gfortran -c ~/programs/functions/precision.f08 ~/programs/functions/functions.f08 ~/programs/functions/stringlib.f90 compute_avgForce.3D.hist.f08 $(flags)
	gfortran precision.o functions.o stringlib.o compute_avgForce.3D.hist.o $(flags) -o compute_avgForce.3D.hist.x

ellipse_2D: compute_avgForce.2D.ellipse.f08 stringlib.f90
	gfortran compute_avgForce.2D.ellipse.f08 stringlib.f90 $(flags) -o compute_avgForce.2D.ellipse.x

sphere_2D: compute_avgForce.2D.sphere.f08 stringlib.f90
	gfortran compute_avgForce.2D.sphere.f08 stringlib.f90 $(flags) -o compute_avgForce.2D.sphere.x

hist_2D: compute_avgForce.2D.hist.f08 stringlib.f90
	gfortran compute_avgForce.2D.hist.f08 stringlib.f90 $(flags) -o compute_avgForce.2D.hist.x

sphere_1D: compute_avgForce.1D.sphere.f08 stringlib.f90
	gfortran compute_avgForce.1D.sphere.f08 stringlib.f90 $(flags) -o compute_avgForce.1D.sphere.x

ellipse_1D: compute_avgForce.1D.ellipse.f08 stringlib.f90
	gfortran compute_avgForce.1D.ellipse.f08 stringlib.f90 $(flags) -o compute_avgForce.1D.ellipse.x

crd: crd_list.f08 stringlib.f90
	gfortran crd_list.f08 stringlib.f90 $(flags) -o crd_list.x

test: test.f08
	gfortran test.f08 $(flags) -o test.x

clean : 
	rm -f *.o *.mod *.x
