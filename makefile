
dir = ~/programs
func = $(dir)/functions

flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none

test: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 test.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 test.f90 $(flags)
	gfortran precision.o functions.o constants.o stringlib.o test.o $(flags) -o test.x

spline_3D: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.3D.hist.spline.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.3D.hist.spline.f90 $(flags)
	gfortran precision.o functions.o constants.o ideal_CL3.o stringlib.o compute_avgForce.3D.hist.spline.o $(flags) -o compute_avgForce.3D.hist.spline.x

spline_1D_parallel: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.1D.hist.spline.parallel.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.1D.hist.spline.parallel.f90 $(flags)
	gfortran precision.o functions.o constants.o ideal_CL3.o stringlib.o compute_avgForce.1D.hist.spline.parallel.o $(flags) -o compute_avgForce.1D.hist.spline.parallel.x

spline_2D_parallel: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.2D.hist.spline.parallel.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.2D.hist.spline.parallel.f90 $(flags)
	gfortran precision.o functions.o constants.o ideal_CL3.o stringlib.o compute_avgForce.2D.hist.spline.parallel.o $(flags) -o compute_avgForce.2D.hist.spline.parallel.x

spline_3D_parallel: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.3D.hist.spline.parallel.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.3D.hist.spline.parallel.f90 $(flags)
	gfortran precision.o functions.o constants.o ideal_CL3.o stringlib.o compute_avgForce.3D.hist.spline.parallel.o $(flags) -o compute_avgForce.3D.hist.spline.parallel.x

hist_3D: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 compute_avgForce.3D.hist.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 compute_avgForce.3D.hist.f90 $(flags)
	gfortran precision.o functions.o constants.o stringlib.o compute_avgForce.3D.hist.o $(flags) -o compute_avgForce.3D.hist.x

ellipse_2D: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 compute_avgForce.2D.ellipse.f90 $(func)/stringlib.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 compute_avgForce.2D.ellipse.f90 $(flags)
	gfortran precision.o functions.o constants.o stringlib.o compute_avgForce.2D.ellipse.o $(flags) -o compute_avgForce.2D.ellipse.x

sphere_2D: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 compute_avgForce.2D.sphere.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 compute_avgForce.2D.sphere.f90 $(flags)
	gfortran precision.o functions.o constants.o stringlib.o compute_avgForce.2D.sphere.o $(flags) -o compute_avgForce.2D.sphere.x

hist_2D: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 compute_avgForce.2D.hist.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 compute_avgForce.2D.hist.f90 $(flags)
	gfortran precision.o functions.o constants.o stringlib.o compute_avgForce.2D.hist.o $(flags) -o compute_avgForce.2D.hist.x

sphere_1D: compute_avgForce.1D.sphere.f90 $(func)/stringlib.f90
	gfortran -c compute_avgForce.1D.sphere.f90 $(func)/stringlib.f90 $(flags)
	gfortran compute_avgForce.1D.sphere.o stringlib.o $(flags) -o compute_avgForce.1D.sphere.x

ellipse_1D: compute_avgForce.1D.ellipse.f90 $(func)/stringlib.f90
	gfortran -c compute_avgForce.1D.ellipse.f90 $(func)/stringlib.f90 $(flags)
	gfortran compute_avgForce.1D.ellipse.o stringlib.o $(flags) -o compute_avgForce.1D.ellipse.x

crd: crd_list.f90 $(func)/stringlib.f90
	gfortran -c crd_list.f90 $(func)/stringlib.f90 $(flags)
	gfortran crd_list.o stringlib.o $(flags) -o crd_list.x

clean : 
	rm -f *.o *.mod *.x
