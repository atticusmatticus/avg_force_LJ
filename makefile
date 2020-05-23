
dir = ~/programs
func = $(dir)/functions
#flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none -Wall -fcheck=all	# debug flags
#flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none
#flags = -ftree-vectorize -O3 -ffree-line-length-none
flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none -g -Wall -fcheck=all -fbacktrace	# debug flags

test: $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/test.f90
	gfortran -c $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/test.f90 $(flags)
	gfortran precision.o constants.o functions.o stringlib.o test.o $(flags) -o test.x

isspa2_1D_big: $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/is-spa2.1D.hist.omp.big.f90
	gfortran -c $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/is-spa2.1D.hist.omp.big.f90 $(flags)
	gfortran precision.o constants.o functions.o ideal_CL3.o stringlib.o is-spa2.1D.hist.omp.big.o $(flags) -o is-spa2.1D.hist.omp.big.x

isspa2_1D: $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/is-spa2.1D.hist.omp.f90
	gfortran -c $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/is-spa2.1D.hist.omp.f90 $(flags)
	gfortran precision.o constants.o functions.o ideal_CL3.o stringlib.o is-spa2.1D.hist.omp.o $(flags) -o is-spa2.1D.hist.omp.x

map: $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/map_isspa.f90
	gfortran -c $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/map_isspa.f90 $(flags)
	gfortran precision.o constants.o functions.o stringlib.o map_isspa.o $(flags) -o map_isspa.x

spline_3D: $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.3D.hist.spline.f90
	gfortran -c $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.3D.hist.spline.f90 $(flags)
	gfortran precision.o constants.o functions.o ideal_CL3.o stringlib.o compute_avgForce.3D.hist.spline.o $(flags) -o compute_avgForce.3D.hist.spline.x

spline_1D_parallel: $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.1D.hist.spline.parallel.f90
	gfortran -c $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.1D.hist.spline.parallel.f90 $(flags)
	gfortran precision.o constants.o functions.o ideal_CL3.o stringlib.o compute_avgForce.1D.hist.spline.parallel.o $(flags) -o compute_avgForce.1D.hist.spline.parallel.x

spline_2D_parallel: $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.2D.hist.spline.parallel.f90
	gfortran -c $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.2D.hist.spline.parallel.f90 $(flags)
	gfortran precision.o constants.o functions.o ideal_CL3.o stringlib.o compute_avgForce.2D.hist.spline.parallel.o $(flags) -o compute_avgForce.2D.hist.spline.parallel.x

spline_3D_parallel: $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.3D.hist.spline.parallel.f90
	gfortran -c $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/ideal_CL3.f90 $(func)/stringlib.f90 $(dir)/avg_force_LJ/compute_avgForce.3D.hist.spline.parallel.f90 $(flags)
	gfortran precision.o constants.o functions.o ideal_CL3.o stringlib.o compute_avgForce.3D.hist.spline.parallel.o $(flags) -o compute_avgForce.3D.hist.spline.parallel.x

hist_3D: $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/stringlib.f90 compute_avgForce.3D.hist.f90
	gfortran -c $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/stringlib.f90 compute_avgForce.3D.hist.f90 $(flags)
	gfortran precision.o constants.o functions.o stringlib.o compute_avgForce.3D.hist.o $(flags) -o compute_avgForce.3D.hist.x

ellipse_2D: $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/stringlib.f90 compute_avgForce.2D.ellipse.f90 $(func)/stringlib.f90
	gfortran -c $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/stringlib.f90 compute_avgForce.2D.ellipse.f90 $(flags)
	gfortran precision.o constants.o functions.o stringlib.o compute_avgForce.2D.ellipse.o $(flags) -o compute_avgForce.2D.ellipse.x

sphere_2D: $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/stringlib.f90 compute_avgForce.2D.sphere.f90
	gfortran -c $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/stringlib.f90 compute_avgForce.2D.sphere.f90 $(flags)
	gfortran precision.o constants.o functions.o stringlib.o compute_avgForce.2D.sphere.o $(flags) -o compute_avgForce.2D.sphere.x

hist_2D: $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/stringlib.f90 compute_avgForce.2D.hist.f90
	gfortran -c $(func)/precision.f90 $(func)/constants.f90 $(func)/functions.f90 $(func)/stringlib.f90 compute_avgForce.2D.hist.f90 $(flags)
	gfortran precision.o constants.o functions.o stringlib.o compute_avgForce.2D.hist.o $(flags) -o compute_avgForce.2D.hist.x

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
