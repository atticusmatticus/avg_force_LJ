
dir = ~/programs
func = $(dir)/functions

flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none

test: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 test.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 test.f90 $(flags)
	gfortran precision.o functions.o constants.o stringlib.o test.o $(flags) -o test.x

testlin: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 test.lin.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 test.lin.f90 $(flags)
	gfortran precision.o functions.o constants.o stringlib.o test.lin.o $(flags) -o test.lin.x

testlog: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 test.log.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 test.log.f90 $(flags)
	gfortran precision.o functions.o constants.o stringlib.o test.log.o $(flags) -o test.log.x

single2D: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 singleSolute.2D.hist.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 singleSolute.2D.hist.f90 $(flags)
	gfortran precision.o functions.o constants.o stringlib.o singleSolute.2D.hist.o $(flags) -o singleSolute.2D.hist.x

single3D: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 singleSolute.3D.hist.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 singleSolute.3D.hist.f90 $(flags)
	gfortran precision.o functions.o constants.o stringlib.o singleSolute.3D.hist.o $(flags) -o singleSolute.3D.hist.x

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

ideal3D: $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 compute_avgForce.ideal.f90
	gfortran -c $(func)/precision.f90 $(func)/functions.f90 $(func)/constants.f90 $(func)/stringlib.f90 compute_avgForce.ideal.f90 $(flags)
	gfortran precision.o functions.o constants.o stringlib.o compute_avgForce.ideal.o $(flags) -o compute_avgForce.ideal.x

clean : 
	rm -f *.o *.mod *.x
