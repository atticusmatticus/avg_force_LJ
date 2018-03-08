! USAGE: this_file.x -cfg [CONFIGURATION FILE] -frc [FORCE FILE]
!
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!   Modules   !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! data for the force table.
module frcData
	double precision, allocatable :: fDist(:), fDir(:)
	double precision, allocatable :: R_axis(:), fAvg(:), u_dir(:)
	double precision, allocatable :: x_axis(:), z_axis(:)
	double precision :: fDist_step_size
	integer :: num_R_bins, numFrcLines

endmodule frcData

! data from the config file.
module cfgData
	double precision :: x1, x2, x0, R_step_size, R_min, R_max, xz_range, phi_step_size, cosTheta_step_size, y0, offset

endmodule cfgData

! data for calculating alpha value.
module alphaData
	double precision :: rSolv1(3), rSolv2(3)
	double precision :: rSolv1n, rSolv2n, y02, alp1, alp2

endmodule alphaData


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!  Main Program  !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program compute_avgForce
	implicit none
	character*64 frcFile, cfgFile, outFile
	double precision omp_get_wtime, ti, tf

	ti = omp_get_wtime()

	! make list of average direct force from 'collapsed' file.
	call parse_command_line(frcFile, cfgFile) !, outFile)

	! read config file
	call read_cfg(cfgFile, outFile)

	! make list of average direct force from 'collapsed' file.
	call make_force_table(frcFile)
	
	! compute average force integral.
	call compute_avg_force
	
	! integrate average force to get PMF.
	call integrate_force

	! write PMF output file
	call write_output(outFile)
	
	! Print time taken to finish calculation.
	tf = omp_get_wtime()
	write(*,*) "Total time elapsed: ", tf-ti, "seconds"

endprogram compute_avgForce


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse commandline for relevant files.
subroutine parse_command_line(frcFile,cfgFile) !,outFile)
	implicit none
	character*64 frcFile, cfgFile !, outFile
	character*16 arg
	integer i
	logical frcFileFlag, cfgFileFlag !, outFileFlag

	frcFileFlag = .false.
	cfgFileFlag = .false.
	i=1
	do
		call get_command_argument(i,arg)
		select case (arg)

		case ('-frc')
			i = i+1
			call get_command_argument(i,frcFile)
			frcFileFlag=.true.
			print*, "frc File: ", frcFile
		case ('-cfg')
			i = i+1
			call get_command_argument(i,cfgFile)
			cfgFileFlag=.true.
			print*, "cfg File: ", cfgFile
		case default
			print*, 'Unrecognized command-line option: ', arg
			print*, 'Usage: compute_avgForce.x -frc [frc file]'
			print*, 'Usage: compute_avgForce.x -cfg [cfg file]'
			stop

		end select
		i = i+1
		if (i.ge.command_argument_count()) exit
	enddo

	if (frcFileFlag.eqv..false.) then
		write(*,*) "Must provide a frc file using command line argument -frc [frc file name]"
		stop
	endif

	if (cfgFileFlag.eqv..false.) then
		write(*,*) "Must provide a cfg file using command line argument -cfg [cfg file name]"
		stop
	endif

endsubroutine parse_command_line


! read python cfg file for g(r) parameters
subroutine read_cfg(cfgFile,outFile)
	use cfgData
	implicit none
	character*64 cfgFile, outFile
	character*128 line
	character*32 firstWord, sep
	integer ios
	logical x1Flag, x2Flag, x0Flag, RstepSizeFlag, outFileFlag, RmaxFlag, RminFlag, xzRangeFlag, phiStepSizeFlag
	logical cosThetaStepFlag, y0Flag, offsetFlag

	x1Flag = .false.
	x2Flag = .false.
	x0Flag = .false.
	RstepSizeFlag = .false.
	outFileFlag = .false.
	RmaxFlag = .false.
	RminFlag = .false.
	xzRangeFlag = .false.
	phiStepSizeFlag = .false.
	cosThetaStepFlag = .false.
	y0Flag = .false.
	offsetFlag = .false.

	ios = 0

	open(20,file=cfgFile)
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		call split(line,'=',firstWord, sep)
		if (line .ne. "") then
			if (firstWord .eq. "x1") then
				read(line,*) x1
				write(*,*) "x1 Parameter:	", x1
				x1Flag = .true.
			else if (firstWord .eq. "x2") then
				read(line,*) x2
				write(*,*) "x2 Parameter:	", x2
				x2Flag = .true.
			else if (firstWord .eq. "x0") then
				read(line,*) x0
				write(*,*) "x0 Parameter:	", x0
				x0Flag = .true.
			else if (firstWord .eq. "R_step_size") then
				read(line,*) R_step_size
				write(*,*) "PMF Step Size:	", R_step_size
				RstepSizeFlag = .true.
			else if (firstWord .eq. "out_file") then
				read(line,*) outFile
				write(*,*) "Output File:	", outFile
				outFileFlag = .true.
			else if (firstWord .eq. "R_max") then
				read(line,*) R_max
				write(*,*) "R Maximum Value:	", R_max
				RmaxFlag = .true.
			else if (firstWord .eq. "R_min") then
				read(line,*) R_min
				write(*,*) "R Minimum Value:	", R_min
				RminFlag = .true.
			else if (firstWord .eq. "xz_range") then
				read(line,*) xz_range
				write(*,*) "XZ - Range:	", xz_range
				xzRangeFlag = .true.
			else if (firstWord .eq. "phi_step_size") then
				read(line,*) phi_step_size
				write(*,*) "Phi Step Size:	", phi_step_size
				phiStepSizeFlag= .true.
			else if (firstWord .eq. "cosTheta_step_size") then
				read(line,*) y0
				write(*,*) "Cosine Theta Step Size:	", cosTheta_step_size
				cosThetaStepFlag = .true.
			else if (firstWord .eq. "offset") then
				read(line,*) offset
				write(*,*) "Offset Along Minor Axis:	", offset
				offsetFlag = .true.
			else if (firstWord .eq. "y0") then
				read(line,*) y0
				write(*,*) "Spheroid Minor Axis:	", y0
				y0Flag = .true.
			endif
		endif
	enddo
	close(20)

	if (x1Flag.eqv..false.) then
		write(*,*) "Config file must have a 'x1' value"
		stop
	endif
	if (x2Flag.eqv..false.) then
		write(*,*) "Config file must have a 'x2' value"
		stop
	endif
	if (x0Flag.eqv..false.) then
		write(*,*) "Config file must have a 'x0' value"
		stop
	endif
	if (RstepSizeFlag.eqv..false.) then
		write(*,*) "Config file must have a 'R_step_size' value"
		stop
	endif
	if (outFileFlag.eqv..false.) then
		write(*,*) "Config file must have a 'out_file' value"
		stop
	endif
	if (RmaxFlag.eqv..false.) then
		write(*,*) "Config file must have a 'R_max' value"
		stop
	endif
	if (RminFlag.eqv..false.) then
		write(*,*) "Config file must have a 'R_min' value"
		stop
	endif
	if (xzRangeFlag.eqv..false.) then
		write(*,*) "Config file must have a 'xz_range' value"
		stop
	endif
	if (phiStepSizeFlag.eqv..false.) then
		write(*,*) "Config file must have a 'phi_step_size' value"
		stop
	endif
	if (cosThetaStepFlag.eqv..false.) then
		write(*,*) "Config file must have a 'cosTheta_step_size' value"
		stop
	endif
	if (offsetFlag.eqv..false.) then
		write(*,*) "Config file must have a 'offset' value"
		stop
	endif
	if (y0Flag.eqv..false.) then
		write(*,*) "Config file must have a 'y0' value"
		stop
	endif

endsubroutine read_cfg


! read force file and make a lookup table.
subroutine make_force_table(frcFile)
	use frcData
	implicit none
	character*64 frcFile
	character*32 junk
	character*256 line
	integer ios, i


	! read number of lines in frcFile and allocate that many points in force list.
	ios = 0
	numFrcLines = 0
	open(20,file=frcFile)
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		if (line(1:1) .ne. "#") then
			numFrcLines = numFrcLines + 1
		endif
	enddo
	close(20)

	allocate( fDist(numFrcLines), fDir(numFrcLines) )


	! populate force list
	ios = 0
	i = 1
	open(20,file=frcFile)
	! read file ignoring comment lines at the beginning
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		if (line(1:1) .ne. "#") then
			read(line,*) fDist(i), junk, junk, fDir(i)
!			write(6,*) fDist(i), fDir(i)
			i = i + 1
		endif
	enddo
	close(20)
	
	fDist_step_size = fDist(2) - fDist(1)

endsubroutine make_force_table


! do the average force integral
subroutine compute_avg_force
	use frcData
	use cfgData
	use alphaData
	implicit none
	integer num_xz_bins, r, i, j, thLF, phLF
	double precision :: gx, lin_out, pi

	pi = 3.1415926535
	y02 = y0*y0

	num_R_bins = int( (R_max - R_min)/R_step_size )
	write(*,*) "Number of R Bins: ", num_R_bins
	num_xz_bins = int( (xz_range + xz_range)/R_step_size )
	write(*,*) "Number of XZ Bins: ", num_xz_bins

	allocate( R_axis(num_R_bins), fAvg(num_R_bins), x_axis(num_xz_bins), z_axis(num_xz_bins) )
	fAvg = 0

	do r = 1, num_R_bins
		R_axis(r) = r * R_step_size + R_min
	enddo
	do i = 1, num_xz_bins
		x_axis(i) = i * R_step_size - xz_range
	enddo
	do j = 1, int(num_xz_bins/2.)
		z_axis(j) = j * R_step_size
	enddo

	! Calculate the average force integral for top half of cylinder
	do r = 1, num_R_bins ! loop lj--lj distances
		do i = 1, num_xz_bins ! loop solvent locations
			do j = 1, num_xz_bins ! loop solvent locations
				rSolv1(1) = x_axis(i)-R_axis(r)/2.0
				rSolv1(2) = 0.0
				rSolv1(3) = z_axis(j)
				rSolv1n = sqrt( rSolv1(1)**2 + rSolv1(2)**2 + rSolv1(3)**2 )
				!rSolv1n = norm2(rSolv1)

				rSolv2(1) = x_axis(i)+R_axis(r)/2.0
				rSolv2(2) = 0.0
				rSolv2(3) = z_axis(j)
				rSolv2n = sqrt( rSolv2(1)**2 + rSolv2(2)**2 + rSolv2(3)**2 )
				!rSolv2n = norm2(rSolv2)

				! note: loop through orientations of solvent at x(i) and z(j)
				do thLF = 0, 49 ! if the limit is exclusive. 49 if inclusive.
					do phLF = 0, 100 ! same as th1

						call alpha(phLF, thLF)

						call g12(alp1, alp2, gx)
						! if gx = 0 then don't waste time with the rest of the calculation
						!write(*,*) 'alpha1: ', alp1, ' alpha2: ', alp2, ' gx: ', gx
						if (gx .gt. 1d-6) then
							call lin_interpolate(alp1, lin_out)

							! NOTES : 	'gx' is Kirkwood Super Position Approximation of g(r)
							! 			'lin_out' is |fs(x,z)|, force from solvent at (x,z)
							!			Now we need cos(theta) and z and we should have the integral.

							! FIXME:	what happens to this integral when we change to alpha1 and alpha2 ???
							fAvg(r) = fAvg(r) + ( gx * (-1)*lin_out * z_axis(j) * ( (x_axis(i)-(R_axis(r)/2.0)) / rSolv1n ) )
						endif
					enddo
				enddo


			enddo
		enddo

		! NOTE : After the fact multiply all elements by 2*pi*density*r**2
		! 		Number density of chloroform per Angstrom**2 == 0.00750924
		fAvg(r) = fAvg(r)*2*pi*0.00750924*R_step_size**2
	enddo

endsubroutine compute_avg_force


! compute the superposition g1(r1)*g2(r2) depending on what range the solvent is in.
subroutine g12(rSolv1n, rSolv2n, gx)
	use cfgData
	implicit none
	double precision :: rSolv1n, rSolv2n, parabola
	double precision :: gx

	if ( (rSolv1n .le. x1) .or. (rSolv2n .le. x1) ) then
		gx = 0
	elseif ( (rSolv1n .gt. x1) .and. (rSolv1n .lt. x2) .and. (rSolv2n .ge. x2) ) then
		gx = parabola(rSolv1n)
	elseif ( (rSolv1n .gt. x1) .and. (rSolv1n .lt. x2) .and. (rSolv2n .gt. x1) .and. (rSolv2n .lt. x2) ) then
		gx = parabola(rSolv1n) * parabola(rSolv2n)
	elseif ( (rSolv1n .ge. x2) .and. (rSolv2n .gt. x1) .and. (rSolv2n .lt. x2) ) then
		gx = parabola(rSolv2n)
	elseif ( (rSolv1n .ge. x2) .and. (rSolv2n .ge. x2) ) then
		gx = 1
	endif
	
endsubroutine g12


! parabola function for g(r) fit.
function parabola(x) result(p)
	use cfgData
	implicit none
	double precision, intent(in)	:: x	! inputs not to be changed
	double precision				:: p	! ouput

	p = ( (x1-x0)**2 - (x-x0)**2 )/ ( (x1-x0)**2 - (x2-x0)**2 )

endfunction parabola


! rotate solvent vector 'p' and calculate alpha for lj particles 1 and 2.
subroutine alpha(phLF, thLF)
	use cfgData
	use alphaData
	implicit none
	integer				:: phLF, thLF
	double precision	:: phiLF, cosThetaLF
	double precision	:: p(3)
	double precision	:: cosTh1, cosTh2

	! Lab Frame angles
	phiLF = phLF*phi_step_size
	cosThetaLF = thLF*cosTheta_step_size

	! make rotated solvent vector at origin
	p(1) = sin(phiLF)*sqrt((1.0-cosThetaLF**2))
	p(2) = cos(phiLF)*sqrt((1.0-cosThetaLF**2))
	p(3) = cosThetaLF

	! calculate cos(theta1) and cos(theta2) relative to lj-spheres 1 and 2.
	cosTh1 = dot_product(rSolv1, p) / rSolv1n
	cosTh2 = dot_product(rSolv2, p) / rSolv2n

	alp1 = sqrt( rSolv1n**2 * (1.0-cosTh1**2) + (rSolv1n*cosTh1 - offset)**2 / y02 )
	alp2 = sqrt( rSolv2n**2 * (1.0-cosTh2**2) + (rSolv2n*cosTh2 - offset)**2 / y02 )

endsubroutine alpha


! linearly interpolate between force at 'a' and 'b' to give 'output' force.
subroutine lin_interpolate(rSolv, lin_out)
	use frcData
	implicit none
	integer i1, i2	! the left and right flanking indices
	double precision float_index, rSolv, lin_out

	float_index = rSolv / fDist_step_size
	i1 = floor(float_index)
	i2 = i1 + 1

	if (i2 .ge. numFrcLines) then
		lin_out = 0
	else
		lin_out = (( fDir(i2) - fDir(i1) ) / ( fDist(i2) - fDist(i1) )) * (rSolv - fDist(i1)) + fDir(i1)
	endif

	!write(*,*) 'rSolv: ', rSolv, '   lin_out: ', lin_out
endsubroutine lin_interpolate


! integrate the average force from 'compute_avg_force' to get the PMF.
subroutine integrate_force
	use frcData
	use cfgData
	implicit none
	integer d

	allocate( u_dir(num_R_bins) )
	u_dir = 0

	do d = 0, num_R_bins
		if (d .eq. 0) then
			u_dir(num_R_bins) = fAvg(num_R_bins) * R_step_size
		else
			u_dir(num_R_bins-d) = u_dir(num_R_bins-(d-1)) + fAvg(num_R_bins-d) * R_step_size
		endif
	enddo

endsubroutine integrate_force


! write output file
subroutine write_output(outFile)
	use frcData
	implicit none
	character*64 outFile
	integer i

	open(25,file=outFile)
	write(6,*) "Writing output file:	", outFile
	write(25,*) "# 1.	R Distance"
	write(25,*) "# 2.	Avg Force"
	write(25,*) "# 3.	PMF"
	do i = 1, num_R_bins
		write(25,899) R_axis(i), fAvg(i), u_dir(i)
	enddo
	close(25)

899		format (f20.5,2(1x,f20.5))

endsubroutine write_output
