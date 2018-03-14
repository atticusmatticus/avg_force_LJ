! USAGE: ./this_file.x -cfg [CONFIGURATION FILE] -frc [FORCE FILE]
!
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!   Modules   !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! data for the force table.
module frcData
	double precision, allocatable	:: fDist(:), fDir(:)
	double precision, allocatable	:: R_axis(:), fAvg(:), u_dir(:)
	double precision, allocatable	:: x_axis(:), z_axis(:)
	double precision				:: fDist_step_size
	integer							:: num_R_bins, numFrcLines

endmodule frcData

! data from the config file.
module cfgData
	double precision :: x1, x2, x0, R_step_size, xz_step_size, R_min, R_max, xz_range, offset

endmodule cfgData

! testing arrays for force and g(r)
module testData
	double precision, allocatable :: linAvg(:,:), grSPA(:,:)

endmodule testData


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!  Main Program  !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program compute_avgForce
	implicit none
	character(len=64) 		:: frcFile, cfgFile, outFile
	double precision 		:: omp_get_wtime, ti, tf

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
subroutine parse_command_line(frcFile, cfgFile) !,outFile)
	implicit none
	character(len=64) 	:: frcFile, cfgFile !, outFile
	character(len=16) 	:: arg
	integer 			:: i
	logical 			:: frcFileFlag, cfgFileFlag !, outFileFlag

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
subroutine read_cfg(cfgFile, outFile)
	use cfgData
	implicit none
	character(len=64) 	:: cfgFile, outFile
	character(len=128) 	:: line
	character(len=32) 	:: firstWord, sep
	integer 			:: ios
	logical 			:: x1Flag, x2Flag, x0Flag, RstepSizeFlag, xzStepSizeFlag, outFileFlag, RmaxFlag, RminFlag
	logical 			:: xzRangeFlag, offsetFlag

	x1Flag = .false.
	x2Flag = .false.
	x0Flag = .false.
	RstepSizeFlag = .false.
	xzStepSizeFlag = .false.
	outFileFlag = .false.
	RmaxFlag = .false.
	RminFlag = .false.
	xzRangeFlag = .false.
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
			else if (firstWord .eq. "xz_step_size") then
				read(line,*) xz_step_size
				write(*,*) "Solvent Grid Step Size:	", xz_step_size
				xzStepSizeFlag = .true.
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
			else if (firstWord .eq. "offset") then
				read(line,*) offset
				write(*,*) "Offset Along Minor Axis:	", offset
				offsetFlag = .true.
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
	if (xzStepSizeFlag.eqv..false.) then
		write(*,*) "Config file must have a 'xz_step_size' value"
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
	if (offsetFlag.eqv..false.) then
		write(*,*) "Config file must have a 'offset' value"
		stop
	endif

endsubroutine read_cfg


! read force file and make a lookup table.
subroutine make_force_table(frcFile)
	use frcData
	implicit none
	character(len=64) 	:: frcFile
	character(len=32) 	:: junk
	character(len=256) 	:: line
	integer 			:: ios, i


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
	use testData
	implicit none
	integer 			:: num_x_bins, num_z_bins, r, i, j
	double precision	:: gx, lin_out, pi, rSolv1n, rSolv2n, rSolv1(3), rSolv2(3)

	pi = 3.1415926535

	num_R_bins = int( (R_max - R_min)/R_step_size )
	write(*,*) "Number of R Bins: ", num_R_bins
	num_x_bins = int( (xz_range + xz_range)/xz_step_size )
	write(*,*) "Number of X Bins: ", num_x_bins
	num_z_bins = int( (xz_range)/xz_step_size )
	write(*,*) "Number of Z Bins: ", num_z_bins

	! allocate array sizes for axes and average force
	allocate( R_axis(num_R_bins), fAvg(num_R_bins), x_axis(num_x_bins), z_axis(num_z_bins) )
	R_axis = 0.d0
	fAvg = 0.d0
	x_axis = 0.d0
	z_axis = 0.d0

	! allocate arrays for testing arrays
	allocate( linAvg(num_x_bins, num_z_bins), grSPA(num_x_bins, num_z_bins) )
	linAvg = 0.d0
	grSPA = 0.d0

	! Distance Axes
	do r = 1, num_R_bins
		R_axis(r) = (r-1) * R_step_size + R_min
	enddo
	do i = 1, num_x_bins
		x_axis(i) = (i-1) * xz_step_size - xz_range ! fixme xz_step_size
	enddo
	do j = 1, num_z_bins
		z_axis(j) = (j-1) * xz_step_size ! fixme xz_step_size
	enddo

	! Calculate the average force integral for top half of cylinder
	do r = 1, num_R_bins ! loop lj--lj distances
		linAvg = 0.d0
		grSPA = 0.d0
		do i = 1, num_x_bins ! full length of cylinder
			do j = 1, num_z_bins ! top half of bisecting plane of cylinder
				rSolv1(1) = x_axis(i)-R_axis(r)/2.0
				rSolv1(2) = 0.0
				rSolv1(3) = z_axis(j)
				rSolv1n = norm2(rSolv1)

				rSolv2(1) = x_axis(i)+R_axis(r)/2.0
				rSolv2(2) = 0.0
				rSolv2(3) = z_axis(j)
				rSolv2n = norm2(rSolv2)

				call g12(rSolv1n, rSolv2n, gx)
				!write(*,*) 'rSolv1n: ', rSolv1n, ' rSolve2n: ', rSolv2n, ' gx: ', gx

				! if gx == 0 then don't waste time with the rest of the calculation
				if (gx .gt. 1d-6) then
					call lin_interpolate(rSolv1n, lin_out)

					! NOTES : 	'gx' is Kirkwood Super Position Approximation of g(r)
					! 			'lin_out' is ||fs(x,z)||, force from solvent at (x,z)
					!			Now we need cos(theta) and z and we should have the integral.

					fAvg(r) = fAvg(r) + ( gx * (-1)*lin_out * z_axis(j) * ( rSolv1(1) / rSolv1n ) )
					linAvg(i,j) = linAvg(i,j) + lin_out ! fixme print one of these to file for every r
					grSPA(i,j) = grSPA(i,j) + gx ! fixme same with gx
				endif
			enddo !z
		enddo !x
		call write_test_out(r, num_x_bins, num_z_bins) ! write grSPA and lin_out arrays

		! NOTE : After the fact multiply all elements by 2*pi*density*r**2
		! 		Number density of chloroform per Angstrom**3 == 0.00750924
		fAvg(r) = fAvg(r)*2*pi*0.00750924*xz_step_size**2
	enddo !r

endsubroutine compute_avg_force


! compute the superposition g1(r1)*g2(r2) depending on what range the solvent is in.
subroutine g12(dist1, dist2, gx)
	use cfgData
	implicit none
	double precision :: parabola, dist1, dist2
	double precision :: gx

	if ( (dist1 .le. x1) .or. (dist2 .le. x1) ) then
		gx = 0
	elseif ( (dist1 .gt. x1) .and. (dist1 .lt. x2) .and. (dist2 .ge. x2) ) then
		gx = parabola(dist1)
	elseif ( (dist1 .gt. x1) .and. (dist1 .lt. x2) .and. (dist2 .gt. x1) .and. (dist2 .lt. x2) ) then
		gx = parabola(dist1) * parabola(dist2)
	elseif ( (dist1 .ge. x2) .and. (dist2 .gt. x1) .and. (dist2 .lt. x2) ) then
		gx = parabola(dist2)
	elseif ( (dist1 .ge. x2) .and. (dist2 .ge. x2) ) then
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


! linearly interpolate between force at 'a' and 'b' to give 'output' force.
subroutine lin_interpolate(alp, lin_out)
	use frcData
	implicit none
	integer 			:: i1, i2	! the left and right flanking indices
	double precision 	:: float_index, alp, lin_out

	float_index = alp / fDist_step_size + 0.5 ! fixme forces come from half bin position. print lin_out(alp) compare to MD sim.
	i1 = floor(float_index)
	i2 = i1 + 1

	if (i2 .ge. numFrcLines) then
		lin_out = 0
	else
		lin_out = (( fDir(i2) - fDir(i1) ) / ( fDist(i2) - fDist(i1) )) * (alp - fDist(i1)) + fDir(i1)
	endif

	!write(*,*) 'alp: ', alp, '   lin_out: ', lin_out
endsubroutine lin_interpolate


! integrate the average force from 'compute_avg_force' to get the PMF.
subroutine integrate_force
	use frcData
	use cfgData
	implicit none
	integer 		:: d

	allocate( u_dir(num_R_bins) )
	u_dir = 0.d0

	do d = 1, num_R_bins
		if (d .eq. 1) then
			u_dir(num_R_bins) = fAvg(num_R_bins) * R_step_size
		else
			u_dir(num_R_bins-(d-1)) = u_dir(num_R_bins-(d-2)) + fAvg(num_R_bins-(d-1)) * R_step_size
		endif
	enddo

endsubroutine integrate_force


! write output file
subroutine write_output(outFile)
	use frcData
	implicit none
	character(len=64) 		:: outFile
	integer 				:: r

	open(35,file=outFile)
	write(6,*) "Writing output file:	", outFile
	write(35,*) "# 1.	R Distance"
	write(35,*) "# 2.	Avg Force"
	write(35,*) "# 3.	PMF"
	do r = 1, num_R_bins
		write(35,899) R_axis(r), fAvg(r), u_dir(r)
	enddo
	close(35)

899		format (3(1x,f16.12))
!899		format (3(1x,e20.10)) ! scientific format

endsubroutine write_output


! write force out and g(r) out to compare against explicit
subroutine write_test_out(r, num_x_bins, num_z_bins)
	use frcData
	use testData
	implicit none
	integer				:: r, i, j, num_x_bins, num_z_bins
	character(len=32)	:: temp, filename
	character(len=8)	:: frmt

	frmt = '(I3.3)' ! an integer of width 3 with zeroes on the left
	write(temp,frmt) r ! converting integer to string using 'internal file'
	filename='sphere_output.'//trim(temp)//'.dat'


	open(35,file=filename)
	write(6,*) "Writing test file:	", filename
	write(35,*) "# 1.	X Distance"
	write(35,*) "# 2.	Z Distance"
	write(35,*) "# 3.	Avg Force List"
	write(35,*) "# 4.	g(r) List"
	do i = 1, num_x_bins
		do j = 1, num_z_bins
			write(35,898) x_axis(i), z_axis(j), linAvg(i,j), grSPA(i,j)
		enddo
	enddo
	close(35)

898		format (4(1x,f16.12))

endsubroutine write_test_out
