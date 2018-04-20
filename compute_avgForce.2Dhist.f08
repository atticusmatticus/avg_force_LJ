! USAGE: ./this_file.x -cfg [CONFIGURATION FILE] -hist [2D HIST FILE]
!
! using 2D histrograms of spherically binned ie. g(r,cos(th)) for both PDF and Force
!
!					    ^
!					   z|
!						|
!						|
!						| /y
!						|/
!		<_______O_______|_______O_______>
!		-x		2				1	+x
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!   Modules   !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! precision module
module prec
	integer, parameter	:: dp = kind(1.0d0)
endmodule prec

! data for the density and force tables.
module histData
	use prec
	real(kind=dp), allocatable	:: histDist(:), histAng(:), histFrc(:,:), histGr(:,:), R_axis(:), fAvg(:), u_dir(:), &
		& x_axis(:), z_axis(:)
	real(kind=dp)				:: histDist_step_size
	integer						:: num_R_bins, numFrcLines

endmodule histData

! data from the config file.
module cfgData
	use prec
	real(kind=dp)	:: R_step_size, xz_step_size, R_min, R_max, xz_range, phi_step_size, cosTheta_step_size
	integer			:: nThetaBins, phi_bins

endmodule cfgData

! data for calculating cosTh value.
module thetaData
	use prec
	real(kind=dp), allocatable	:: cosThetaLF(:), sinThetaLF(:), sinPhiLF(:), cosPhiLF(:)
	real(kind=dp)				:: rSolv1(3), rSolv2(3), rSolv1n, rSolv2n

endmodule thetaData

! testing arrays for force and g(r)
module ctrlData
	use prec
	real(kind=dp), allocatable :: linAvg(:,:), grSPA(:,:), fxArray(:)

endmodule ctrlData


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!  Main Program  !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program compute_avgForce
	use prec
	implicit none
	character(len=64) 	:: histFile, cfgFile, outFile
	real(kind=dp) 		:: omp_get_wtime, ti, tf

	ti = omp_get_wtime()

	! make list of average direct force from 'collapsed' file.
	call parse_command_line(histFile, cfgFile) !, outFile)

	! read config file
	call read_cfg(cfgFile, outFile)

	! make list of average direct force from 'collapsed' file.
	call make_hist_table(histFile)
	
	! compute average force integral.
	!call compute_avg_force
	
	! integrate average force to get PMF.
	!call integrate_force

	! write PMF output file
	!call write_output(outFile)
	
	! Print time taken to finish calculation.
	tf = omp_get_wtime()
	write(*,*) "Total time elapsed: ", tf-ti, "seconds"

endprogram compute_avgForce


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse commandline for relevant files.
subroutine parse_command_line(histFile, cfgFile) !,outFile)
	implicit none
	character(len=64) 	:: histFile, cfgFile !, outFile
	character(len=16) 	:: arg
	integer 			:: i
	logical 			:: histFileFlag, cfgFileFlag !, outFileFlag

	histFileFlag = .false.
	cfgFileFlag = .false.
	i=1
	do
		call get_command_argument(i,arg)
		select case (arg)

		case ('-hist')
			i = i+1
			call get_command_argument(i,histFile)
			histFileFlag=.true.
			print*, "hist File: ", histFile
		case ('-cfg')
			i = i+1
			call get_command_argument(i,cfgFile)
			cfgFileFlag=.true.
			print*, "cfg File: ", cfgFile
		case default
			print*, 'Unrecognized command-line option: ', arg
			print*, 'Usage: compute_avgForce.x -hist [hist file]'
			print*, 'Usage: compute_avgForce.x -cfg [cfg file]'
			stop

		end select
		i = i+1
		if (i.ge.command_argument_count()) exit
	enddo

	if (histFileFlag.eqv..false.) then
		write(*,*) "Must provide a hist file using command line argument -hist [hist file name]"
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
	logical 			:: RstepSizeFlag, xzStepSizeFlag, outFileFlag, RmaxFlag, RminFlag, xzRangeFlag, phiBinsFlag

	RstepSizeFlag = .false.
	xzStepSizeFlag = .false.
	outFileFlag = .false.
	RmaxFlag = .false.
	RminFlag = .false.
	xzRangeFlag = .false.
	phiBinsFlag = .false.

	ios = 0

	open(20,file=cfgFile)
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		call split(line,'=',firstWord, sep)
		if (line .ne. "") then
			if (firstWord .eq. "R_step_size") then
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
			else if (firstWord .eq. "phi_bins") then
				read(line,*) phi_bins
				write(*,*) "Phi Bins:	", phi_bins
				phiBinsFlag= .true.
			endif
		endif
	enddo
	close(20)

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
	if (phiBinsFlag.eqv..false.) then
		write(*,*) "Config file must have a 'phi_bins' value"
		stop
	endif

endsubroutine read_cfg


! read force file and make a lookup table.
subroutine make_hist_table(histFile)
	use histData
	implicit none
	character(len=64)			:: histFile
	character(len=32)			:: junk
	character(len=256)			:: line
	integer						:: ios, i
	real(kind=dp), allocatable	:: histTmp(:,:)

	! read number of lines in histFile and allocate that many points in temporary histogram list.
	ios = 0
	nHistLines = 0
	open(20,file=histFile)
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		if (line(1:1) .ne. "#") then
			nHistLines = nHistLines + 1
		endif
	enddo
	close(20)

!	nDistBins = int(dble(numHistLines) / dble(nThetaBins))
!	print*, "Number of distance bins in hist table: should be 250!!: ", nDistBins
	! XXX: I wonder if i could read in each r value from the hist table and see if it's within .000001 of the previous value read
	! in. If so then ignore it, if not then add it to array and increment nDistBins variable. This would be a good way for getting
	! the nDistBins and nThetaBins from the table itself rather than the config file or hard coding.

	allocate( histTmp(4,nHistLines) )

	! populate hist arrays
	ios = 0
	i = 1
	open(20,file=histFile)
	! read file ignoring comment lines at the beginning
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		if (line(1:1) .ne. "#") then
			!			 dist		   cos(Th)		 g(r)+		   g(r)- force(r)+
			read(line,*) histTmp(1,i), histTmp(2,i), histTmp(3,i), junk, histTmp(4,i) ! FIXME: best way to populate these arrays?
			i = i + 1
		endif
	enddo
	close(20)

	! XXX: Unique value determination -- TESTING
	do i = 1, nHistLines
		if (i .eq. 1) then
			nDistBins = 1
			ios = 0
		else ! i = 2, nHistLines
			if ( (histTmp(1,i) .lt. (histTmp(1,i-1)-1d-6)) .or. (histTmp(1,i) .gt. (histTmp(1,i-1)+1d-6)) ) then
				nDistBins = nDistBins + 1
				print *, "nDistBins: ", nDistBins
			endif
			if ( (histTmp(2,i) .gt. (histTmp(2,1)-1d-6)) .and. (histTmp(2,i) .lt. (histTmp(2,1)+1d-6)) .and. (ios .eq. 0) ) then
				! xxx: this statement will trigger when i = nThetaBins + 1 because it finds the first repeated element
				nThetaBins = i - 1
				print *, "nThetaBins: ", nThetaBins
				ios = 1
			endif
		endif
	enddo
	
	allocate( histDist(nDistBins), histAng(nThetaBins), histFrc(nDistBins,nThetaBins), histGr(nDistBins,nThetaBins) )
	
	! populate arrays that will be used in the rest of the calculation from temp array
	do i = 1, nDistBins		! the values printed out from python script are at half-bin distances
		histDist(i) = histTmp(1,nThetaBins*i-1)
	enddo
	do i = 1, nThetaBins
		histAng(i) = histTmp(2,i)
	enddo
	do i = 1, nDistBins
		do j = 1, nThetaBins
			histGr(i,j) = histTmp(3,i+j)
			histFrc(i,j) = histTmp(4,i+j)
		enddo
	enddo

	histDist_step_size = histDist(2) - histDist(1)

endsubroutine make_hist_table


! do the average force integral
subroutine compute_avg_force
	use histData
	use cfgData
	use thetaData
	use ctrlData
	implicit none
	integer 		:: num_x_bins, num_z_bins, r, i, j, ithLF, iphiLF
	real(kind=dp)	:: gx, fx, pi, phiLF, phi_max, phi_min

	pi = 3.1415926535_dp

	num_R_bins = int( (R_max - R_min)/R_step_size )
	write(*,*) "Number of R Bins: ", num_R_bins
	num_x_bins = int( (2_dp * xz_range)/xz_step_size )
	write(*,*) "Number of X Bins: ", num_x_bins
	num_z_bins = int( (xz_range)/xz_step_size )
	write(*,*) "Number of Z Bins: ", num_z_bins

	! allocate array sizes for axes and average force
	allocate( R_axis(num_R_bins), fAvg(num_R_bins), x_axis(num_x_bins), z_axis(num_z_bins) )
	R_axis = 0_dp
	fAvg = 0_dp
	x_axis = 0_dp
	z_axis = 0_dp

	! allocate arrays for testing arrays
	allocate( linAvg(num_x_bins, num_z_bins), grSPA(num_x_bins, num_z_bins), fxArray(num_x_bins) )
	fxArray = 0_dp
	linAvg = 0_dp
	grSPA = 0_dp

	! Distance Axes
	do r = 1, num_R_bins
		R_axis(r) = (r-1) * R_step_size + R_min
	enddo
	do i = 1, num_x_bins
		x_axis(i) = (i-1) * xz_step_size - xz_range
	enddo
	do j = 1, num_z_bins
		z_axis(j) = (j-1) * xz_step_size
	enddo

	! Angles
	allocate( cosThetaLF(nThetaBins), sinThetaLF(nThetaBins), sinPhiLF(phi_bins), cosPhiLF(phi_bins) )
	cosTheta_step_size = (cosTh_max - cosTh_min) / real(nThetaBins, dp)
	do ithLF = 1, nThetaBins
		cosThetaLF(ithLF) = (ithLF-0.5_dp)*cosTheta_step_size - 1_dp
		sinThetaLF(ithLF) = dsqrt(abs(1_dp-cosThetaLF(ithLF)**2))
	enddo
	write(*,*) "Cos(Theta) Step Size: ", cosTheta_step_size
	phi_max = 2_dp*pi
	phi_min = 0_dp
	phi_step_size = (phi_max - phi_min) / real(phi_bins, dp)
	do iphiLF = 1, phi_bins
		phiLF = (iphiLF+0.5_dp)*phi_step_size
		sinPhiLF(iphiLF) = dsin(phiLF)
		cosPhiLF(iphiLF) = dcos(phiLF)
	enddo
	write(*,*) "Phi Step Size: ", phi_step_size

	! Calculate the average force integral for top half of cylinder
	do r = 1, num_R_bins ! loop lj--lj distances
		linAvg = 0_dp
		grSPA = 0_dp
		do i = 1, num_x_bins ! full length of cylinder
			do j = 1, num_z_bins ! top half of bisecting plane of cylinder
				rSolv1(1) = x_axis(i)-R_axis(r)/2_dp
				rSolv1(2) = 0_dp
				rSolv1(3) = z_axis(j)
				rSolv1n = norm2(rSolv1)

				rSolv2(1) = x_axis(i)+R_axis(r)/2_dp
				rSolv2(2) = 0_dp
				rSolv2(3) = z_axis(j)
				rSolv2n = norm2(rSolv2)

				! loop through orientations of solvent at x(i) and z(j)
				do ithLF = 1, nThetaBins
					do iphiLF = 1, phi_bins

						if ((rSolv1n .lt. 1d-6) .or. (rSolv2n .lt. 1d-6)) then
							gx = 0_dp ! avoid NaNs in calc_cosTh
						else
							call calc_cosTh(iphiLF, ithLF)

							! FIXME: bilinearly interpolate g(r, cos)
							!call g12(gx)
							call bilin_interpolate(histGr, gx)
						endif

						! if gx == 0 then don't waste time with the rest of the calculation
						if (gx .gt. 1d-6) then
							! FIXME: bilinearly interpolate frc(r, cos)
							!call lin_interpolate(alp1, lin_out)
							call bilin_interpolate(histFrc, fx)

							!			'gx' is Kirkwood Super Position Approximation of g(r)
							! 			'fx' is ||fs(x,z)||, force from solvent at (x,z). Now we
							!			need cos(theta) and z and we should have the integral.

							fAvg(r) = fAvg(r) + ( gx * (-1)*fx * z_axis(j) * ( rSolv1(1) / rSolv1n ) )
							linAvg(i,j) = linAvg(i,j) + fx
							fxArray(i) = fxArray(i) + fx*xz_step_size
							grSPA(i,j) = grSPA(i,j) + gx
						endif
					enddo !phi
				enddo !theta
			enddo !z
		enddo !x
		call write_test_out(r, num_x_bins, num_z_bins) ! write grSPA and fx arrays

		! NOTE : After the fact multiply all elements by 2*pi*density*/4/pi (4pi steradians from orientations)
		! 		Number density of chloroform per Angstrom**3 == 0.00750924
		fAvg(r) = fAvg(r)/2_dp*0.00750924_dp*xz_step_size*xz_step_size*phi_step_size*cosTheta_step_size
	enddo !r

endsubroutine compute_avg_force


! rotate solvent vector 'p' and calculate cosTh for lj particles 1 and 2.
subroutine calc_cosTh(iphiLF, ithLF)
	use cfgData
	use thetaData
	implicit none
	integer			:: iphiLF, ithLF
	real(kind=dp)	:: p(3), cosTh1, cosTh2

	! make rotated solvent vector at origin
	p(1) = sinPhiLF(iphiLF)*sinThetaLF(ithLF)
	p(2) = cosPhiLF(iphiLF)*sinThetaLF(ithLF)
	p(3) = cosThetaLF(ithLF)

	! calculate cos(theta1) and cos(theta2) relative to lj-spheres 1 and 2.
	cosTh1 = dot_product(rSolv1, p) / rSolv1n
	cosTh2 = dot_product(rSolv2, p) / rSolv2n ! fixme is this needed?

endsubroutine calc_cosTh


! bilinearly interpolate
subroutine bilin_interpolate(hist_array, bi_out)
	use thetaData
	implicit none
	integer			:: ir1, ir2, ic1, ic2
	real(kind=dp)	:: float_index, ra, rb, rd, cd, fR1, fR2, fP

	float_index = rSolv1n / r_step_size + 0.5_dp ! forces come from half bin position.
	ir1 = floor(f_index) ! get flanking r indicies
	ir2 = ir1 + 1

	float_index = cosTh1 / cosTheta_step_size + 0.5_dp
	ic1 = floor(f_index) ! get flanking cosTh indicies
	ic2 = ic1 + 1

	f11 = hist_array(ir1, ic1)
	f12 = hist_array(ir1, ic2)
	f21 = hist_array(ir2, ic1)
	f22 = hist_array(ir2, ic2)

	! fixme: define r1, r2, c1, and c2
	ra = r2-rSolv1n
	rb = rSolv1n-r1
	rd = r2-r1

	!fR1 = ra/rd*hist_array(ir1,ic1) + rb/rd*hist_array(ir2,ic1)
	!fR2 = ra/rd*hist_array(ir1,ic2) + rb/rd*hist_array(ir2,ic2)
	cd = c2-c1
	bi_out = (c2-c)/cd*(ra/rd*hist_array(ir1,ic1) + rb/rd*hist_array(ir2,ic1)) + (c-c1)/cd*(ra/rd*hist_array(ir1,ic2) + &
		& rb/rd*hist_array(ir2,ic2))

endsubroutine bilin_interpolate


! integrate the average force from 'compute_avg_force' to get the PMF.
subroutine integrate_force
	use histData
	use cfgData
	implicit none
	integer 		:: d

	allocate( u_dir(num_R_bins) )
	u_dir = 0_dp

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
	use histData
	implicit none
	character(len=64) 	:: outFile
	integer 			:: r

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
	use histData
	use cfgData
	use ctrlData
	implicit none
	integer				:: r, i, j, num_x_bins, num_z_bins
	real(kind=dp)		:: norm_const
	character(len=32)	:: temp, filename
	character(len=8)	:: frmt

	frmt = '(I3.3)' ! an integer of width 3 with zeroes on the left
	write(temp,frmt) r ! converting integer to string using 'internal file'
	filename='ellipse_output.'//trim(temp)//'.dat'

	norm_const = 1_dp / real(nThetaBins, dp) / real(phi_bins, dp)

	open(35,file=filename)
	write(6,*) "Writing test file:	", filename
	write(35,*) "# 1.	X Distance"
	write(35,*) "# 2.	Z Distance"
	write(35,*) "# 3.	Avg Force List"
	write(35,*) "# 4.	g(r) List"
	do i = 1, num_x_bins
		do j = 1, num_z_bins
			write(35,898) x_axis(i), z_axis(j), linAvg(i,j)*norm_const, grSPA(i,j)*norm_const
		enddo
	enddo
	close(35)
	
898		format (4(1x,f16.12))

endsubroutine write_test_out
