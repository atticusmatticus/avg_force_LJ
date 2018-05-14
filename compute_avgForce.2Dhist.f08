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
!		-x		1				2	+x
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
	real(kind=dp), allocatable	:: histDist(:), histAng(:), hist2D(:,:,:)
	real(kind=dp)				:: histDist_step_size, histCosTh_step_size, hHistDist_step_size
	integer						:: histDistBins, histCosThBins

endmodule histData

! data from the config file.
module cfgData
	use prec
	real(kind=dp), allocatable	:: x_axis(:), z_axis(:), R_axis(:), fAvg(:), u_dir(:)
	real(kind=dp)				:: R_step_size, xz_step_size, R_min, R_max, xz_range, cosTh_max, cosTh_min, cfgCosTh_step_size, &
		phi_step_size
	character(len=8)			:: c_explicit_R
	integer						:: num_R_bins, cfgCosTh_bins, cfgPhiBins

endmodule cfgData

! data for calculating cosTh value.
module thetaData
	use prec
	real(kind=dp), allocatable	:: cosThetaLF(:), sinThetaLF(:), sinPhiLF(:), cosPhiLF(:)
	real(kind=dp)				:: rSolv1(3), rSolv2(3), rSolv1n, rSolv2n, cosTh1, cosTh2

endmodule thetaData

! testing arrays for force and g(r)
module ctrlData
	use prec
	real(kind=dp), allocatable	:: frcSPA(:,:), grSPA(:,:), explicitDist(:)
	integer						:: crdLines
	logical						:: explicit_R

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

	! read in LJ--LJ dist array from file
	call R_list

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
	logical 			:: outFileFlag, RstepSizeFlag, xzStepSizeFlag, RmaxFlag, RminFlag, xzRangeFlag, phiBinsFlag, thetaBinsFlag,&
		c_explicit_RFlag

	outFileFlag = .false.
	RstepSizeFlag = .false.
	xzStepSizeFlag = .false.
	c_explicit_RFlag = .false.
	RmaxFlag = .false.
	RminFlag = .false.
	xzRangeFlag = .false.
	phiBinsFlag = .false.
	thetaBinsFlag = .false.

	ios = 0

	open(20,file=cfgFile)
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		call split(line,'=',firstWord, sep)
		if (line .ne. "") then
			if (firstWord .eq. "out_file") then
				read(line,*) outFile
				write(*,*) "Output File:	", outFile
				outFileFlag = .true.
			else if (firstWord .eq. "R_step_size") then
				read(line,*) R_step_size
				write(*,*) "PMF Step Size:	", R_step_size
				RstepSizeFlag = .true.
			else if (firstWord .eq. "xz_step_size") then
				read(line,*) xz_step_size
				write(*,*) "Solvent Grid Step Size:	", xz_step_size
				xzStepSizeFlag = .true.
			else if (firstWord .eq. "explicit_R") then
				read(line,*) c_explicit_R
				write(*,*) "Use Explicit R Values:		", c_explicit_R
				c_explicit_RFlag = .true.
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
			else if (firstWord .eq. "theta_bins") then
				read(line,*) cfgCosTh_bins
				write(*,*) "Theta Bins:	", cfgCosTh_bins
				thetaBinsFlag= .true.
			else if (firstWord .eq. "phi_bins") then
				read(line,*) cfgPhiBins
				write(*,*) "Phi Bins:	", cfgPhiBins
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
	if (c_explicit_RFlag.eqv..false.) then
		write(*,*) "Config file must have a 'explicit_R' value"
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
	if (thetaBinsFlag.eqv..false.) then
		write(*,*) "Config file must have a 'cfgCosTh_bins' value"
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
	integer						:: ios, i, j, nHistLines
	real(kind=dp), allocatable	:: histTmp(:,:)

	! read number of lines in histFile and allocate that many points in temporary histogram list.
	ios = 0; nHistLines = -1
	open(20,file=histFile)
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		if (line(1:1) .ne. "#") then
			nHistLines = nHistLines + 1
		endif
	enddo
	close(20)

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
			histDistBins = 1
			ios = 0
		else ! i = 2, nHistLines
			if ( (histTmp(1,i) .lt. (histTmp(1,i-1)-1d-6)) .or. (histTmp(1,i) .gt. (histTmp(1,i-1)+1d-6)) ) then
				histDistBins = histDistBins + 1
			endif
			if ( (histTmp(2,i) .gt. (histTmp(2,1)-1d-6)) .and. (histTmp(2,i) .lt. (histTmp(2,1)+1d-6)) .and. (ios .eq. 0) ) then
				! note: this statement will trigger when i = histCosThBins + 1 because it finds the first repeated element
				histCosThBins = i - 1
				ios = 1
			endif
		endif
	enddo
	print *, "Histogram Distance Bins: ", histDistBins
	print *, "Histogram Cosine Theta Bins: ", histCosThBins
	
	allocate( histDist(histDistBins), histAng(histCosThBins), hist2D(2,histDistBins,histCosThBins) )
	
	! populate arrays that will be used in the rest of the calculation from temp array
	do i = 1, histDistBins		! the values printed out from python script are at half-bin distances
		histDist(i) = histTmp(1,histCosThBins*i-1)
	enddo
	do i = 1, histCosThBins
		histAng(i) = histTmp(2,i)
	enddo

	do i = 1, histDistBins
		do j = 1, histCosThBins
			hist2D(1,i,j) = histTmp(3,(i-1)*histCosThBins+j) ! g(r,cos)
			hist2D(2,i,j) = histTmp(4,(i-1)*histCosThBins+j) ! frc(r,cos)
		enddo
	enddo

	histDist_step_size = histDist(2) - histDist(1)
	hHistDist_step_size = histDist_step_size / 2_dp
	write(*,*) "Histogram Distance Step Size: ", histDist_step_size
	histCosTh_step_size = histAng(2) - histAng(1)
	write(*,*) "Histogram Angle Step Size: ", histCosTh_step_size

endsubroutine make_hist_table


! read LJ--LJ displacements from file
subroutine R_list
	use cfgData
	use ctrlData
	implicit none
	integer				:: ios, i
	character(len=16)	:: junk
	character(len=64)	:: line

	if (c_explicit_R .eq. 'no') then
		explicit_R = .false.
	else if (c_explicit_R .eq. 'yes') then
		explicit_R = .true.

		ios = 0; crdLines = -1
		open(20,file='crd_list.out',status='old')
		do while(ios>=0)
			read(20,'(a)',IOSTAT=ios) line
			crdLines = crdLines + 1
		enddo
		close(20)

		allocate( explicitDist(crdLines) )

		ios = 0
		open(20,file='crd_list.out',status='old')
		do i = 1, crdLines
			read(20,*,iostat=ios) junk, explicitDist(i)
		enddo
		close(20)
	endif

endsubroutine R_list


! do the average force integral
subroutine compute_avg_force
	use cfgData
	use thetaData
	use ctrlData
	implicit none
	integer 		:: num_x_bins, num_z_bins, r, i, j, ithLF, iphiLF
	real(kind=dp)	:: pi, phiLF, phi_max, phi_min, gx, fx, fNew

	pi = 3.1415926535_dp

	if (explicit_R .eqv. .true.) then
		num_R_bins = crdLines
		write(*,*) "Number of R Bins: ", num_R_bins
	else if (explicit_R .eqv. .false.) then
		num_R_bins = int( (R_max - R_min)/R_step_size )
		write(*,*) "Number of R Bins: ", num_R_bins
	endif
	num_x_bins = int( (2_dp * xz_range)/xz_step_size )
	write(*,*) "Number of X Bins: ", num_x_bins
	num_z_bins = int( (xz_range)/xz_step_size )
	write(*,*) "Number of Z Bins: ", num_z_bins

	! allocate array sizes for axes and average force
	allocate( R_axis(num_R_bins), fAvg(num_R_bins), x_axis(num_x_bins), z_axis(num_z_bins) )
	R_axis = 0_dp; fAvg = 0_dp; x_axis = 0_dp; z_axis = 0_dp

	! allocate arrays for control arrays
	allocate( frcSPA(num_x_bins, num_z_bins), grSPA(num_x_bins, num_z_bins) )

	! Distance Axes
	do r = 1, num_R_bins
		if (explicit_R .eqv. .true.) then
			R_axis(r) = explicitDist(r)
		else if (explicit_r .eqv. .false.) then
			R_axis(r) = (r-1) * R_step_size + R_min
		endif
	enddo
	do r = 1, crdLines
		R_axis(r) = explicitDist(r)
	enddo
	do i = 1, num_x_bins
		x_axis(i) = (i-1) * xz_step_size - xz_range + xz_step_size/2_dp
	enddo
	do j = 1, num_z_bins
		z_axis(j) = (j-1) * xz_step_size + xz_step_size/2_dp
	enddo

	! Angles
	allocate( cosThetaLF(cfgCosTh_bins), sinThetaLF(cfgCosTh_bins), sinPhiLF(cfgPhiBins), cosPhiLF(cfgPhiBins) )
	cosTh_max = 1_dp
	cosTh_min = -1_dp
	cfgCosTh_step_size = (cosTh_max - cosTh_min) / real(cfgCosTh_bins, dp)
	do ithLF = 1, cfgCosTh_bins
		cosThetaLF(ithLF) = (ithLF-0.5_dp)*cfgCosTh_step_size - 1_dp
		sinThetaLF(ithLF) = dsqrt(abs(1_dp-cosThetaLF(ithLF)**2))
	enddo
	write(*,*) "Cos(Theta) Step Size: ", cfgCosTh_step_size
	phi_max = 2_dp*pi
	phi_min = 0_dp
	phi_step_size = (phi_max - phi_min) / real(cfgPhiBins, dp)
	do iphiLF = 1, cfgPhiBins
		phiLF = (iphiLF+0.5_dp)*phi_step_size
		sinPhiLF(iphiLF) = dsin(phiLF)
		cosPhiLF(iphiLF) = dcos(phiLF)
	enddo
	write(*,*) "Phi Step Size: ", phi_step_size

	! Calculate the average force integral for top half of bisecting plane of cylinder
	do r = 1, num_R_bins ! loop lj--lj distances
		frcSPA = 0_dp; grSPA = 0_dp
		do i = 1, num_x_bins ! full length of cylinder
			do j = 1, num_z_bins ! top half of bisecting plane of cylinder
				rSolv1(1) = x_axis(i)+R_axis(r)/2_dp
				rSolv1(2) = 0_dp
				rSolv1(3) = z_axis(j)
				rSolv1n = norm2(rSolv1)

				rSolv2(1) = x_axis(i)-R_axis(r)/2_dp
				rSolv2(2) = 0_dp
				rSolv2(3) = z_axis(j)
				rSolv2n = norm2(rSolv2)

				! loop through orientations of solvent at x(i) and z(j)
				do ithLF = 1, cfgCosTh_bins
					do iphiLF = 1, cfgPhiBins
						if ((rSolv1n .lt. 1d-6) .or. (rSolv2n .lt. 1d-6)) then
							gx = 0_dp ! avoid NaNs in calc_cosTh
						else
							call calc_cosTh(iphiLF, ithLF)
							call bilin_interpolate(1, gx) ! 1 is the index for the g(r,cos)
						endif

						! if gx == 0 then don't waste time with the rest of the calculation
						if (gx .gt. 1d-6) then
							call bilin_interpolate(2, fx) ! 2 is the index for the frc(r,cos)
							!		'gx' is the g(r) value taken from a histogram of g(r,cosTh)
							! 		'fx' is ||fs(x,z)||, force from solvent at (x,z). Now we
							!		need cos(theta) and z and we should have the integral.
							fNew = ( gx * fx * z_axis(j) * ( rSolv1(1) / rSolv1n ) )
							fAvg(r) = fAvg(r) + fNew
							frcSPA(i,j) = frcSPA(i,j) + fNew
							grSPA(i,j) = grSPA(i,j) + gx
						endif
					enddo !phi
				enddo !theta
			enddo !z
		enddo !x
		do i = 1, num_x_bins
			do j = 1, num_z_bins
				!frcSPA(i,j) = frcSPA(i,j) / 2_dp / pi / z_axis(j) ! get rid of jacobian that puts force field to 0 at z=0
				! normalize
				frcSPA(i,j) = frcSPA(i,j)/2_dp*0.00750924_dp*xz_step_size*xz_step_size*phi_step_size*cfgCosTh_step_size
				grSPA(i,j) = grSPA(i,j)/cfgPhiBins/cfgCosTh_bins
			enddo !z again
		enddo !x again
		call write_test_out(r, num_x_bins, num_z_bins) ! write grSPA and fx arrays

		! NOTE : After the fact multiply all elements by 2*pi*density/4/pi (4pi steradians from orientations)
		! 		Number density of chloroform per Angstrom**3 == 0.00750924
		fAvg(r) = fAvg(r)/2_dp*0.00750924_dp*xz_step_size*xz_step_size*phi_step_size*cfgCosTh_step_size
	enddo !r

endsubroutine compute_avg_force


! rotate solvent vector 'p' and calculate cosTh for lj particles 1 and 2.
subroutine calc_cosTh(iphiLF, ithLF)
	use cfgData
	use thetaData
	implicit none
	integer			:: iphiLF, ithLF
	real(kind=dp)	:: p(3)

	! make rotated solvent vector at origin
	p(1) = sinPhiLF(iphiLF)*sinThetaLF(ithLF)
	p(2) = cosPhiLF(iphiLF)*sinThetaLF(ithLF)
	p(3) = cosThetaLF(ithLF)

	! calculate cos(theta1) and cos(theta2) relative to lj-spheres 1 and 2.
	cosTh1 = dot_product(rSolv1, p) / rSolv1n
	cosTh2 = dot_product(rSolv2, p) / rSolv2n

endsubroutine calc_cosTh


! bilinearly interpolate
subroutine bilin_interpolate(ihist, fP)
	use cfgData
	use histData
	use thetaData
	implicit none
	integer			:: i, iloop, ihist, ir1, ir2, ic1, ic2
	real(kind=dp)	:: rSolv(2), cosTh(2), f_index, r1, r2, c1, c2, ra, rb, fP

	rSolv(1) = rSolv1n
	rSolv(2) = rSolv2n
	cosTh(1) = cosTh1
	cosTh(2) = cosTh2
	fP = 1_dp
	if (ihist .eq. 1) then
		iloop = 2
	else if (ihist .eq. 2) then
		iloop = 1
	endif
	do i = 1, iloop
		f_index = (rSolv(i) - histDist_step_size) / histDist_step_size + 1.5_dp ! take into account half-bin positions
		ir1 = floor(f_index) ! get flanking r indicies
		if (ir1 .ge. histDistBins) then
			ir1 = histDistBins - 1 ! put larger solvent positon values in the last bin of the histogram
		else if (ir1 .lt. 1) then
			ir1 = 1
		endif
		ir2 = ir1 + 1
		r1 = histDist(ir1)
		r2 = histDist(ir2)

		f_index = (cosTh(i) - cosTh_min) / histCosTh_step_size + 1_dp ! so index values start at 1
		ic1 = floor(f_index) ! get flanking cosTh indicies
		if (ic1 .ge. histCosThBins) then
			ic1 = histCosThBins - 1
		else if (ic1 .lt. 1) then
			ic1 = 1
		endif
		ic2 = ic1 + 1
		c1 = histAng(ic1)
		c2 = histAng(ic2)

		if (rSolv(i) .gt. r2) then
			rSolv(i) = r2 ! if the solvent was far away enough to get repositioned to the last bin, set the distance to the last bin
		endif
		ra = r2 - rSolv(i)
		rb = rSolv(i) - r1

		fP = fP * ( &
			+ (c2-cosTh(i)) / (c2 - c1) * (ra / (r2 - r1) * hist2D(ihist,ir1,ic1) &
			+ rb / (r2 - r1) * hist2D(ihist,ir2,ic1)) &
			+ (cosTh(i)-c1) / (c2 - c1) * (ra / (r2 - r1) * hist2D(ihist,ir1,ic2) &
			+ rb / (r2 - r1) * hist2D(ihist,ir2,ic2)) )
	enddo

endsubroutine bilin_interpolate


! integrate the average force from 'compute_avg_force' to get the PMF.
subroutine integrate_force
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


! write force out and g(r) out to compare against explicit
subroutine write_test_out(i_f, num_x_bins, num_z_bins)
	use cfgData
	use ctrlData
	implicit none
	integer				:: i_f, i, j, num_x_bins, num_z_bins
	character(len=32)	:: temp, filename
	character(len=8)	:: frmt

	frmt = '(I3.3)' ! an integer of width 3 with zeroes on the left
	write(temp,frmt) i_f ! converting integer to string using 'internal file'
	filename='hist2D_output.'//trim(temp)//'.dat'


	open(35,file=filename,status='replace')
	write(6,*) "Writing test file:	", filename
	write(35,*) "# 1.	X Distance"
	write(35,*) "# 2.	Z Distance"
	write(35,*) "# 3.	Avg Force List"
	write(35,*) "# 4.	g(r) List"
	do i = 1, num_x_bins
		do j = 1, num_z_bins
			write(35,898) x_axis(i), z_axis(j), frcSPA(i,j), grSPA(i,j)
		enddo
	enddo
	close(35)
	
898		format (4(1x,es14.7))

endsubroutine write_test_out


! write output file
subroutine write_output(outFile)
	use cfgData
	implicit none
	character(len=64) 	:: outFile
	integer 			:: r

	open(35,file=outFile,status='replace')
	write(6,*) "Writing output file:	", outFile
	write(35,*) "# 1.	R Distance"
	write(35,*) "# 2.	Avg Force"
	write(35,*) "# 3.	PMF"
	do r = 1, num_R_bins
		write(35,899) R_axis(r), fAvg(r), u_dir(r)
	enddo
	close(35)

899		format (3(1x,es14.7)) ! scientific format

endsubroutine write_output
