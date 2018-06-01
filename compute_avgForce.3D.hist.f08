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

! data for the density and force tables.
module histData
	use prec
	real(kind=dp), allocatable	:: histDist(:), histCosTh(:), histPhi(:), hist3D(:,:,:,:)
	real(kind=dp)				:: histDistStepSize, histCosThStepSize, histPhiStepSize
	integer						:: histDistBins, histCosThBins, histPhiBins

end module histData

! data from the config file.
module cfgData
	use prec
	real(kind=dp), allocatable	:: x_axis(:), z_axis(:), R_axis(:), fAvg(:), u_dir(:)
	real(kind=dp)				:: RStepSize, xzStepSize, R_min, R_max, xz_range, cosTh_max, cosTh_min, phi_max, phi_min, &
		cfgCosThStepSize, cfgPhiStepSize, cfgPsiStepSize
	character(len=8)			:: c_explicit_R
	integer						:: cfgRBins, cfgCosThBins, cfgPhiBins, cfgPsiBins

end module cfgData

! data for calculating cosTh value.
module angleData
	use prec
	real(kind=dp), allocatable	:: sinThetaLF(:), cosThetaLF(:), sinPhiLF(:), cosPhiLF(:), sinPsiLF(:), cosPsiLF(:)
	real(kind=dp)				:: rSolv1(3), rSolv2(3), rSolv1n, rSolv2n, cosTh1, cosTh2, phi1, phi2, psi1, psi2

end module angleData

! testing arrays for force and g(r)
module ctrlData
	use prec
	real(kind=dp), allocatable	:: frcSPA(:,:), grSPA(:,:), explicitDist(:)
	integer						:: crdLines
	logical						:: explicit_R

end module ctrlData


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
	
	! Write time taken to finish calculation.
	tf = omp_get_wtime()
	write(*,*) "Total time elapsed: ", tf-ti, "seconds"

end program compute_avgForce


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
			write(*,*) "hist File: ", histFile
		case ('-cfg')
			i = i+1
			call get_command_argument(i,cfgFile)
			cfgFileFlag=.true.
			write(*,*) "cfg File: ", cfgFile
		case default
			write(*,*) 'Unrecognized command-line option: ', arg
			write(*,*) 'Usage: compute_avgForce.x -hist [hist file]'
			write(*,*) 'Usage: compute_avgForce.x -cfg [cfg file]'
			stop

		end select
		i = i+1
		if (i.ge.command_argument_count()) exit
	end do

	if (histFileFlag.eqv..false.) then
		write(*,*) "Must provide a hist file using command line argument -hist [hist file name]"
		stop
	end if

	if (cfgFileFlag.eqv..false.) then
		write(*,*) "Must provide a cfg file using command line argument -cfg [cfg file name]"
		stop
	end if

end subroutine parse_command_line


! read python cfg file for g(r) parameters
subroutine read_cfg(cfgFile, outFile)
	use cfgData
	implicit none
	character(len=64) 	:: cfgFile, outFile
	character(len=128) 	:: line
	character(len=32) 	:: firstWord, sep
	integer 			:: ios
	logical 			:: outFileFlag, RstepSizeFlag, xzStepSizeFlag, RmaxFlag, RminFlag, xzRangeFlag, thetaBinsFlag, phiBinsFlag,&
		psiBinsFlag, c_explicit_RFlag

	outFileFlag = .false.
	RstepSizeFlag = .false.
	xzStepSizeFlag = .false.
	c_explicit_RFlag = .false.
	RmaxFlag = .false.
	RminFlag = .false.
	xzRangeFlag = .false.
	thetaBinsFlag = .false.
	phiBinsFlag = .false.
	psiBinsFlag = .false.

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
			else if (firstWord .eq. "RStepSize") then
				read(line,*) RStepSize
				write(*,*) "PMF Step Size:	", RStepSize
				RstepSizeFlag = .true.
			else if (firstWord .eq. "xzStepSize") then
				read(line,*) xzStepSize
				write(*,*) "Solvent Grid Step Size:	", xzStepSize
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
				read(line,*) cfgCosThBins
				write(*,*) "Theta Bins:	", cfgCosThBins
				thetaBinsFlag= .true.
			else if (firstWord .eq. "phi_bins") then
				read(line,*) cfgPhiBins
				write(*,*) "Phi Bins:	", cfgPhiBins
				phiBinsFlag= .true.
			else if (firstWord .eq. "psi_bins") then
				read(line,*) cfgPsiBins
				write(*,*) "Psi Bins:	", cfgPsiBins
				psiBinsFlag= .true.
			end if
		end if
	end do
	close(20)

	if (RstepSizeFlag.eqv..false.) then
		write(*,*) "Config file must have a 'RStepSize' value"
		stop
	end if
	if (xzStepSizeFlag.eqv..false.) then
		write(*,*) "Config file must have a 'xzStepSize' value"
		stop
	end if
	if (outFileFlag.eqv..false.) then
		write(*,*) "Config file must have a 'out_file' value"
		stop
	end if
	if (c_explicit_RFlag.eqv..false.) then
		write(*,*) "Config file must have a 'explicit_R' value"
		stop
	end if
	if (RmaxFlag.eqv..false.) then
		write(*,*) "Config file must have a 'R_max' value"
		stop
	end if
	if (RminFlag.eqv..false.) then
		write(*,*) "Config file must have a 'R_min' value"
		stop
	end if
	if (xzRangeFlag.eqv..false.) then
		write(*,*) "Config file must have a 'xz_range' value"
		stop
	end if
	if (thetaBinsFlag.eqv..false.) then
		write(*,*) "Config file must have a 'theta_bins' value"
		stop
	end if
	if (phiBinsFlag.eqv..false.) then
		write(*,*) "Config file must have a 'phi_bins' value"
		stop
	end if
	if (psiBinsFlag.eqv..false.) then
		write(*,*) "Config file must have a 'psi_bins' value"
		stop
	end if

end subroutine read_cfg


! read force file and make a lookup table.
subroutine make_hist_table(histFile)
	use histData
	implicit none
	character(len=64)			:: histFile
	character(len=32)			:: junk
	character(len=256)			:: line
	integer						:: ios, ios2, i, j, k, nHistLines
	real(kind=dp), allocatable	:: histTmp(:,:)

	! read number of lines in histFile and allocate that many points in temporary histogram list, histTmp.
	ios = 0; nHistLines = -1
	open(20,file=histFile)
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		if (line(1:1) .ne. "#") then
			nHistLines = nHistLines + 1
		end if
	end do
	close(20)

	allocate( histTmp(5,nHistLines) )

	! populate hist arrays
	ios = 0
	i = 1
	open(20,file=histFile)
	! read file ignoring comment lines at the beginning
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		if (line(1:1) .ne. "#") then
			!			 dist		   cos(Th)		 phi/3			g(r)+		 g(r)- force(r)+
			read(line,*) histTmp(1,i), histTmp(2,i), histTmp(3,i), histTmp(4,i), junk, histTmp(5,i)
			i = i + 1
		end if
	end do
	close(20)

	! XXX: Unique value determination -- TESTING
	do i = 1, nHistLines
		if (i .eq. 1) then
			histDistBins = 1
			ios = 0; ios2 = 0
		else ! i = 2, nHistLines
			if (( histTmp(1,i) .lt. (histTmp(1,i-1)-1d-6) ) .or. ( histTmp(1,i) .gt. (histTmp(1,i-1)+1d-6) )) then
				! note: this statement will trigger when a value in the first column (dist) is different than the value in the row
				! before it.
				histDistBins = histDistBins + 1
			end if
			if (( histTmp(2,i) .gt. (histTmp(2,1)-1d-6) ) .and. ( histTmp(2,i) .lt. (histTmp(2,1)+1d-6) ) .and. ( ios .eq. 0 ) &
				.and. ( ios2 .eq. 1 )) then
				! note: this statement will trigger when i = histCosThBins+1 because it finds the first repeated element
				histCosThBins = (i - 1)/histPhiBins
				ios = 1
			end if
			if (( histTmp(3,i) .gt. (histTmp(3,1)-1d-6) ) .and. ( histTmp(3,i) .lt. (histTmp(3,1)+1d-6) ) .and. ( ios2 .eq. 0 )) &
				then
				! note: this statement will trigger when i = histPhiBins+1 because it finds the first repeated element
				histPhiBins = i - 1
				ios2 = 1
			end if
		end if
	end do
	write(*,*) "Histogram Distance Bins: ", histDistBins
	write(*,*) "Histogram Cosine Theta Bins: ", histCosThBins
	write(*,*) "Histogram Phi Bins: ", histPhiBins
	
	allocate( histDist(histDistBins), histCosTh(histCosThBins), histPhi(histPhiBins), &
		hist3D(2,histDistBins,histCosThBins,histPhiBins) )
	
	! populate arrays that will be used in the rest of the calculation from temp array
	do i = 1, histDistBins		! the values written out from python script are at half-bin distances
		histDist(i) = histTmp(1,histCosThBins*i-1)
	end do
	do i = 1, histCosThBins
		histCosTh(i) = histTmp(2,i)
	end do
	do i = 1, histPhiBins
		histPhi(i) = histTmp(3,i)
	end do

	do i = 1, histDistBins
		do j = 1, histCosThBins
			do k = 1, histPhiBins
				hist3D(1,i,j,k) = histTmp(4, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k) ! g(r,cos,phi)
				hist3D(2,i,j,k) = histTmp(5, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k) ! frc(r,cos,phi)
			end do
		end do
	end do

	histDistStepSize = histDist(2) - histDist(1)
	write(*,*) "Histogram Distance Step Size: ", histDistStepSize
	histCosThStepSize = histCosTh(2) - histCosTh(1)
	write(*,*) "Histogram Cosine Theta Step Size: ", histCosThStepSize
	histPhiStepSize = histPhi(2) - histPhi(1)
	write(*,*) "Histogram Phi Step Size: ", histPhiStepSize 

end subroutine make_hist_table


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
		end do
		close(20)

		allocate( explicitDist(crdLines) )

		ios = 0
		open(20,file='crd_list.out',status='old')
		do i = 1, crdLines
			read(20,*,iostat=ios) junk, explicitDist(i)
		end do
		close(20)
	end if

end subroutine R_list


! do the average force integral
subroutine compute_avg_force
	use cfgData
	use angleData
	use ctrlData
	implicit none
	integer 		:: xBins, zBins, r, i, j, ithLF, iphiLF, ipsiLF
	real(kind=dp)	:: pi, density, phiLF, psiLF, psi_max, psi_min, gx, fx, fNew

	write(*,*) "Setting up for average force iteration..."
	pi = 3.1415926535_dp
	density = 0.00750924_dp

	if (explicit_R .eqv. .true.) then
		cfgRBins = crdLines
		write(*,*) "Number of R Bins: ", cfgRBins
	else if (explicit_R .eqv. .false.) then
		cfgRBins = int( (R_max - R_min)/RStepSize )
		write(*,*) "Number of R Bins: ", cfgRBins
	end if
	xBins = int( (2_dp * xz_range)/xzStepSize )
	write(*,*) "Number of X Bins: ", xBins
	zBins = int( (xz_range)/xzStepSize )
	write(*,*) "Number of Z Bins: ", zBins

	! allocate array sizes for axes and average force
	allocate( R_axis(cfgRBins), fAvg(cfgRBins), x_axis(xBins), z_axis(zBins) )
	R_axis = 0_dp; fAvg = 0_dp; x_axis = 0_dp; z_axis = 0_dp

	! allocate arrays for control arrays
	allocate( frcSPA(xBins, zBins), grSPA(xBins, zBins) )

	! Distance Axes
	do r = 1, cfgRBins
		if (explicit_R .eqv. .true.) then
			R_axis(r) = explicitDist(r)
		else if (explicit_r .eqv. .false.) then
			R_axis(r) = (r-1) * RStepSize + R_min
		end if
	end do
	do r = 1, crdLines
		R_axis(r) = explicitDist(r)
	end do
	do i = 1, xBins
		x_axis(i) = (i-1) * xzStepSize - xz_range + xzStepSize/2_dp
	end do
	do j = 1, zBins
		z_axis(j) = (j-1) * xzStepSize + xzStepSize/2_dp
	end do

	! Angles
	allocate( cosThetaLF(cfgCosThBins), sinThetaLF(cfgCosThBins), sinPhiLF(cfgPhiBins), cosPhiLF(cfgPhiBins), &
		sinPsiLF(cfgPsiBins), cosPsiLF(cfgPsiBins) )

	cosTh_max = 1_dp
	cosTh_min = -1_dp
	cfgCosThStepSize = (cosTh_max - cosTh_min) / real(cfgCosThBins, dp)
	do ithLF = 1, cfgCosThBins
		cosThetaLF(ithLF) = (ithLF-0.5_dp)*cfgCosThStepSize - 1_dp
		sinThetaLF(ithLF) = dsqrt(abs(1_dp-cosThetaLF(ithLF)**2))
	end do
	write(*,*) "Config Cos(Theta) Step Size: ", cfgCosThStepSize

	phi_max = 2_dp*pi/3_dp ! FIXME: this should be sufficient for a molecule with C3 symmetry.
	phi_min = 0_dp
	cfgPhiStepSize = (phi_max - phi_min) / real(cfgPhiBins, dp)
	do iphiLF = 1, cfgPhiBins
		phiLF = (iphiLF+0.5_dp)*cfgPhiStepSize
		sinPhiLF(iphiLF) = dsin(phiLF)
		cosPhiLF(iphiLF) = dcos(phiLF)
	end do
	write(*,*) "Config Phi Step Size: ", cfgPhiStepSize

	psi_max = 2_dp*pi
	psi_min = 0_dp
	cfgPsiStepSize = (psi_max - psi_min) / real(cfgPsiBins, dp)
	do ipsiLF = 1, cfgPsiBins
		psiLF = (ipsiLF+0.5_dp)*cfgPsiStepSize
		sinPsiLF(ipsiLF) = dsin(psiLF)
		cosPsiLF(ipsiLF) = dcos(psiLF)
	end do
	write(*,*) "Config Psi Step Size: ", cfgPsiStepSize

	write(*,*) "Computing average force..."
	! Calculate the average force integral for top half of bisecting plane of cylinder
	do r = 1, cfgRBins ! loop lj--lj distances
		frcSPA = 0_dp; grSPA = 0_dp
		do i = 1, xBins ! full length of cylinder
			do j = 1, zBins ! top half of bisecting plane of cylinder
				rSolv1(1) = x_axis(i)+R_axis(r)/2_dp
				rSolv1(2) = 0_dp
				rSolv1(3) = z_axis(j)
				rSolv1n = norm2(rSolv1)

				rSolv2(1) = x_axis(i)-R_axis(r)/2_dp
				rSolv2(2) = 0_dp
				rSolv2(3) = z_axis(j)
				rSolv2n = norm2(rSolv2)

				! loop through orientations of solvent at x(i) and z(j)
				do ithLF = 1, cfgCosThBins
					do iphiLF = 1, cfgPhiBins
						do ipsiLF = 1, cfgPsiBins
							if ((rSolv1n .lt. 1d-6) .or. (rSolv2n .lt. 1d-6)) then
								gx = 0_dp ! avoid NaNs in calc_angles
							else
								call calc_angles(ipsiLF, ithLF, iphiLF)
								call trilin_interpolate(1, gx) ! 1 is the index for the g(r,cos,phi)
							end if

							! if gx == 0 then don't waste time with the rest of the calculation
							if (gx .gt. 1d-6) then
								call trilin_interpolate(2, fx) ! 2 is the index for the frc(r,cos,phi)
								!		'gx' is the g(r) value taken from a histogram of g(r,cosTh,phi)
								! 		'fx' is ||fs(x,z)||, force from solvent at (x,z). Now we
								!		need cos(theta) and z and we should have the integral.
								fNew = ( gx * fx * z_axis(j) * ( rSolv1(1) / rSolv1n ) )
								fAvg(r) = fAvg(r) + fNew
								frcSPA(i,j) = frcSPA(i,j) + fNew
								grSPA(i,j) = grSPA(i,j) + gx
							end if
						end do !psi
					end do !phi
				end do !theta
			end do !z
		end do !x
		do i = 1, xBins
			do j = 1, zBins
				!frcSPA(i,j) = frcSPA(i,j) / 2_dp / pi / z_axis(j) ! get rid of jacobian that puts force field to 0 at z=0
				! normalize
				frcSPA(i,j) = frcSPA(i,j)/2_dp*density*xzStepSize*xzStepSize*cfgCosThStepSize*cfgPhiStepSize*cfgPsiStepSize
				grSPA(i,j) = grSPA(i,j)/cfgCosThBins/cfgPhiBins/cfgPsiBins
			end do !z again
		end do !x again
		call write_test_out(r, xBins, zBins) ! write grSPA and fx arrays

		! NOTE : After the fact multiply all elements by 2*pi*density/4/pi (4pi steradians from orientations)
		! 		Number density of chloroform per Angstrom**3 == 0.00750924
		fAvg(r) = fAvg(r)/2_dp*density*xzStepSize*xzStepSize*cfgCosThStepSize*cfgPhiStepSize*cfgPsiStepSize
	end do !r

end subroutine compute_avg_force


! rotate two solvent vectors 'h' for the dipole and 'l' for the Cl1 via a twist 'phi', tilt 'theta', and procession about z 'psi'
! for lj particles 1 and 2.
subroutine calc_angles(ipsiLF, ithLF, iphiLF)
	use cfgData
	use angleData
	use functions
	implicit none
	integer						:: iphiLF, ithLF, ipsiLF
	real(kind=dp),dimension(3)	:: h, l, na1, nb1, na2, nb2

	! make rotated solvent dipole vector at origin
	h(1) = sinPsiLF(ipsiLF)*sinThetaLF(ithLF)
	h(2) = -cosPsiLF(ipsiLF)*sinThetaLF(ithLF)
	h(3) = cosThetaLF(ithLF)

	! calculate cos(theta1) and cos(theta2) of the solvent to lj-spheres 1 and 2 respectively.
	cosTh1 = dot_product(rSolv1, h) / rSolv1n
	cosTh2 = dot_product(rSolv2, h) / rSolv2n

	! make rotated vector that represents the x-axis projection of one of the Cl vectors.
	l(1) = cosPhiLF(iphiLF)*cosPsiLF(ipsiLF) - sinPhiLF(iphiLF)*cosThetaLF(ithLF)*sinPsiLF(ipsiLF)
	l(2) = cosPhiLF(iphiLF)*sinPsiLF(ipsiLF) + sinPhiLF(iphiLF)*cosThetaLF(ithLF)*cosPsiLF(ipsiLF)
	l(3) = sinPhiLF(iphiLF)*sinThetaLF(ithLF)

	! calculate plane-normal vectors for LJ-C-H and LJ-C-Cl1
	! note: use the cross product of the lj-c (rsolv1) and c-h (h) vectors, and the lj-c (rsolv1) and c-cl (l) vectors
	na1 = cross_product(rSolv1,h)
	nb1 = cross_product(rSolv1,l)
	na2 = cross_product(rSolv2,h)
	nb2 = cross_product(rSolv2,l)

	! calculate phi1 and phi2 of the solvent to lj-spheres 1 and 2 respectively.
	phi1 = acos( dot_product(na1,nb1) / ( norm2(na1)*norm2(nb1) ) )
	phi2 = acos( dot_product(na2,nb2) / ( norm2(na2)*norm2(nb2) ) )

	if (phi1 .ge. phi_max) then
		phi1 = phi1 - phi_max
	end if
	if (phi2 .ge. phi_max) then
		phi2 = phi2 - phi_max
	end if

end subroutine calc_angles


! trilinearly interpolate
subroutine trilin_interpolate(ihist, fP)
	use cfgData
	use histData
	use angleData
	implicit none
	integer			:: i, iloop, ihist, ir0, ir1, ic0, ic1, ip0, ip1
	real(kind=dp)	:: rSolv(2), cosTh(2), phi(2), f_index, r0, r1, c0, c1, p0, p1, rd, cd, pd, f00, f01, f10, f11, f0, f1, fP

	rSolv(1) = rSolv1n
	rSolv(2) = rSolv2n
	cosTh(1) = cosTh1
	cosTh(2) = cosTh2
	phi(1) = phi1
	phi(2) = phi2
	fP = 1_dp
	if (ihist .eq. 1) then ! for g(r) both particles need to be taken into account
		iloop = 2
	else if (ihist .eq. 2) then ! for frc(r) only particle 1 needs to be accounted for
		iloop = 1
	end if
	do i = 1, iloop
		! Note: get the bounding values in each dimension.
		! r
		f_index = (rSolv(i) - histDistStepSize) / histDistStepSize + 1.5_dp ! take into account half-bin positions
		ir0 = floor(f_index) ! get flanking r indicies
		if (ir0 .ge. histDistBins) then
			ir0 = histDistBins - 1 ! put larger solvent positon values in the last bin of the histogram
		else if (ir0 .lt. 1) then
			ir0 = 1
		end if
		ir1 = ir0 + 1
		r0 = histDist(ir0)
		r1 = histDist(ir1)
		! FIXME: there might be an issue here with only the ith element of rSolv() getting affected. or maybe it ends up working out
		! with g(r) needing both so i->1-2 but force(r) only needing one so i->1
		if (rSolv(i) .gt. r1) then
			rSolv(i) = r1 ! if the solvent was far away enough to get repositioned to the last bin, set the distance to the last bin
		end if

		! cosTheta
		f_index = (cosTh(i) - cosTh_min) / histCosThStepSize + 1_dp ! so index values start at 1
		ic0 = floor(f_index) ! get flanking cosTh indicies
		if (ic0 .ge. histCosThBins) then
			ic0 = histCosThBins - 1
		else if (ic0 .lt. 1) then
			ic0 = 1
		end if
		ic1 = ic0 + 1
		c0 = histCosTh(ic0)
		c1 = histCosTh(ic1)

		! phi
		f_index = (phi(i) - phi_min) / histPhiStepSize + 1_dp ! so index values start at 1
		ip0 = floor(f_index) ! get flanking phi indicies
		if (ip0 .ge. histPhiBins) then
			ip0 = histPhiBins - 1
		else if (ip0 .lt. 1) then
			ip0 = 1
		end if
		ip1 = ip0 + 1
		p0 = histPhi(ip0)
		p1 = histPhi(ip1)

		! Note: interpolate
		rd = (rSolv(i)-r0)/(r1-r0)
		cd = (cosTh(i)-c0)/(c1-c0)
		pd = (phi(i)-p0)/(p1-p0)

		f00 = hist3D(ihist,ir0,ic0,ip0)*(1-rd) + hist3D(ihist,ir1,ic0,ip0)*rd
		f01 = hist3D(ihist,ir0,ic0,ip1)*(1-rd) + hist3D(ihist,ir1,ic0,ip1)*rd
		f10 = hist3D(ihist,ir0,ic1,ip0)*(1-rd) + hist3D(ihist,ir1,ic1,ip0)*rd
		f11 = hist3D(ihist,ir0,ic1,ip1)*(1-rd) + hist3D(ihist,ir1,ic1,ip1)*rd

		f0 = f00*(1-cd) + f10*cd
		f1 = f01*(1-cd) + f11*cd
		
		fP = f0*(1-pd) + f1*pd
	end do

end subroutine trilin_interpolate


! integrate the average force from 'compute_avg_force' to get the PMF.
subroutine integrate_force
	use cfgData
	implicit none
	integer 		:: d

	allocate( u_dir(cfgRBins) )
	u_dir = 0_dp

	do d = 1, cfgRBins
		if (d .eq. 1) then
			u_dir(cfgRBins) = fAvg(cfgRBins) * RStepSize
		else
			u_dir(cfgRBins-(d-1)) = u_dir(cfgRBins-(d-2)) + fAvg(cfgRBins-(d-1)) * RStepSize
		end if
	end do

end subroutine integrate_force


! write force out and g(r) out to compare against explicit
subroutine write_test_out(i_f, xBins, zBins)
	use cfgData
	use ctrlData
	implicit none
	integer				:: i_f, i, j, xBins, zBins
	character(len=32)	:: temp, filename
	character(len=8)	:: frmt

	frmt = '(I3.3)' ! an integer of width 3 with zeroes on the left
	write(temp,frmt) i_f ! converting integer to string using 'internal file'
	filename='hist3D_output.'//trim(temp)//'.dat'


	open(35,file=filename,status='replace')
	write(6,*) "Writing test file:	", filename
	write(35,*) "# 1.	X Distance"
	write(35,*) "# 2.	Z Distance"
	write(35,*) "# 3.	g(r)"
	write(35,*) "# 4.	Force"
	do i = 1, xBins
		do j = 1, zBins
			write(35,898) x_axis(i), z_axis(j), grSPA(i,j), frcSPA(i,j)
		end do
	end do
	close(35)
	
898		format (4(1x,es14.7))

end subroutine write_test_out


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
	do r = 1, cfgRBins
		write(35,899) R_axis(r), fAvg(r), u_dir(r)
	end do
	close(35)

899		format (3(1x,es14.7)) ! scientific format

end subroutine write_output
