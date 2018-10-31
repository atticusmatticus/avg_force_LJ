! USAGE: ./this_file.x -cfg [CONFIGURATION FILE] -hist [3D HIST FILE]
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
	real(kind=dp)				:: RStepSize, xzStepSize, R_min, R_max, xz_range, cosTh_max, cosTh_min, phi_max, phi_min, phi_hmax,&
		cfgCosThStepSize, cfgPhiStepSize, cfgPsiStepSize
	character(len=8)			:: c_explicit_R, g_fit
	integer						:: cfgRBins, cfgCosThBins, cfgPhiBins, cfgPsiBins

end module cfgData

! data for calculating cosTh value.
module angleData
	use prec
	real(kind=dp), allocatable	:: sinThetaLF(:), cosThetaLF(:), sinPhiLF(:), cosPhiLF(:), sinPsiLF(:), cosPsiLF(:)
	real(kind=dp)				:: rSolv1(3), rSolv2(3), rSolvn(2), cosTh(2), phi(2), sSolv1(3), tSolv1(3), sSolv1n, tSolv1n

end module angleData

! testing arrays for force and g(r)
module ctrlData
	use prec
	real(kind=dp), allocatable	:: frcSPA(:,:,:), grSPA(:,:), explicitDist(:)
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

	flush(6)

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

	flush(6)

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
		psiBinsFlag, c_explicit_RFlag, g_fitFlag

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
	g_fitFlag = .false.

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
			else if (firstWord .eq. "g_fit") then
				read(line,*) g_fit
				write(*,*) "Use g(r) fitting:		", g_fit
				g_fitFlag = .true.
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
	if (g_fitFlag.eqv..false.) then
		write(*,*) "Config file must have a 'g_fit' value"
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

	flush(6)

end subroutine read_cfg


! read force file and make a lookup table.
subroutine make_hist_table(histFile)
	use histData
	use cfgData
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
	print*, "nHistLines", nHistLines

	allocate( histTmp(7,nHistLines) )

	! populate hist arrays
	ios = 0; i = 1
	open(20,file=histFile)
	! read file ignoring comment lines at the beginning
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		if ((line(1:1) .ne. "#") .and. (ios .ge. 0)) then
			!			 dist		   cos(Th)		 phi/3			g(r)+		 g(r)-	<f.r>+			<f.s>+		<f.t>+
			read(line,*) histTmp(1,i), histTmp(2,i), histTmp(3,i), histTmp(4,i), junk, histTmp(5,i), histTmp(6,i), histTmp(7,i)
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
		hist3D(4,histDistBins,histCosThBins,histPhiBins) )
	
	! populate arrays that will be used in the rest of the calculation from temp array
	do i = 1, histDistBins		! the values written out from python script are at half-bin distances
		histDist(i) = histTmp(1,histCosThBins*histPhiBins*(i-1)+1)
	end do
	do i = 1, histCosThBins
		histCosTh(i) = histTmp(2,histPhiBins*(i-1)+1)
	end do
	do i = 1, histPhiBins
		histPhi(i) = histTmp(3,i)
	end do

	do i = 1, histDistBins
		do j = 1, histCosThBins
			do k = 1, histPhiBins
				if (g_fit .eq. 'log') then
					! Note: interpolate the log of g(r) because it will be smoother and better behaved facilitating linear fitting.
					if (histTmp(4, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k) .lt. 1d-6) then ! g(r,cos,phi)
						hist3D(1,i,j,k) = -1d12
					else
						hist3D(1,i,j,k) = dlog(histTmp(4, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k))
					end if
				else if (g_fit .eq. 'lin') then
					hist3D(1,i,j,k) = histTmp(4, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k) ! g(r,cos,phi)
				end if
				hist3D(2,i,j,k) = histTmp(5, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k) ! <f.r>(r,cos,phi)
				hist3D(3,i,j,k) = histTmp(6, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k) ! <f.s>(r,cos,phi)
				hist3D(4,i,j,k) = histTmp(7, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k) ! <f.t>(r,cos,phi)
			end do
		end do
	end do

	histDistStepSize = histDist(2) - histDist(1)
	write(*,*) "Histogram Distance Step Size: ", histDistStepSize
	histCosThStepSize = histCosTh(2) - histCosTh(1)
	write(*,*) "Histogram Cosine Theta Step Size: ", histCosThStepSize
	histPhiStepSize = histPhi(2) - histPhi(1)
	write(*,*) "Histogram Phi Step Size: ", histPhiStepSize 

	flush(6)

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
	use constants
	implicit none
	integer 		:: xBins, zBins, r, i, j, ithLF, iphiLF, ipsiLF
	real(kind=dp)	:: density, phiLF, psiLF, psi_max, psi_min, gx(3), fx(3), fNew1, fNew2, fNew3

	write(*,*) "Setting up for average force iteration..."
	density = 0.00750924_dp ! numerical density of chloroforms per Angstrom**3

	if (explicit_R .eqv. .true.) then
		cfgRBins = crdLines
		write(*,*) "Number of R Bins: ", cfgRBins
	else if (explicit_R .eqv. .false.) then
		cfgRBins = int( (R_max - R_min)/RStepSize )
		if (cfgRBins .eq. 0) then
			cfgRBins = 1
		end if
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
	allocate( frcSPA(3, xBins, zBins), grSPA(xBins, zBins) )

	! Distance Axes
	do r = 1, cfgRBins
		if (explicit_R .eqv. .true.) then
			R_axis(r) = explicitDist(r)
		else if (explicit_r .eqv. .false.) then
			R_axis(r) = (r-1) * RStepSize + R_min
		end if
	end do
	do i = 1, xBins
		x_axis(i) = (i-1) * xzStepSize - xz_range + xzStepSize/2_dp
	end do
	do j = 1, zBins
		z_axis(j) = (j-1) * xzStepSize + xzStepSize/2_dp
	end do

	! ANGLES
	! Theta
	! tilt off of z
	allocate( cosThetaLF(cfgCosThBins), sinThetaLF(cfgCosThBins), sinPhiLF(cfgPhiBins), cosPhiLF(cfgPhiBins), &
		sinPsiLF(cfgPsiBins), cosPsiLF(cfgPsiBins) )

	cosTh_max = 1_dp
	cosTh_min = -1_dp
	cfgCosThStepSize = (cosTh_max - cosTh_min) / real(cfgCosThBins, dp)
	do ithLF = 1, cfgCosThBins
		cosThetaLF(ithLF) = (ithLF-0.5_dp)*cfgCosThStepSize - cosTh_max
		sinThetaLF(ithLF) = dsqrt(abs(1_dp-cosThetaLF(ithLF)**2))
	end do
	write(*,*) "Config Cos(Theta) Step Size: ", cfgCosThStepSize

	! Phi
	! twist about z
	phi_max = 2_dp*pi/3_dp ! note: this is sufficient for a molecule with C3 symmetry.
	phi_hmax = pi/3_dp
	phi_min = 0_dp
	cfgPhiStepSize = (phi_max - phi_min) / real(cfgPhiBins, dp)
	do iphiLF = 1, cfgPhiBins
		phiLF = (iphiLF+0.5_dp)*cfgPhiStepSize
		sinPhiLF(iphiLF) = dsin(phiLF)
		cosPhiLF(iphiLF) = dcos(phiLF)
	end do
	write(*,*) "Config Phi Step Size: ", cfgPhiStepSize

	! Psi
	! processison about z
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

	flush(6)

	! Calculate the average force integral for top half of bisecting plane of cylinder
	do r = 1, cfgRBins ! loop lj--lj distances
		frcSPA = 0_dp; grSPA = 0_dp
		do i = 1, xBins ! full length of cylinder
			do j = 1, zBins ! top half of bisecting plane of cylinder
				rSolv1(1) = x_axis(i)+R_axis(r)/2_dp
				rSolv1(2) = 0_dp
				rSolv1(3) = z_axis(j)
				rSolvn(1) = norm2(rSolv1)

				rSolv2(1) = x_axis(i)-R_axis(r)/2_dp
				rSolv2(2) = 0_dp
				rSolv2(3) = z_axis(j)
				rSolvn(2) = norm2(rSolv2)

				! loop through orientations of solvent at x(i) and z(j)
				do ithLF = 1, cfgCosThBins
					do iphiLF = 1, cfgPhiBins
						do ipsiLF = 1, cfgPsiBins
							if ((rSolvn(1) .lt. 1d-6) .or. (rSolvn(2) .lt. 1d-6)) then
								gx = 0_dp ! avoid NaNs in calc_angles
							else
								call calc_angles(ipsiLF, ithLF, iphiLF)
								call trilin_interpolate(1, gx) ! 1 is the index for the g(r,theta,phi)
								if (g_fit .eq. 'log') then
									gx(1) = dexp(gx(1))
								end if
							end if

							if (gx(1) .gt. 1d-6) then ! if gx(1) == 0 then don't waste time with the rest of the calculation
								call trilin_interpolate(2, fx) ! 2 is the index for the frc(r,cos,phi)

								fNew1 = gx(1) * fx(1) * rSolv1(1)/rSolvn(1) ! f*g.x^{hat} for r
								fNew2 = gx(1) * fx(2) * (sSolv1(1)/sSolv1n) ! f*g.x^{hat} for s
								fNew3 = gx(1) * fx(3) * (tSolv1(1)/tSolv1n) ! f*g.x^{hat} for t

								fAvg(r) = fAvg(r) + ((fNew1 + fNew2 + fNew3) * z_axis(j)) ! note: added r here instead of in fNew

								frcSPA(1,i,j) = frcSPA(1,i,j) + fNew1
								frcSPA(2,i,j) = frcSPA(2,i,j) + fNew2
								frcSPA(3,i,j) = frcSPA(3,i,j) + fNew3
								grSPA(i,j) = grSPA(i,j) + gx(1)
							end if
						end do !psi
					end do !phi
				end do !theta
			end do !z
		end do !x
		! Normalize
		do i = 1, xBins
			do j = 1, zBins
				frcSPA(1,i,j) = frcSPA(1,i,j)/grSPA(i,j)
				frcSPA(2,i,j) = frcSPA(2,i,j)/grSPA(i,j)
				frcSPA(3,i,j) = frcSPA(3,i,j)/grSPA(i,j)
				grSPA(i,j) = grSPA(i,j)/cfgCosThBins/cfgPhiBins/cfgPsiBins
			end do !z again
		end do !x again
		call write_test_out(r, xBins, zBins) ! write grSPA and frcSPA arrays

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
	use constants
	implicit none
	integer						:: iphiLF, ithLF, ipsiLF
	real(kind=dp),dimension(3)	:: h, x, y

	! make rotated solvent dipole vector at origin
	h(1) = sinPsiLF(ipsiLF)*sinThetaLF(ithLF)
	h(2) = -cosPsiLF(ipsiLF)*sinThetaLF(ithLF)
	h(3) = cosThetaLF(ithLF)

	! calculate cos(theta1) and cos(theta2) of the solvent to lj-spheres 1 and 2 respectively.
	cosTh(1) = dot_product(rSolv1, h) / rSolvn(1)
	cosTh(2) = dot_product(rSolv2, h) / rSolvn(2)

	! make rotated vector that represents the x-axis projection of one of the Cl vectors.
	x(1) = cosPhiLF(iphiLF)*cosPsiLF(ipsiLF) - sinPhiLF(iphiLF)*cosThetaLF(ithLF)*sinPsiLF(ipsiLF)
	x(2) = cosPhiLF(iphiLF)*sinPsiLF(ipsiLF) + sinPhiLF(iphiLF)*cosThetaLF(ithLF)*cosPsiLF(ipsiLF)
	x(3) = sinPhiLF(iphiLF)*sinThetaLF(ithLF)

	! calculate plane-normal vectors for LJ-C-H and LJ-C-Cl1
	! note: use the cross product of the lj-c (rsolv1) and c-h (h) vectors, and the lj-c (rsolv1) and c-cl (l) vectors
	y = cross_product(h,x)

	phi(1) = datan2( dot_product(y,-rSolv1) / (norm2(y)*rSolvn(1)), dot_product(x,-rSolv1) / (norm2(y)*rSolvn(1)))
	phi(2) = datan2( dot_product(y,-rSolv2) / (norm2(y)*rSolvn(2)), dot_product(x,-rSolv2) / (norm2(y)*rSolvn(2)))

	! phi [-pi,pi] --> phi'' [0,pi/3]
	if ((phi(1) .gt. phi_hmax) .and. (phi(1) .lt. phi_max)) then
		phi(1) = phi_max - phi(1)
	else if ((phi(1) .gt. phi_max) .and. (phi(1) .lt. pi)) then
		phi(1) = phi(1) - phi_max
	else if ((phi(1) .gt. -pi) .and. (phi(1) .lt. -phi_max)) then
		phi(1) = -(phi(1) + phi_max)
	else if ((phi(1) .gt. -phi_max) .and. (phi(1) .lt. -phi_hmax)) then
		phi(1) = phi(1) + phi_max
	else if ((phi(1) .gt. -phi_hmax) .and. (phi(1) .lt. phi_min)) then
		phi(1) = -phi(1)
	end if
	if ((phi(2) .gt. phi_hmax) .and. (phi(2) .lt. phi_max)) then
		phi(2) = phi_max - phi(2)
	else if ((phi(2) .gt. phi_max) .and. (phi(2) .lt. pi)) then
		phi(2) = phi(2) - phi_max
	else if ((phi(2) .gt. -pi) .and. (phi(2) .lt. -phi_max)) then
		phi(2) = -(phi(2) + phi_max)
	else if ((phi(2) .gt. -phi_max) .and. (phi(2) .lt. -phi_hmax)) then
		phi(2) = phi(2) + phi_max
	else if ((phi(2) .gt. -phi_hmax) .and. (phi(2) .lt. phi_min)) then
		phi(2) = -phi(2)
	end if

	tSolv1 = cross_product(rSolv1,h) ! this is (r1 x p) ie. the t vector
	sSolv1 = cross_product(tSolv1,rSolv1)
	tSolv1n = norm2(tSolv1)
	sSolv1n = norm2(sSolv1)

end subroutine calc_angles


! trilinearly interpolate
subroutine trilin_interpolate(ihist, fP)
	use cfgData
	use histData
	use angleData
	implicit none
	integer			:: i, j, iloop, ihist, ir0, ir1, ic0, ic1, ip0, ip1
	real(kind=dp)	:: f_index, r0, r1, c0, c1, p0, p1, rd, cd, pd, f00, f01, f10, f11, f0, f1, fP(3), rSolvnInt(2), cosThInt(2), &
		phiInt(2)

	if (g_fit .eq. 'lin') then
		fP = 1_dp
	else if (g_fit .eq. 'log') then
		fP = 0_dp
	end if
	if (ihist .eq. 1) then ! for g(r) both particles need to be taken into account
		iloop = 2
	else if (ihist .eq. 2) then ! for frc(r) only particle 1 needs to be accounted for
		iloop = 1
	end if
	do i = 1, iloop ! Note: get the bounding values in each dimension.
		! r
		f_index = (rSolvn(i)) / histDistStepSize + 0.5_dp ! take into account half-bin positions
		ir0 = floor(f_index) ! get flanking r indicies
		if (ir0 .ge. histDistBins) then
			ir0 = histDistBins
			ir1 = histDistBins
		else if (ir0 .lt. 1) then
			ir0 = 1
			ir1 = 1
		else
			ir1 = ir0 + 1
		end if
		r0 = histDist(ir0)
		r1 = histDist(ir1)
		if ((rSolvn(i) .lt. r0) .or. (rSolvn(i) .gt. r1)) then
			rSolvnInt(i) = r0 ! when rSolvn is outside the bounds it gets set to r0=r1.
		else
			rSolvnInt(i) = rSolvn(i)
		end if

		! cosTheta
		f_index = (cosTh(i) - cosTh_min) / histCosThStepSize + 0.5_dp ! so index values start at 1
		ic0 = floor(f_index) ! get flanking cosTh indicies
		if (ic0 .ge. histCosThBins) then
			ic0 = histCosThBins
			ic1 = histCosThBins
		else if (ic0 .lt. 1) then
			ic0 = 1
			ic1 = 1
		else
			ic1 = ic0 + 1
		end if
		c0 = histCosTh(ic0)
		c1 = histCosTh(ic1)
		if ((cosTh(i) .lt. c0) .or. (cosTh(i) .gt. c1)) then
			cosThInt(i) = c0
		else
			cosThInt(i) = cosTh(i)
		end if

		! phi
		f_index = (phi(i) - phi_min) / histPhiStepSize + 0.5_dp ! so index values start at 1
		ip0 = floor(f_index) ! get flanking phi indicies
		if (ip0 .ge. histPhiBins) then
			ip0 = histPhiBins
			ip1 = histPhiBins
		else if (ip0 .lt. 1) then
			ip0 = 1
			ip1 = 1
		else
			ip1 = ip0 + 1
		end if
		p0 = histPhi(ip0)
		p1 = histPhi(ip1)
		if ((phi(i) .lt. p0) .or. (phi(i) .gt. p1)) then
			phiInt(i) = p0
		else
			phiInt(i) = phi(i)
		end if

		! Note: set fractional distances
		if ( ir0 .eq. ir1) then ! if r1=r0 then rd would become a NaN.
			rd = 1_dp
		else
			rd = (rSolvnInt(i)-r0)/(r1-r0)
		end if
		if ( ic0 .eq. ic1 ) then ! if c1=c0 then cd would become a NaN.
			cd = 1_dp
		else
			cd = (cosThInt(i)-c0)/(c1-c0)
		end if
		if ( ip0 .eq. ip1 ) then ! if p1=p0 then pd would become a NaN.
			pd = 1_dp
		else
			pd = (phiInt(i)-p0)/(p1-p0)
		end if

		f00 = hist3D(ihist,ir0,ic0,ip0)*(1-rd) + hist3D(ihist,ir1,ic0,ip0)*rd
		f01 = hist3D(ihist,ir0,ic0,ip1)*(1-rd) + hist3D(ihist,ir1,ic0,ip1)*rd
		f10 = hist3D(ihist,ir0,ic1,ip0)*(1-rd) + hist3D(ihist,ir1,ic1,ip0)*rd
		f11 = hist3D(ihist,ir0,ic1,ip1)*(1-rd) + hist3D(ihist,ir1,ic1,ip1)*rd

		f0 = f00*(1-cd) + f10*cd
		f1 = f01*(1-cd) + f11*cd

		if (g_fit .eq. 'lin') then
			fP(1) = fP(1) * (f0*(1-pd) + f1*pd) ! forces are additive and gr is multiplicative
		else if (g_fit .eq. 'log') then
			fP(1) = fP(1) + (f0*(1-pd) + f1*pd) ! forces are additive and FE is additive.
		end if

		if (ihist .eq. 2) then ! it's a force calculation
			do j = 3, 4
				f00 = hist3D(j,ir0,ic0,ip0)*(1-rd) + hist3D(j,ir1,ic0,ip0)*rd
				f01 = hist3D(j,ir0,ic0,ip1)*(1-rd) + hist3D(j,ir1,ic0,ip1)*rd
				f10 = hist3D(j,ir0,ic1,ip0)*(1-rd) + hist3D(j,ir1,ic1,ip0)*rd
				f11 = hist3D(j,ir0,ic1,ip1)*(1-rd) + hist3D(j,ir1,ic1,ip1)*rd

				f0 = f00*(1-cd) + f10*cd
				f1 = f01*(1-cd) + f11*cd

				if (g_fit .eq. 'lin') then
					fP(j-1) = fP(j-1) * (f0*(1-pd) + f1*pd) ! forces are additive and g(r) is mulitiplicative
				else if (g_fit .eq. 'log') then
					fP(j-1) = fP(j-1) + (f0*(1-pd) + f1*pd) ! forces are additive and FE is additive.
				end if
			end do
		end if
	end do

end subroutine trilin_interpolate


! integrate the average force from 'compute_avg_force' to get the PMF.
subroutine integrate_force
	use cfgData
	use ctrlData
	implicit none
	integer 		:: d

	allocate( u_dir(cfgRBins) )
	u_dir = 0_dp

	if (explicit_R .eqv. .false.) then
		do d = 1, cfgRBins
			if (d .eq. 1) then
				u_dir(cfgRBins) = fAvg(cfgRBins) * RStepSize
			else
				u_dir(cfgRBins-(d-1)) = u_dir(cfgRBins-(d-2)) + fAvg(cfgRBins-(d-1)) * RStepSize
			end if
		end do
	else if (explicit_R .eqv. .true.) then
		do d = 1, cfgRBins
			if (d .eq. 1) then
				u_dir(cfgRBins) = fAvg(cfgRBins) * (R_axis(cfgRBins)-R_axis(cfgRBins-1))
				print*, (R_axis(cfgRBins)-R_axis(cfgRBins-1))
			else
				! FIXME: is the delta R part of this correct?
				u_dir(cfgRBins-(d-1)) = u_dir(cfgRBins-(d-2)) + fAvg(cfgRBins-(d-1)) * &
					(R_axis(cfgRBins-(d-1))-R_axis(cfgRBins-d))
				! it looks like the first value is getting printed twice. Also, the values might be wrong. Should it be (d-1) and
				! (d-0)? instead of -2 and -1?
				print*, (R_axis(cfgRBins-(d-1))-R_axis(cfgRBins-d))
			end if
		end do
	end if

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
	write(35,*) "# 4.	Force.r"
	write(35,*) "# 5.	Force.s"
	write(35,*) "# 6.	Force.t"
	do j = 1, zBins
		do i = 1, xBins
			write(35,898) x_axis(i), z_axis(j), grSPA(i,j), frcSPA(1,i,j), frcSPA(2,i,j), frcSPA(3,i,j)
		end do
	end do
	close(35)

	flush(6)
	
898		format (6(1x,es14.7))

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

	flush(6)

899		format (3(1x,es14.7)) ! scientific format

end subroutine write_output
