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
!					<-----R
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!   Modules   !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! data for the density and force tables.
module histData
	use prec
	real(kind=dp),allocatable	:: histDist(:), histCosTh(:), histPhi(:), g(:,:,:), fr(:,:,:), fs(:,:,:), ft(:,:,:), g2(:,:,:), &
		& fr2(:,:,:), fs2(:,:,:), ft2(:,:,:), idealHist(:,:,:,:), gc(:,:,:), gTmp1(:,:,:), gTmp2(:,:,:), frTmp1(:,:,:), &
		& fsTmp1(:,:,:), ftTmp1(:,:,:)
	real(kind=dp)	:: histDistStepSize, histCosThStepSize, histPhiStepSize
	integer	:: histDistBins, histCosThBins, histPhiBins

end module histData

! data from the config file.
module cfgData
	use prec
	use constants
	real(kind=dp),allocatable :: x_axis(:), z_axis(:), R_axis(:), fAvg(:), u_dir(:)
	real(kind=dp) :: RStepSize, xzStepSize, R_min, R_max, xz_range, cfgCosThStepSize, cfgPhiStepSize, cfgPsiStepSize, T
	character(len=8) :: c_explicit_R
	integer :: cfgRBins, cfgCosThBins, cfgPhiBins, cfgPsiBins
!
	real(kind=dp)	:: cosTh_max = 1_dp
	real(kind=dp)	:: cosTh_min = -1_dp
	real(kind=dp)	:: phi_max = 2_dp*pi/3_dp ! This is sufficient for a molecule with C3 symmetry.
	real(kind=dp)	:: phi_hmax = pi/3_dp
	real(kind=dp)	:: phi_min = 0_dp
end module cfgData

! data for calculating cosTh value.
module angleData
	use prec
	real(kind=dp),allocatable :: sinThetaLF(:), cosThetaLF(:), sinPhiLF(:), cosPhiLF(:), sinPsiLF(:), cosPsiLF(:)
	real(kind=dp) :: rSolv1(3), rSolv2(3), rSolvn(2), cosTh(2), phi(2), sSolv1(3), tSolv1(3), sSolv1n, tSolv1n

end module angleData

! testing arrays for force and g(r)
module ctrlData
	use prec
	real(kind=dp),allocatable	:: frcSPA(:,:,:), grSPA(:,:), explicitDist(:)
	integer :: crdLines
	logical :: explicit_R

end module ctrlData


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!  Main Program  !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program compute_avgForce
	use prec
	implicit none
	character(len=64) :: histFile, cfgFile, outFile
	real(kind=dp) 		:: omp_get_wtime, ti, tf, seconds
	integer				:: hours, minutes

	ti = omp_get_wtime()

	! make list of average direct force from 'collapsed' file.
	call parse_command_line(histFile, cfgFile) !, outFile)

	! read config file
	call read_cfg(cfgFile, outFile)

	! make list of average direct force from 'collapsed' file.
	call make_hist_table(histFile)
	
	! Now that we have the relevant information spline the g and f arrays along r.
	call spline_hist_array

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

	hours = (tf-ti)/3600
	minutes = mod((tf-ti),3600d0)/60
	seconds = mod(mod((tf-ti),3600d0),60d0)

	write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	write(*,'(a,i4,a,i2,a,f6.3,a)') "Total time elapsed:	", hours, "h ", minutes, "m ", seconds, "s"

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
			write(*,*) "Histogram File:			", histFile
		case ('-cfg')
			i = i+1
			call get_command_argument(i,cfgFile)
			cfgFileFlag=.true.
			write(*,*) "Config File:				", cfgFile
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
	logical 			:: outFileFlag, RstepSizeFlag, xzStepSizeFlag, RmaxFlag, RminFlag, xzRangeFlag, thetaBinsFlag, phiBinsFlag, &
		& psiBinsFlag, c_explicit_RFlag, TFlag

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
	TFlag = .false.

	ios = 0

	open(20,file=cfgFile)
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		call split(line,'=',firstWord, sep)
		if (line .ne. "") then
			if (firstWord .eq. "out_file") then
				read(line,*) outFile
				write(*,*) "Output File:				", outFile
				outFileFlag = .true.
			else if (firstWord .eq. "RStepSize") then
				read(line,*) RStepSize
				write(*,*) "PMF Step Size:				", RStepSize
				RstepSizeFlag = .true.
			else if (firstWord .eq. "xzStepSize") then
				read(line,*) xzStepSize
				write(*,*) "Solvent Grid Step Size:		", xzStepSize
				xzStepSizeFlag = .true.
			else if (firstWord .eq. "explicit_R") then
				read(line,*) c_explicit_R
				write(*,*) "Use Explicit R Values:			", c_explicit_R
				c_explicit_RFlag = .true.
			else if (firstWord .eq. "R_max") then
				read(line,*) R_max
				write(*,*) "R Maximum Value:			", R_max
				RmaxFlag = .true.
			else if (firstWord .eq. "R_min") then
				read(line,*) R_min
				write(*,*) "R Minimum Value:			", R_min
				RminFlag = .true.
			else if (firstWord .eq. "xz_range") then
				read(line,*) xz_range
				write(*,*) "XZ - Range:				", xz_range
				xzRangeFlag = .true.
			else if (firstWord .eq. "theta_bins") then
				read(line,*) cfgCosThBins
				write(*,*) "Theta Bins:			", cfgCosThBins
				thetaBinsFlag= .true.
			else if (firstWord .eq. "phi_bins") then
				read(line,*) cfgPhiBins
				write(*,*) "Phi Bins:			", cfgPhiBins
				phiBinsFlag= .true.
			else if (firstWord .eq. "psi_bins") then
				read(line,*) cfgPsiBins
				write(*,*) "Psi Bins:			", cfgPsiBins
				psiBinsFlag= .true.
			else if (firstWord .eq. "temperature") then
				read(line,*) T
				write(*,*) "Temperature (K):			", T
				TFlag= .true.
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
	if (TFlag.eqv..false.) then
		write(*,*) "Config file must have a 'temperature' value"
		stop
	end if

	flush(6)

end subroutine read_cfg


! read force file and make a lookup table.
subroutine make_hist_table(histFile)
	use histData
	use cfgData
	implicit none
	character(len=64)				:: histFile
	character(len=32)				:: junk
	character(len=256)			:: line
	integer							:: ios, ios2, i, j, k, nHistLines
	real(kind=dp),allocatable	:: histTmp(:,:)

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
	!write(*,*) "nHistLines", nHistLines

	allocate( histTmp(8,nHistLines) )

	! populate hist arrays
	ios = 0; i = 1
	open(20,file=histFile)
	! read file ignoring comment lines at the beginning
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		if ((line(1:1) .ne. "#") .and. (ios .ge. 0)) then
			!						dist		   cos(Th)			phi/3			g(r)+			g(r)-		<f.r>+			<f.s>+		<f.t>+
			read(line,*) histTmp(1,i), histTmp(2,i), histTmp(3,i), histTmp(4,i), junk, histTmp(5,i), histTmp(6,i), histTmp(7,i), &
				!		gc(r)+	 gc(r)-
				& histTmp(8,i), junk
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
	write(*,*) "Histogram Distance Bins:	", histDistBins
	write(*,*) "Histogram Cosine Theta Bins:	", histCosThBins
	write(*,*) "Histogram Phi Bins:		", histPhiBins
	
	allocate( histDist(histDistBins), histCosTh(histCosThBins), histPhi(histPhiBins), g(histDistBins,histCosThBins,histPhiBins), &
		& fr(histDistBins,histCosThBins,histPhiBins), fs(histDistBins,histCosThBins,histPhiBins), &
		& ft(histDistBins,histCosThBins,histPhiBins), gc(histDistBins,histCosThBins,histPhiBins) )
	
	! Allocate temporary wrapped angular arrays for bicubic interpolation
	allocate( gTmp1(4,0:histCosThBins+1,0:histPhiBins+1), gTmp2(4,0:histCosThBins+1,0:histPhiBins+1), &
		& frTmp1(4,0:histCosThBins+1,0:histPhiBins+1), fsTmp1(4,0:histCosThBins+1,0:histPhiBins+1), &
		& ftTmp1(4,0:histCosThBins+1,0:histPhiBins+1) )
	
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
				g(i,j,k) = histTmp(4, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k) ! g(r,cos,phi)
				fr(i,j,k) = histTmp(5, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k) ! <f.r>(r,cos,phi)
				fs(i,j,k) = histTmp(6, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k) ! <f.s>(r,cos,phi)
				ft(i,j,k) = histTmp(7, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k) ! <f.t>(r,cos,phi)
				gc(i,j,k) = histTmp(8, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k) ! gc(r,cos,phi)
			end do
		end do
	end do

	histDistStepSize = histDist(2) - histDist(1)
	write(*,*) "Histogram Distance Step Size:		", histDistStepSize
	histCosThStepSize = histCosTh(2) - histCosTh(1)
	write(*,*) "Histogram Cosine Theta Step Size:	", histCosThStepSize
	histPhiStepSize = histPhi(2) - histPhi(1)
	write(*,*) "Histogram Phi Step Size:		", histPhiStepSize 

	flush(6)

end subroutine make_hist_table


! spline the r dimension of each theta phi stack.
subroutine spline_hist_array
	use constants
	use functions
	use histData
	use cfgData
	use idealSolv
	implicit none
	integer				:: i, ir, ith, iphi, imin
	real(kind=dp)	:: x, y, norm_factor

	! Calculate a 4D array idealHist(g/f,r,th,phi)
	allocate( idealHist(4,histDistBins,histCosThBins,histPhiBins) )
	idealHist = 0_dp
	call ideal_CL3(histDistBins,histDistStepSize,histCosThBins,cosTh_min,cosTh_max,histPhiBins,phi_min,phi_max, idealHist)

	! Edit the input hist arrays to more smoothly transition to +/- infinity with the help of idealHist.
	do ith = 1, histCosThBins
		do iphi = 1, histPhiBins
			imin = 0
			! Normalization factor or each theta phi array.
			norm_factor = gc(histDistBins,ith,iphi)/(g(histDistBins,ith,iphi)*4*pi*histDist(histDistBins)**2)
find:		do ir = 1, histDistBins
				if (g(ir,ith,iphi).gt.1d-6) then
					imin = ir
					exit find
				end if
			end do find

			! NOTE: Add the ideal values to bins with no sampling.
			do ir = histDistBins, 1, -1
				if (ir.ge.imin) then
					if (g(ir,ith,iphi).gt.1d-6) then
						g(ir,ith,iphi) = log(g(ir,ith,iphi))
					else	! note: This is a zero bin where there probably should have been something. So put a single count in.
						g(ir,ith,iphi) = log(real(1,dp)/(4*pi*histDist(ir)**2)/norm_factor)
						fr(ir,ith,iphi) = idealHist(2,ir,ith,iphi)
						fs(ir,ith,iphi) = idealHist(3,ir,ith,iphi)
						ft(ir,ith,iphi) = idealHist(4,ir,ith,iphi)
					end if
				else if (ir.lt.imin) then ! .lt.imin ==> in the region of no sampling.
					! log(g(r<r0)) = -u_pmf(r)/T - ( u_pmf(r0)/T + u_dir(r0)/T )
					g(ir,ith,iphi) = ( -idealHist(1,ir,ith,iphi) - ( g(imin,ith,iphi) - idealHist(1,imin,ith,iphi) ) )/kb_kcalmolK/T
					if (g(ir,ith,iphi).lt.-1d4) then
						g(ir,ith,iphi) = -1d4
					end if
					fr(ir,ith,iphi) = idealHist(2,ir,ith,iphi)
					if (fr(ir,ith,iphi).gt.1d4) then
						fr(ir,ith,iphi) = 1d4
					end if
					fs(ir,ith,iphi) = idealHist(3,ir,ith,iphi)
					if (fs(ir,ith,iphi).gt.1d4) then
						fs(ir,ith,iphi) = 1d4
					end if
					ft(ir,ith,iphi) = idealHist(4,ir,ith,iphi)
					if (ft(ir,ith,iphi).gt.1d4) then
						ft(ir,ith,iphi) = 1d4
					end if
				end if
			end do
		end do
	end do

	allocate( g2(histDistBins,histCosThBins,histPhiBins), fr2(histDistBins,histCosThBins,histPhiBins), & 
		& fs2(histDistBins,histCosThBins,histPhiBins), ft2(histDistBins,histCosThBins,histPhiBins) )

	! todo: edited g and g2
	open(45)
	write(45,*) '# just g and its second derivatives'
	do ith = 1, 1!histcosthbins
		do iphi = 3, 3!histphibins
			call spline(histDist,g(:,ith,iphi),histDistBins,dble(0),dble(0),g2(:,ith,iphi))
			do ir = 1, histdistbins
				write(45,*) histDist(ir), histcosth(ith), histphi(iphi), gc(ir,ith,iphi)/4/pi/histDist(ir)**2/norm_factor, &
					& g(ir,ith,iphi), g2(ir,ith,iphi), gc(ir,ith,iphi)
			end do
		end do
	end do
	close(45)

	! todo: test interpolation 
	open(55)
	write(55,*) '# test interpolation of g'
	do ith = 1, 1!histcosthbins
		do iphi = 3, 3!histphibins
			do i = 1, 2500
				x = i*histDistStepSize/10_dp
				call splint(histDist,g(:,ith,iphi),g2(:,ith,iphi),histDistBins,x,y)
				!print*, x, y
				write(55,*) x, histcosth(ith), histphi(iphi), y
			end do
		end do
	end do
	close(55)

	! todo: write out ideal hist to make sure it's okay.
	open(65)
	write(65,*) '# test idealhist output'
	do ith = 1, histcosthbins
		do iphi = 1, histphibins
			do ir = 1, histdistbins
				write(65,*) histdist(ir), histcosth(ith), histphi(iphi), idealHist(1,ir,ith,iphi), idealHist(2,ir,ith,iphi), &
					& idealHist(3,ir,ith,iphi), idealHist(4,ir,ith,iphi)
			end do
		end do
	end do
	close(65)

end subroutine spline_hist_array


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


! Populate and wrap the temporary (cosTh,phi) arrays into a torus.
subroutine set_tmp_arrays
	use histData
	use angleData
	use functions
	implicit none
	integer			:: j, k
	real(kind=dp)	:: x1a(0:histCosThBins+1), x2a(0:histPhiBins+1)

	! todo: Evaluate the splines for every theta,phi point at these distances and save those arrays at gTmp1, gTmp2, and fTmp1. The
	! numbers refer to which solute they correpsond to. These arrays are saved and used until the next distance.
	do j = 1, histCosThBins
		do k = 1, histPhiBins
			if ((rSolvn(1).gt.histDist(histDistBins)).and.(rSolvn(2).gt.histDist(histDistBins))) then
				call splint(histDist,g(:,j,k),g2(:,j,k),histDistBins,histDist(histDistBins-1),gTmp1(1,j,k))
				call splint(histDist,g(:,j,k),g2(:,j,k),histDistBins,histDist(histDistBins-1),gTmp2(1,j,k))
				call splint(histDist,fr(:,j,k),fr2(:,j,k),histDistBins,histDist(histDistBins-1),frTmp1(1,j,k))
				call splint(histDist,fs(:,j,k),fs2(:,j,k),histDistBins,histDist(histDistBins-1),fsTmp1(1,j,k))
				call splint(histDist,ft(:,j,k),ft2(:,j,k),histDistBins,histDist(histDistBins-1),ftTmp1(1,j,k))
			else if (rSolvn(1).gt.histDist(histDistBins)) then
				call splint(histDist,g(:,j,k),g2(:,j,k),histDistBins,histDist(histDistBins-1),gTmp1(1,j,k))
				call splint(histDist,g(:,j,k),g2(:,j,k),histDistBins,rSolvn(2),gTmp2(1,j,k))
				call splint(histDist,fr(:,j,k),fr2(:,j,k),histDistBins,histDist(histDistBins-1),frTmp1(1,j,k))
				call splint(histDist,fs(:,j,k),fs2(:,j,k),histDistBins,histDist(histDistBins-1),fsTmp1(1,j,k))
				call splint(histDist,ft(:,j,k),ft2(:,j,k),histDistBins,histDist(histDistBins-1),ftTmp1(1,j,k))
			else if (rSolvn(2).gt.histDist(histDistBins)) then
				call splint(histDist,g(:,j,k),g2(:,j,k),histDistBins,rSolvn(1),gTmp1(1,j,k))
				call splint(histDist,g(:,j,k),g2(:,j,k),histDistBins,histDist(histDistBins-1),gTmp2(1,j,k))
				call splint(histDist,fr(:,j,k),fr2(:,j,k),histDistBins,rSolvn(1),frTmp1(1,j,k))
				call splint(histDist,fs(:,j,k),fs2(:,j,k),histDistBins,rSolvn(1),fsTmp1(1,j,k))
				call splint(histDist,ft(:,j,k),ft2(:,j,k),histDistBins,rSolvn(1),ftTmp1(1,j,k))
			else
				call splint(histDist,g(:,j,k),g2(:,j,k),histDistBins,rSolvn(1),gTmp1(1,j,k))
				call splint(histDist,g(:,j,k),g2(:,j,k),histDistBins,rSolvn(2),gTmp2(1,j,k))
				call splint(histDist,fr(:,j,k),fr2(:,j,k),histDistBins,rSolvn(1),frTmp1(1,j,k))
				call splint(histDist,fs(:,j,k),fs2(:,j,k),histDistBins,rSolvn(1),fsTmp1(1,j,k))
				call splint(histDist,ft(:,j,k),ft2(:,j,k),histDistBins,rSolvn(1),ftTmp1(1,j,k))
			end if
		end do
	end do
	! todo: calculate derivatives and store them in the tmp arrays
	! make independent variable arrays that have the extra 2 elements for taking derivatives.
	x1a(1:histCosThBins) = histCosTh(:)
	x1a(0) = histCosTh(1)
	x1a(histCosThBins+1) = histCosTh(histCosThBins)
	x2a(1:histPhiBins) = histPhi(:)
	x2a(0) = histPhi(1)
	x2a(histPhiBins+1) = histPhi(histPhiBins)
	! Centered difference method to find derivatives.
	do j = 1, histCosThBins
		do k = 1, histPhiBins
			gTmp1(2,j,k) = (gTmp1(1,j+1,k)-gTmp1(1,j-1,k))/(x1a(j+1)-x1a(j-1))
			gTmp1(3,j,k) = (gTmp1(1,j,k+1)-gTmp1(1,j,k-1))/(x2a(k+1)-x2a(k-1))
			gTmp1(4,j,k) = (gTmp1(1,j+1,k+1)-gTmp1(1,j+1,k-1)-gTmp1(1,j-1,k+1)+gTmp1(1,j-1,k-1))/((x1a(j+1)-x1a(j-1))*(x2a(k+1)-x2a(k-1)))
			gTmp2(2,j,k) = (gTmp2(1,j+1,k)-gTmp2(1,j-1,k))/(x1a(j+1)-x1a(j-1))
			gTmp2(3,j,k) = (gTmp2(1,j,k+1)-gTmp2(1,j,k-1))/(x2a(k+1)-x2a(k-1))
			gTmp2(4,j,k) = (gTmp2(1,j+1,k+1)-gTmp2(1,j+1,k-1)-gTmp2(1,j-1,k+1)+gTmp2(1,j-1,k-1))/((x1a(j+1)-x1a(j-1))*(x2a(k+1)-x2a(k-1)))
			frTmp1(2,j,k) = (frTmp1(1,j+1,k)-frTmp1(1,j-1,k))/(x1a(j+1)-x1a(j-1))
			frTmp1(3,j,k) = (frTmp1(1,j,k+1)-frTmp1(1,j,k-1))/(x2a(k+1)-x2a(k-1))
			frTmp1(4,j,k) = (frTmp1(1,j+1,k+1)-frTmp1(1,j+1,k-1)-frTmp1(1,j-1,k+1)+frTmp1(1,j-1,k-1))/((x1a(j+1)-x1a(j-1))*(x2a(k+1)-x2a(k-1)))
			fsTmp1(2,j,k) = (fsTmp1(1,j+1,k)-fsTmp1(1,j-1,k))/(x1a(j+1)-x1a(j-1))
			fsTmp1(3,j,k) = (fsTmp1(1,j,k+1)-fsTmp1(1,j,k-1))/(x2a(k+1)-x2a(k-1))
			fsTmp1(4,j,k) = (fsTmp1(1,j+1,k+1)-fsTmp1(1,j+1,k-1)-fsTmp1(1,j-1,k+1)+fsTmp1(1,j-1,k-1))/((x1a(j+1)-x1a(j-1))*(x2a(k+1)-x2a(k-1)))
			ftTmp1(2,j,k) = (ftTmp1(1,j+1,k)-ftTmp1(1,j-1,k))/(x1a(j+1)-x1a(j-1))
			ftTmp1(3,j,k) = (ftTmp1(1,j,k+1)-ftTmp1(1,j,k-1))/(x2a(k+1)-x2a(k-1))
			ftTmp1(4,j,k) = (ftTmp1(1,j+1,k+1)-ftTmp1(1,j+1,k-1)-ftTmp1(1,j-1,k+1)+ftTmp1(1,j-1,k-1))/((x1a(j+1)-x1a(j-1))*(x2a(k+1)-x2a(k-1)))
		end do
	end do

	! xxx: Wrap values:
	! Wrap all 4 corners
	gTmp1(1,0,0) = gTmp1(1,1,1)
	gTmp1(1,histCosThBins+1,0) = gTmp1(1,histCosThBins,1)
	gTmp1(1,0,histPhiBins+1) = gTmp1(1,1,histPhiBins)
	gTmp1(1,histCosThBins+1,histPhiBins+1) = gTmp1(1,histCosThBins,histPhiBins)
	gTmp2(1,0,0) = gTmp2(1,1,1)
	gTmp2(1,histCosThBins+1,0) = gTmp2(1,histCosThBins,1)
	gTmp2(1,0,histPhiBins+1) = gTmp2(1,1,histPhiBins)
	gTmp2(1,histCosThBins+1,histPhiBins+1) = gTmp2(1,histCosThBins,histPhiBins)
	frTmp1(1,0,0) = frTmp1(1,1,1)
	frTmp1(1,histCosThBins+1,0) = frTmp1(1,histCosThBins,1)
	frTmp1(1,0,histPhiBins+1) = frTmp1(1,1,histPhiBins)
	frTmp1(1,histCosThBins+1,histPhiBins+1) = frTmp1(1,histCosThBins,histPhiBins)
	fsTmp1(1,0,0) = fsTmp1(1,1,1)
	fsTmp1(1,histCosThBins+1,0) = fsTmp1(1,histCosThBins,1)
	fsTmp1(1,0,histPhiBins+1) = fsTmp1(1,1,histPhiBins)
	fsTmp1(1,histCosThBins+1,histPhiBins+1) = fsTmp1(1,histCosThBins,histPhiBins)
	ftTmp1(1,0,0) = ftTmp1(1,1,1)
	ftTmp1(1,histCosThBins+1,0) = ftTmp1(1,histCosThBins,1)
	ftTmp1(1,0,histPhiBins+1) = ftTmp1(1,1,histPhiBins)
	ftTmp1(1,histCosThBins+1,histPhiBins+1) = ftTmp1(1,histCosThBins,histPhiBins)
	! Wrap all 4 edges
	gTmp1(1,0 , 1:histPhiBins) = gTmp1(1,1 , 1:histPhiBins)
	gTmp1(1,histCosThBins+1 , 1:histPhiBins) = gTmp1(1,histCosThBins , 1:histPhiBins)
	gTmp1(1,1:histCosThBins , 0) = gTmp1(1,1:histCosThBins , 1)
	gTmp1(1,1:histCosThBins , histPhiBins+1) = gTmp1(1,1:histCosThBins , histPhiBins)
	gTmp2(1,0 , 1:histPhiBins) = gTmp2(1,1 , 1:histPhiBins)
	gTmp2(1,histCosThBins+1 , 1:histPhiBins) = gTmp2(1,histCosThBins , 1:histPhiBins)
	gTmp2(1,1:histCosThBins , 0) = gTmp2(1,1:histCosThBins , 1)
	gTmp2(1,1:histCosThBins , histPhiBins+1) = gTmp2(1,1:histCosThBins , histPhiBins)
	frTmp1(1,0 , 1:histPhiBins) = frTmp1(1,1 , 1:histPhiBins)
	frTmp1(1,histCosThBins+1 , 1:histPhiBins) = frTmp1(1,histCosThBins , 1:histPhiBins)
	frTmp1(1,1:histCosThBins , 0) = frTmp1(1,1:histCosThBins , 1)
	frTmp1(1,1:histCosThBins , histPhiBins+1) = frTmp1(1,1:histCosThBins , histPhiBins)
	fsTmp1(1,0 , 1:histPhiBins) = fsTmp1(1,1 , 1:histPhiBins)
	fsTmp1(1,histCosThBins+1 , 1:histPhiBins) = fsTmp1(1,histCosThBins , 1:histPhiBins)
	fsTmp1(1,1:histCosThBins , 0) = fsTmp1(1,1:histCosThBins , 1)
	fsTmp1(1,1:histCosThBins , histPhiBins+1) = fsTmp1(1,1:histCosThBins , histPhiBins)
	ftTmp1(1,0 , 1:histPhiBins) = ftTmp1(1,1 , 1:histPhiBins)
	ftTmp1(1,histCosThBins+1 , 1:histPhiBins) = ftTmp1(1,histCosThBins , 1:histPhiBins)
	ftTmp1(1,1:histCosThBins , 0) = ftTmp1(1,1:histCosThBins , 1)
	ftTmp1(1,1:histCosThBins , histPhiBins+1) = ftTmp1(1,1:histCosThBins , histPhiBins)

	! xxx: Wrap x1 derivatives
	! Wrap all 4 corners
	gTmp1(2,0,0) = -gTmp1(2,1,1)
	gTmp1(2,histCosThBins+1,0) = -gTmp1(2,histCosThBins,1)
	gTmp1(2,0,histPhiBins+1) = -gTmp1(2,1,histPhiBins)
	gTmp1(2,histCosThBins+1,histPhiBins+1) = -gTmp1(2,histCosThBins,histPhiBins)
	gTmp2(2,0,0) = -gTmp2(2,1,1)
	gTmp2(2,histCosThBins+1,0) = -gTmp2(2,histCosThBins,1)
	gTmp2(2,0,histPhiBins+1) = -gTmp2(2,1,histPhiBins)
	gTmp2(2,histCosThBins+1,histPhiBins+1) = -gTmp2(2,histCosThBins,histPhiBins)
	frTmp1(2,0,0) = -frTmp1(2,1,1)
	frTmp1(2,histCosThBins+1,0) = -frTmp1(2,histCosThBins,1)
	frTmp1(2,0,histPhiBins+1) = -frTmp1(2,1,histPhiBins)
	frTmp1(2,histCosThBins+1,histPhiBins+1) = -frTmp1(2,histCosThBins,histPhiBins)
	fsTmp1(2,0,0) = -fsTmp1(2,1,1)
	fsTmp1(2,histCosThBins+1,0) = -fsTmp1(2,histCosThBins,1)
	fsTmp1(2,0,histPhiBins+1) = -fsTmp1(2,1,histPhiBins)
	fsTmp1(2,histCosThBins+1,histPhiBins+1) = -fsTmp1(2,histCosThBins,histPhiBins)
	ftTmp1(2,0,0) = -ftTmp1(2,1,1)
	ftTmp1(2,histCosThBins+1,0) = -ftTmp1(2,histCosThBins,1)
	ftTmp1(2,0,histPhiBins+1) = -ftTmp1(2,1,histPhiBins)
	ftTmp1(2,histCosThBins+1,histPhiBins+1) = -ftTmp1(2,histCosThBins,histPhiBins)
	! Wrap all 4 edges
	gTmp1(2,0 , 1:histPhiBins) = -gTmp1(2,1 , 1:histPhiBins)
	gTmp1(2,histCosThBins+1 , 1:histPhiBins) = -gTmp1(2,histCosThBins , 1:histPhiBins)
	gTmp1(2,1:histCosThBins , 0) = gTmp1(2,1:histCosThBins , 1)
	gTmp1(2,1:histCosThBins , histPhiBins+1) = gTmp1(2,1:histCosThBins , histPhiBins)
	gTmp2(2,0 , 1:histPhiBins) = -gTmp2(2,1 , 1:histPhiBins)
	gTmp2(2,histCosThBins+1 , 1:histPhiBins) = -gTmp2(2,histCosThBins , 1:histPhiBins)
	gTmp2(2,1:histCosThBins , 0) = gTmp2(2,1:histCosThBins , 1)
	gTmp2(2,1:histCosThBins , histPhiBins+1) = gTmp2(2,1:histCosThBins , histPhiBins)
	frTmp1(2,0 , 1:histPhiBins) = -frTmp1(2,1 , 1:histPhiBins)
	frTmp1(2,histCosThBins+1 , 1:histPhiBins) = -frTmp1(2,histCosThBins , 1:histPhiBins)
	frTmp1(2,1:histCosThBins , 0) = frTmp1(2,1:histCosThBins , 1)
	frTmp1(2,1:histCosThBins , histPhiBins+1) = frTmp1(2,1:histCosThBins , histPhiBins)
	fsTmp1(2,0 , 1:histPhiBins) = -fsTmp1(2,1 , 1:histPhiBins)
	fsTmp1(2,histCosThBins+1 , 1:histPhiBins) = -fsTmp1(2,histCosThBins , 1:histPhiBins)
	fsTmp1(2,1:histCosThBins , 0) = fsTmp1(2,1:histCosThBins , 1)
	fsTmp1(2,1:histCosThBins , histPhiBins+1) = fsTmp1(2,1:histCosThBins , histPhiBins)
	ftTmp1(2,0 , 1:histPhiBins) = -ftTmp1(2,1 , 1:histPhiBins)
	ftTmp1(2,histCosThBins+1 , 1:histPhiBins) = -ftTmp1(2,histCosThBins , 1:histPhiBins)
	ftTmp1(2,1:histCosThBins , 0) = ftTmp1(2,1:histCosThBins , 1)
	ftTmp1(2,1:histCosThBins , histPhiBins+1) = ftTmp1(2,1:histCosThBins , histPhiBins)

	! xxx: Wrap x2 derivatives
	! Wrap all 4 corners
	gTmp1(3,0,0) = -gTmp1(3,1,1)
	gTmp1(3,histCosThBins+1,0) = -gTmp1(3,histCosThBins,1)
	gTmp1(3,0,histPhiBins+1) = -gTmp1(3,1,histPhiBins)
	gTmp1(3,histCosThBins+1,histPhiBins+1) = -gTmp1(3,histCosThBins,histPhiBins)
	gTmp2(3,0,0) = -gTmp2(3,1,1)
	gTmp2(3,histCosThBins+1,0) = -gTmp2(3,histCosThBins,1)
	gTmp2(3,0,histPhiBins+1) = -gTmp2(3,1,histPhiBins)
	gTmp2(3,histCosThBins+1,histPhiBins+1) = -gTmp2(3,histCosThBins,histPhiBins)
	frTmp1(3,0,0) = -frTmp1(3,1,1)
	frTmp1(3,histCosThBins+1,0) = -frTmp1(3,histCosThBins,1)
	frTmp1(3,0,histPhiBins+1) = -frTmp1(3,1,histPhiBins)
	frTmp1(3,histCosThBins+1,histPhiBins+1) = -frTmp1(3,histCosThBins,histPhiBins)
	fsTmp1(3,0,0) = -fsTmp1(3,1,1)
	fsTmp1(3,histCosThBins+1,0) = -fsTmp1(3,histCosThBins,1)
	fsTmp1(3,0,histPhiBins+1) = -fsTmp1(3,1,histPhiBins)
	fsTmp1(3,histCosThBins+1,histPhiBins+1) = -fsTmp1(3,histCosThBins,histPhiBins)
	ftTmp1(3,0,0) = -ftTmp1(3,1,1)
	ftTmp1(3,histCosThBins+1,0) = -ftTmp1(3,histCosThBins,1)
	ftTmp1(3,0,histPhiBins+1) = -ftTmp1(3,1,histPhiBins)
	ftTmp1(3,histCosThBins+1,histPhiBins+1) = -ftTmp1(3,histCosThBins,histPhiBins)
	! Wrap all 4 edges
	gTmp1(3,0 , 1:histPhiBins) = gTmp1(3,1 , 1:histPhiBins)
	gTmp1(3,histCosThBins+1 , 1:histPhiBins) = gTmp1(3,histCosThBins , 1:histPhiBins)
	gTmp1(3,1:histCosThBins , 0) = -gTmp1(3,1:histCosThBins , 1)
	gTmp1(3,1:histCosThBins , histPhiBins+1) = -gTmp1(3,1:histCosThBins , histPhiBins)
	gTmp2(3,0 , 1:histPhiBins) = gTmp2(3,1 , 1:histPhiBins)
	gTmp2(3,histCosThBins+1 , 1:histPhiBins) = gTmp2(3,histCosThBins , 1:histPhiBins)
	gTmp2(3,1:histCosThBins , 0) = -gTmp2(3,1:histCosThBins , 1)
	gTmp2(3,1:histCosThBins , histPhiBins+1) = -gTmp2(3,1:histCosThBins , histPhiBins)
	frTmp1(3,0 , 1:histPhiBins) = frTmp1(3,1 , 1:histPhiBins)
	frTmp1(3,histCosThBins+1 , 1:histPhiBins) = frTmp1(3,histCosThBins , 1:histPhiBins)
	frTmp1(3,1:histCosThBins , 0) = -frTmp1(3,1:histCosThBins , 1)
	frTmp1(3,1:histCosThBins , histPhiBins+1) = -frTmp1(3,1:histCosThBins , histPhiBins)
	fsTmp1(3,0 , 1:histPhiBins) = fsTmp1(3,1 , 1:histPhiBins)
	fsTmp1(3,histCosThBins+1 , 1:histPhiBins) = fsTmp1(3,histCosThBins , 1:histPhiBins)
	fsTmp1(3,1:histCosThBins , 0) = -fsTmp1(3,1:histCosThBins , 1)
	fsTmp1(3,1:histCosThBins , histPhiBins+1) = -fsTmp1(3,1:histCosThBins , histPhiBins)
	ftTmp1(3,0 , 1:histPhiBins) = ftTmp1(3,1 , 1:histPhiBins)
	ftTmp1(3,histCosThBins+1 , 1:histPhiBins) = ftTmp1(3,histCosThBins , 1:histPhiBins)
	ftTmp1(3,1:histCosThBins , 0) = -ftTmp1(3,1:histCosThBins , 1)
	ftTmp1(3,1:histCosThBins , histPhiBins+1) = -ftTmp1(3,1:histCosThBins , histPhiBins)

	! xxx: Wrap cross derivatives
	! Wrap all 4 corners
	gTmp1(4,0,0) = gTmp1(4,1,1)
	gTmp1(4,histCosThBins+1,0) = gTmp1(4,histCosThBins,1)
	gTmp1(4,0,histPhiBins+1) = gTmp1(4,1,histPhiBins)
	gTmp1(4,histCosThBins+1,histPhiBins+1) = gTmp1(4,histCosThBins,histPhiBins)
	gTmp2(4,0,0) = gTmp2(4,1,1)
	gTmp2(4,histCosThBins+1,0) = gTmp2(4,histCosThBins,1)
	gTmp2(4,0,histPhiBins+1) = gTmp2(4,1,histPhiBins)
	gTmp2(4,histCosThBins+1,histPhiBins+1) = gTmp2(4,histCosThBins,histPhiBins)
	frTmp1(4,0,0) = frTmp1(4,1,1)
	frTmp1(4,histCosThBins+1,0) = frTmp1(4,histCosThBins,1)
	frTmp1(4,0,histPhiBins+1) = frTmp1(4,1,histPhiBins)
	frTmp1(4,histCosThBins+1,histPhiBins+1) = frTmp1(4,histCosThBins,histPhiBins)
	fsTmp1(4,0,0) = fsTmp1(4,1,1)
	fsTmp1(4,histCosThBins+1,0) = fsTmp1(4,histCosThBins,1)
	fsTmp1(4,0,histPhiBins+1) = fsTmp1(4,1,histPhiBins)
	fsTmp1(4,histCosThBins+1,histPhiBins+1) = fsTmp1(4,histCosThBins,histPhiBins)
	ftTmp1(4,0,0) = ftTmp1(4,1,1)
	ftTmp1(4,histCosThBins+1,0) = ftTmp1(4,histCosThBins,1)
	ftTmp1(4,0,histPhiBins+1) = ftTmp1(4,1,histPhiBins)
	ftTmp1(4,histCosThBins+1,histPhiBins+1) = ftTmp1(4,histCosThBins,histPhiBins)
	! Wrap all 4 edges
	gTmp1(4,0 , 1:histPhiBins) = -gTmp1(4,1 , 1:histPhiBins)
	gTmp1(4,histCosThBins+1 , 1:histPhiBins) = -gTmp1(4,histCosThBins , 1:histPhiBins)
	gTmp1(4,1:histCosThBins , 0) = -gTmp1(4,1:histCosThBins , 1)
	gTmp1(4,1:histCosThBins , histPhiBins+1) = -gTmp1(4,1:histCosThBins , histPhiBins)
	gTmp2(4,0 , 1:histPhiBins) = -gTmp2(4,1 , 1:histPhiBins)
	gTmp2(4,histCosThBins+1 , 1:histPhiBins) = -gTmp2(4,histCosThBins , 1:histPhiBins)
	gTmp2(4,1:histCosThBins , 0) = -gTmp2(4,1:histCosThBins , 1)
	gTmp2(4,1:histCosThBins , histPhiBins+1) = -gTmp2(4,1:histCosThBins , histPhiBins)
	frTmp1(4,0 , 1:histPhiBins) = -frTmp1(4,1 , 1:histPhiBins)
	frTmp1(4,histCosThBins+1 , 1:histPhiBins) = -frTmp1(4,histCosThBins , 1:histPhiBins)
	frTmp1(4,1:histCosThBins , 0) = -frTmp1(4,1:histCosThBins , 1)
	frTmp1(4,1:histCosThBins , histPhiBins+1) = -frTmp1(4,1:histCosThBins , histPhiBins)
	fsTmp1(4,0 , 1:histPhiBins) = -fsTmp1(4,1 , 1:histPhiBins)
	fsTmp1(4,histCosThBins+1 , 1:histPhiBins) = -fsTmp1(4,histCosThBins , 1:histPhiBins)
	fsTmp1(4,1:histCosThBins , 0) = -fsTmp1(4,1:histCosThBins , 1)
	fsTmp1(4,1:histCosThBins , histPhiBins+1) = -fsTmp1(4,1:histCosThBins , histPhiBins)
	ftTmp1(4,0 , 1:histPhiBins) = -ftTmp1(4,1 , 1:histPhiBins)
	ftTmp1(4,histCosThBins+1 , 1:histPhiBins) = -ftTmp1(4,histCosThBins , 1:histPhiBins)
	ftTmp1(4,1:histCosThBins , 0) = -ftTmp1(4,1:histCosThBins , 1)
	ftTmp1(4,1:histCosThBins , histPhiBins+1) = -ftTmp1(4,1:histCosThBins , histPhiBins)


end subroutine set_tmp_arrays


! do the average force integral
subroutine compute_avg_force
	use cfgData
	use histData
	use angleData
	use ctrlData
	use constants
	use functions
	implicit none
	integer 		:: xBins, zBins, r, i, j, ithLF, iphiLF, ipsiLF
	real(kind=dp)	:: density, phiLF, psiLF, psi_max, psi_min, gx, gx2, fx(3), fNew1, fNew2, fNew3

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

	cfgCosThStepSize = (cosTh_max - cosTh_min) / real(cfgCosThBins, dp)
	do ithLF = 1, cfgCosThBins
		cosThetaLF(ithLF) = (ithLF-0.5_dp)*cfgCosThStepSize - cosTh_max
		sinThetaLF(ithLF) = sqrt(abs(1_dp-cosThetaLF(ithLF)**2))
	end do
	write(*,*) "Config Cos(Theta) Step Size: ", cfgCosThStepSize

	! Phi
	! twist about z
	cfgPhiStepSize = (phi_max - phi_min) / real(cfgPhiBins, dp)
	do iphiLF = 1, cfgPhiBins
		phiLF = (iphiLF+0.5_dp)*cfgPhiStepSize
		sinPhiLF(iphiLF) = sin(phiLF)
		cosPhiLF(iphiLF) = cos(phiLF)
	end do
	write(*,*) "Config Phi Step Size: ", cfgPhiStepSize

	! Psi
	! processison about z
	psi_max = 2_dp*pi
	psi_min = 0_dp
	cfgPsiStepSize = (psi_max - psi_min) / real(cfgPsiBins, dp)
	do ipsiLF = 1, cfgPsiBins
		psiLF = (ipsiLF+0.5_dp)*cfgPsiStepSize
		sinPsiLF(ipsiLF) = sin(psiLF)
		cosPsiLF(ipsiLF) = cos(psiLF)
	end do
	write(*,*) "Config Psi Step Size: ", cfgPsiStepSize

	write(*,*) "Computing average force..."

	flush(6)

	! Calculate the average force integral for top half of bisecting plane of cylinder
	do r = 1, cfgRBins ! loop lj--lj distances
		frcSPA = 0_dp; grSPA = 0_dp
		do i = 1, xBins ! full length of cylinder
			do j = 1, zBins ! top half of bisecting plane of cylinder
				rSolv1(1) = -R_axis(r)/2_dp - x_axis(i)
				rSolv1(2) = 0_dp
				rSolv1(3) = -z_axis(j)
				rSolvn(1) = norm2(rSolv1)

				rSolv2(1) = R_axis(r)/2_dp - x_axis(i)
				rSolv2(2) = 0_dp
				rSolv2(3) = -z_axis(j)
				rSolvn(2) = norm2(rSolv2)

				! Populate and wrap the arrays for taking the derivatives for bicubic interpolation at this distance.
				call set_tmp_arrays

				! Loop through orientations of solvent at x(i) and z(j)
				do ithLF = 1, cfgCosThBins
					do iphiLF = 1, cfgPhiBins
						do ipsiLF = 1, cfgPsiBins
							flush(6)
							if ((rSolvn(1) .lt. 1d-6) .or. (rSolvn(2) .lt. 1d-6)) then
								gx = 0_dp ! avoid NaNs in calc_angles
							else
								call calc_angles(ipsiLF, ithLF, iphiLF)
								! todo: spline along r and then pass that g sub array which is 2D
								call bicubic_int(histCosTh,histPhi,histCosThBins,histPhiBins,histCosThStepSize,histPhiStepSize,cosTh_min, &
									& phi_min,gTmp1,cosTh(1),phi(1), gx)
								call bicubic_int(histCosTh,histPhi,histCosThBins,histPhiBins,histCosThStepSize,histPhiStepSize,cosTh_min, &
									& phi_min,gTmp2,cosTh(2),phi(2), gx2)
								gx = exp(gx+gx2)
							end if

							if (gx .gt. 1d-6) then ! if gx == 0 then don't waste time with the rest of the calculation
								call bicubic_int(histCosTh,histPhi,histCosThBins,histPhiBins,histCosThStepSize,histPhiStepSize,cosTh_min, &
									& phi_min,frTmp1,cosTh(1),phi(1), fx(1))
								call bicubic_int(histCosTh,histPhi,histCosThBins,histPhiBins,histCosThStepSize,histPhiStepSize,cosTh_min, &
									& phi_min,fsTmp1,cosTh(1),phi(1), fx(2))
								call bicubic_int(histCosTh,histPhi,histCosThBins,histPhiBins,histCosThStepSize,histPhiStepSize,cosTh_min, &
									& phi_min,ftTmp1,cosTh(1),phi(1), fx(3))

								fNew1 = gx * fx(1) * (-rSolv1(1)/rSolvn(1))	! f*g.R^{hat} for r
								fNew2 = gx * fx(2) * (-sSolv1(1)/sSolv1n)		! f*g.R^{hat} for s
								fNew3 = gx * fx(3) * (-tSolv1(1)/tSolv1n)		! f*g.R^{hat} for t

								fAvg(r) = fAvg(r) + ((fNew1 + fNew2 + fNew3) * z_axis(j)) ! Added radial distance here instead of in fNew

								frcSPA(1,i,j) = frcSPA(1,i,j) + fNew1
								frcSPA(2,i,j) = frcSPA(2,i,j) + fNew2
								frcSPA(3,i,j) = frcSPA(3,i,j) + fNew3
								grSPA(i,j) = grSPA(i,j) + gx
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

		! NOTE : After the fact multiply all elements by 2*pi*density/8/pi/pi (2*2pi*pi/3 (4pi**2)/3 steradians from orientations)
		! 		Number density of chloroform per Angstrom**3 == 0.00750924
		fAvg(r) = fAvg(r)/real(4,dp)/pi*3_dp*density*xzStepSize*xzStepSize*cfgCosThStepSize*cfgPhiStepSize*cfgPsiStepSize
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

	phi(1) = atan2( dot_product(y,-rSolv1) / (norm2(y)*rSolvn(1)), dot_product(x,-rSolv1) / (norm2(y)*rSolvn(1)))
	phi(2) = atan2( dot_product(y,-rSolv2) / (norm2(y)*rSolvn(2)), dot_product(x,-rSolv2) / (norm2(y)*rSolvn(2)))

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
subroutine trilin_interpolate(f, ispa, fP)
	use cfgData
	use histData
	use angleData
	implicit none
	integer :: i, j, ispa, ir0, ir1, ic0, ic1, ip0, ip1
	real(kind=dp) :: f_index, r0, r1, c0, c1, p0, p1, rd, cd, pd, f00, f01, f10, f11, f0, f1, fP, rSolvnInt(2), cosThInt(2), &
		& phiInt(2), f(histDistBins,histCosThBins,histPhiBins)

	fP = 0_dp

	do i = 1, ispa ! Note: for g do SPA, both solutes. for f don't do SPA, only one solute.
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

		f00 = f(ir0,ic0,ip0)*(1-rd) + f(ir1,ic0,ip0)*rd
		f01 = f(ir0,ic0,ip1)*(1-rd) + f(ir1,ic0,ip1)*rd
		f10 = f(ir0,ic1,ip0)*(1-rd) + f(ir1,ic1,ip0)*rd
		f11 = f(ir0,ic1,ip1)*(1-rd) + f(ir1,ic1,ip1)*rd

		f0 = f00*(1-cd) + f10*cd
		f1 = f01*(1-cd) + f11*cd

		fP = fP + (f0*(1-pd) + f1*pd) ! forces are additive and FE is additive.

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
				!print*, (R_axis(cfgRBins)-R_axis(cfgRBins-1))
			else
				! FIXME: is the delta R part of this correct?
				u_dir(cfgRBins-(d-1)) = u_dir(cfgRBins-(d-2)) + fAvg(cfgRBins-(d-1)) * &
					(R_axis(cfgRBins-(d-1))-R_axis(cfgRBins-d))
				! it looks like the first value is getting printed twice. Also, the values might be wrong. Should it be (d-1) and
				! (d-0)? instead of -2 and -1?
				!print*, (R_axis(cfgRBins-(d-1))-R_axis(cfgRBins-d))
			end if
		end do
	end if

end subroutine integrate_force


! write force out and g(r) out to compare against explicit
subroutine write_test_out(r, xBins, zBins)
	use cfgData
	use ctrlData
	implicit none
	integer				:: r, i_f, i, j, xBins, zBins
	character(len=32)	:: temp, filename
	character(len=8)	:: frmt

	i_f = (r-1) * int(RStepSize*10)
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
