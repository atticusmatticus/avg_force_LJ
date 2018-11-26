! USAGE: ./this_file.x -cfg [CONFIGURATION FILE]
!
! calculate the idealized forces to construct an ideal histogram.
!
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!   Modules   !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! data from the config file.
module cfgData
	use prec
	integer						:: nAtoms=5
	integer						:: rBins, cfgCosThBins, cfgPhiBins
	real(kind=dp)				:: rStepSize, rMin, rMax, cosTh_max, cosTh_min, phi_max, phi_min, cfgCosThStepSize, &
		& cfgPhiStepSize, rEq_CH, rEq_CCl, aEq_HCl, AljH, BljH, AljC, BljC, AljCl, BljCl

end module cfgData


! angle arrays
module arrayData
	use prec
	real(kind=dp), allocatable	:: rAxis(:), pos_i(:,:), pos_r(:,:), pos_t(:,:), Alj(:), Blj(:), theta(:), phi(:), frcHist(:,:,:,:)
end module arrayData


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module subs
	use prec
	implicit none
	contains

! parse commandline for relevant files.
subroutine parse_command_line(cfgFile)
	implicit none
	character(len=64) 	:: cfgFile
	character(len=16) 	:: arg
	integer 			:: i
	logical 			:: cfgFileFlag

	cfgFileFlag = .false.
	i=1
	do
		call get_command_argument(i,arg)
		select case (arg)

		case ('-cfg')
			i = i+1
			call get_command_argument(i,cfgFile)
			cfgFileFlag=.true.
			write(*,*) "cfg File: ", cfgFile
		case default
			write(*,*) 'Unrecognized command-line option: ', arg
			write(*,*) 'Usage: compute_avgForce.x -cfg [cfg file]'
			stop

		end select
		i = i+1
		if (i.ge.command_argument_count()) exit
	end do

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
	logical 			:: outFileFlag, rStepSizeFlag, rMinFlag, rMaxFlag, thetaBinsFlag, phiBinsFlag, rEq_CHFlag, rEq_CClFlag, &
		& aEq_HClFlag, AljHFlag, BljHFlag, AljCFlag, BljCFlag, AljClFlag, BljClFlag

	outFileFlag = .false.
	rStepSizeFlag = .false.
	rMinFlag = .false.
	rMaxFlag = .false.
	thetaBinsFlag = .false.
	phiBinsFlag = .false.
	rEq_CHFlag = .false.
	rEq_CClFlag = .false.
	aEq_HClFlag = .false.
	AljHFlag = .false.
	BljHFlag = .false.
	AljCFlag = .false.
	BljCFlag = .false.
	AljClFlag = .false.
	BljClFlag = .false.

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
			else if (firstWord .eq. "rStepSize") then
				read(line,*) rStepSize
				write(*,*) "Solvent Distance Step Size:	", rStepSize
				rStepSizeFlag = .true.
			else if (firstWord .eq. "rMin") then
				read(line,*) rMin
				write(*,*) "r - Min:	", rMin
				rMinFlag = .true.
			else if (firstWord .eq. "rMax") then
				read(line,*) rMax
				write(*,*) "r - Max:	", rMax
				rMaxFlag = .true.
			else if (firstWord .eq. "theta_bins") then
				read(line,*) cfgCosThBins
				write(*,*) "Theta Bins:	", cfgCosThBins
				thetaBinsFlag = .true.
			else if (firstWord .eq. "phi_bins") then
				read(line,*) cfgPhiBins
				write(*,*) "Phi Bins:	", cfgPhiBins
				phiBinsFlag = .true.
			else if (firstWord .eq. "rEq_CH") then
				read(line,*) rEq_CH
				write(*,*) "Eq. Bond Dist - CH:	", rEq_CH
				rEq_CHFlag = .true.
			else if (firstWord .eq. "rEq_CCl") then
				read(line,*) rEq_CCl
				write(*,*) "Eq. Bond Dist - CCl:	", rEq_CCl
				rEq_CClFlag = .true.
			else if (firstWord .eq. "aEq_HCl") then
				read(line,*) aEq_HCl
				write(*,*) "Eq. Angle - HCl:	", aEq_HCl
				aEq_HClFlag = .true.
			else if (firstWord .eq. "AljH") then
				read(line,*) AljH
				write(*,*) "A lj Coefficient H:	", AljH
				AljHFlag = .true.
			else if (firstWord .eq. "BljH") then
				read(line,*) BljH
				write(*,*) "B lj Coefficient H:	", BljH
				BljHFlag = .true.
			else if (firstWord .eq. "AljC") then
				read(line,*) AljC
				write(*,*) "A lj Coefficient C:	", AljC
				AljCFlag = .true.
			else if (firstWord .eq. "BljC") then
				read(line,*) BljC
				write(*,*) "B lj Coefficient C:	", BljC
				BljCFlag = .true.
			else if (firstWord .eq. "AljCl") then
				read(line,*) AljCl
				write(*,*) "A lj Coefficient Cl:	", AljCl
				AljClFlag = .true.
			else if (firstWord .eq. "BljCl") then
				read(line,*) BljCl
				write(*,*) "B lj Coefficient Cl:	", BljCl
				BljClFlag = .true.
			end if
		end if
	end do
	close(20)

	if (rStepSizeFlag.eqv..false.) then
		write(*,*) "Config file must have a 'rStepSize' value"
		stop
	end if
	if (outFileFlag.eqv..false.) then
		write(*,*) "Config file must have a 'out_file' value"
		stop
	end if
	if (rMinFlag.eqv..false.) then
		write(*,*) "Config file must have a 'rMin' value"
		stop
	end if
	if (rMaxFlag.eqv..false.) then
		write(*,*) "Config file must have a 'rMax' value"
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
	if (rEq_CHFlag.eqv..false.) then
		write(*,*) "Config file must have a 'rEq_CH' value"
		stop
	end if
	if (rEq_CClFlag.eqv..false.) then
		write(*,*) "Config file must have a 'rEq_CCl' value"
		stop
	end if
	if (aEq_HClFlag.eqv..false.) then
		write(*,*) "Config file must have a 'aEq_HCl' value"
		stop
	end if

	flush(6)

end subroutine read_cfg


! generate starting configuration of CL3 @ theta=0, phi=0 and LJ coefficient arrays
subroutine initial_geom
	use constants
	use cfgData
	use arrayData
	implicit none
	integer			:: i, j
	real(kind=dp)	:: psi

	allocate( pos_i(nAtoms,3), pos_r(nAtoms,3), pos_t(nAtoms,3), Alj(nAtoms), Blj(nAtoms) )

	pos_i = 0_dp
	do i = 1, nAtoms
		pos_i(i,3) = 1_dp
	end do

	! 1 is H
	pos_i(1,:) = [ 0_dp, -1_dp, 0_dp ]
	pos_i(1,:) = pos_i(1,:) * rEq_CH
	Alj(1) = AljH
	Blj(1) = BljH

	! 2 is C
	pos_i(2,:) = [ 0_dp, 0_dp, 0_dp ]
	Alj(2) = AljC
	Blj(2) = BljC

	! 3-5 are Cl
	! FIXME, BUG: the Cl-C-Cl angle isn't correct.
	psi = 0_dp
	do i = 3, nAtoms
		pos_i(i,1) = -dsin(psi)*dsin(aEq_HCl)
		pos_i(i,2) = -dcos(aEq_HCl)
		pos_i(i,3) = -dsin(aEq_HCl)*dcos(psi)

		pos_i(i,:) = pos_i(i,:) * rEq_CCl

		psi = psi + 2_dp*pi/3_dp

		Alj(i) = AljCl
		Blj(i) = BljCl
	end do
	print*, pos_i(1,:)
	print*, pos_i(2,:)
	print*, pos_i(3,:)
	print*, pos_i(4,:)
	print*, pos_i(5,:)
	print*, dacos(dot_product(pos_i(1,:),pos_i(3,:))/norm2(pos_i(1,:))/norm2(pos_i(3,:)))
	print*, dacos(dot_product(pos_i(1,:),pos_i(4,:))/norm2(pos_i(1,:))/norm2(pos_i(4,:)))
	print*, dacos(dot_product(pos_i(1,:),pos_i(5,:))/norm2(pos_i(1,:))/norm2(pos_i(5,:)))

	print*, dacos(dot_product(pos_i(3,:),pos_i(4,:))/norm2(pos_i(3,:))/norm2(pos_i(4,:)))
	print*, dacos(dot_product(pos_i(3,:),pos_i(5,:))/norm2(pos_i(3,:))/norm2(pos_i(5,:)))
	print*, dacos(dot_product(pos_i(4,:),pos_i(5,:))/norm2(pos_i(4,:))/norm2(pos_i(5,:)))

end subroutine initial_geom


! compute idealized force histogram
subroutine compute_force_hist
	use constants
	use functions
	use cfgData
	use arrayData
	implicit none
	integer						:: i, ir, ith, iphi
	real(kind=dp),dimension(3)	:: r_vec, s_vec, t_vec, force_vec

	rBins = int( (rMax-rMin)/rStepSize + 1 )
	write(*,*) "Number of r Bins: ", rBins

	! allocate force histogram array
	allocate( frcHist(3, rBins, cfgCosThBins, cfgPhiBins) )
	frcHist = 0_dp

	! Distance Axis
	allocate( rAxis(rBins) , theta(cfgCosThBins), phi(cfgPhiBins) )
	rAxis = 0_dp

	do i = 1, rBins
		rAxis(i) = rMin + ((i-1) * rStepSize + rStepSize/2_dp)
	end do

	! ANGLES
	! Theta
	! tilt off of z
	cosTh_max = 1_dp
	cosTh_min = -1_dp
	cfgCosThStepSize = (cosTh_max - cosTh_min) / real(cfgCosThBins, dp)
	do ith = 1, cfgCosThBins
		theta(ith) = dacos(((ith-1)*cfgCosThStepSize + cfgCosThStepSize/2_dp) + cosTh_min)
	end do
	write(*,*) "Config Cos(Theta) Step Size: ", cfgCosThStepSize

	! Phi
	! twist about z
	phi_max = 2_dp*pi/3_dp ! note: this is sufficient for a molecule with C3 symmetry.
	phi_min = 0_dp
	cfgPhiStepSize = (phi_max - phi_min) / real(cfgPhiBins, dp)
	do iphi = 1, cfgPhiBins
		phi(iphi) = phi_min + ((iphi-1)*cfgPhiStepSize + cfgPhiStepSize/2_dp)
	end do
	write(*,*) "Config Phi Step Size: ", cfgPhiStepSize

	flush(6)

	! Calculate the average force integral for top half of bisecting plane of cylinder
	r_vec = [0_dp,-1_dp,0_dp]
	do ith = 1, cfgCosThBins
		do iphi = 1, cfgPhiBins
			! rotate initial geom
			do i = 1, nAtoms
				pos_r(i,:) = rotate_x(theta(ith), rotate_y(phi(iphi),pos_i(i,:)))
			end do
			t_vec = cross_product(r_vec, pos_r(1,:))
			t_vec = t_vec / norm2(t_vec)
			s_vec = cross_product(t_vec, r_vec)
			pos_t = pos_r ! the position vector that will be translated
			do ir = 1, rBins
				force_vec = 0_dp
				do i = 1, nAtoms
					! translate along y
					pos_t(i,2) = pos_r(i,2) + rAxis(ir)
					! force calculation
					force_vec = force_vec + lj_force(i, pos_t(i,:))
				end do !i atoms
				frcHist(1,ir,ith,iphi) = force_vec(2)					! <f.r>
				frcHist(2,ir,ith,iphi) = dot_product(force_vec, s_vec)	! <f.s>
				frcHist(3,ir,ith,iphi) = dot_product(force_vec, t_vec)	! <f.t>
				! XXX: plot positions instead for troubleshooting 
				!frcHist(1,ir,ith,iphi) = pos_t(1,2)
				!frcHist(2,ir,ith,iphi) = pos_t(1,3)
				!frcHist(3,ir,ith,iphi) = pos_t(2,2)
				!frcHist(4,ir,ith,iphi) = pos_t(2,3)
				!frcHist(5,ir,ith,iphi) = pos_t(3,2)
				!frcHist(6,ir,ith,iphi) = pos_t(3,3)
				!frcHist(7,ir,ith,iphi) = pos_t(4,2)
				!frcHist(8,ir,ith,iphi) = pos_t(4,3)
				!frcHist(9,ir,ith,iphi) = pos_t(5,2)
				!frcHist(10,ir,ith,iphi) = pos_t(5,3)
			end do !r
		end do !phi
	end do !theta
end subroutine compute_force_hist


! calculate LJ force, (-gradient of PE), for a particular atom type on the solute from a distance 'pos'
function lj_force(i, pos)
	use cfgData
	use arrayData
	implicit none
	integer,		intent(in)	:: i								! inputs not to be changed
	real(kind=dp),	intent(in)	:: pos(3)							! inputs not to be changed
	real(kind=dp)				:: posDist2, posDist6, lj_force(3)	! output

	posDist2 = pos(1)**2 + pos(2)**2 + pos(3)**2
	posDist6 = posDist2**(-3)

	lj_force = ( posDist6*( 12_dp*Alj(i)*posDist6 - 6_dp*Blj(i) )/posDist2 ) * pos
end function lj_force


! write output file
subroutine write_output(outFile)
	use cfgData
	use arrayData
	implicit none
	character(len=64) 	:: outFile
	integer 			:: i, j, k

	open(35,file=outFile,status='replace')
	write(6,*) "Writing output file:	", outFile
	write(35,*) "# 1.	r Distance"
	write(35,*) "# 2.	cos(Theta)"
	write(35,*) "# 3.	Phi"
	write(35,*) "# 4.	<f.r>"
	write(35,*) "# 5.	<f.s>"
	write(35,*) "# 6.	<f.t>"
	do i = 1, rBins
		do j = 1, cfgCosThBins
			do k = 1, cfgPhiBins
				write(35,899) rAxis(i), cos(theta(j)), phi(k), frcHist(1,i,j,k), frcHist(2,i,j,k), frcHist(3,i,j,k)!,&
				!&frcHist(4,i,j,k), frcHist(5,i,j,k), frcHist(6,i,j,k), frcHist(7,i,j,k), frcHist(8,i,j,k), frcHist(9,i,j,k),&
				!&frcHist(10,i,j,k)
			end do
		end do
	end do
	close(35)

	flush(6)

899		format (6(1x,es14.7)) ! scientific format
end subroutine write_output

end module subs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!  Main Program  !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program ideal_force
	use subs
	implicit none
	character(len=64) 	:: cfgFile, outFile
	real(kind=dp) 		:: omp_get_wtime, ti, tf

	ti = omp_get_wtime()

	! make list of average direct force from 'collapsed' file.
	call parse_command_line(cfgFile)

	! read config file
	call read_cfg(cfgFile, outFile)

	! calculate initial geometry and initialize arrays
	call initial_geom

	! compute average force integral.
	call compute_force_hist
	
	! write PMF output file
	call write_output(outFile)
	
	! Write time taken to finish calculation.
	tf = omp_get_wtime()
	write(*,*) "Total time elapsed: ", tf-ti, "seconds"

	flush(6)
end program ideal_force
