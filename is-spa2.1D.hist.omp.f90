! USAGE: ./this_file.x -cfg [CONFIGURATION FILE]
!
!
!                      ^
!                     z|
!                      |
!                      |      ^
!                      |   /y
!                      | /
!      <_______O_______|_______O_______>
!      -x        1               2      +x
!               <-----R
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!   Modules   !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! data for the density and force tables.
module histData
   use prec
   real(kind=dp),allocatable :: histDist(:), histCosTh(:), histPhi(:), g(:,:,:,:), fLJr(:,:,:,:), fLJs(:,:,:,:), fLJt(:,:,:,:), &
      & fCr(:,:,:,:), fCs(:,:,:,:), fCt(:,:,:,:), gc(:,:,:,:)
   real(kind=dp),allocatable :: g2D(:,:,:), fLJr2D(:,:,:), fLJs2D(:,:,:), fCr2D(:,:,:), fCs2D(:,:,:)
   real(kind=dp),allocatable :: g1D(:,:), fLJr1D(:,:), fCr1D(:,:), emf1D(:,:), g1D2(:,:), fLJr1D2(:,:), fCr1D2(:,:), emf1D2(:,:), &
      & avgCosTh(:,:)
   real(kind=dp) :: histDistStepSize, histCosThStepSize, histPhiStepSize
   integer :: histDistBins, histCosThBins, histPhiBins
   real(kind=dp),allocatable :: longRange(:) !debug

end module histData

! data from the config file.
module cfgData
   use prec
   use constants
   real(kind=dp),allocatable :: x_axis(:), z_axis(:), R_axis(:), fAvg(:,:), u_dir(:,:)
   real(kind=dp) :: RStepSize, xzStepSize, R_min, R_max, xz_range, cfgCosThStepSize, cfgPsiStepSize, T, cut, offset, soluteChg(2), radius
   character(len=8) :: c_explicit_R
   integer :: cfgRBins, cfgCosThBins, cfgPhiBins, cfgPsiBins, nThreads
!
   integer :: xBins, zBins
   real(kind=dp) :: density = 0.00739753994553948_dp ! numerical density of chloroforms per Angstrom**3
   real(kind=dp) :: dipole_moment = 0.23026679557844498_dp ! magnitude of dipole moment chloroform in q_e * AA
   real(kind=dp) :: dielectric = 2.3507059228494147_dp ! reduced dielectric constant of chloroform
   real(kind=dp) :: cosTh_max = 1_dp
   real(kind=dp) :: cosTh_min = -1_dp
   real(kind=dp) :: phi_max = pi/3_dp
   real(kind=dp) :: phi_min = 0_dp
   real(kind=dp) :: psi_max = 2_dp*pi
   real(kind=dp) :: psi_min = 0_dp
end module cfgData

! data for calculating cosTh value.
module angleData
   use prec
   real(kind=dp),allocatable :: sinThetaLF(:), cosThetaLF(:), sinPhiLF(:), cosPhiLF(:), sinPsiLF(:), cosPsiLF(:)
   real(kind=dp) :: rSolv1(3), rSolv2(3), rSolvn(2), cosTh(2), phi(2), sSolv1(3), tSolv1(3), sSolv1n, tSolv1n

   !$omp THREADPRIVATE( rSolv1, rSolv2, rSolvn, sSolv1, sSolv1n, tSolv1, tSolv1n, cosTh, phi )
end module angleData

! testing arrays for force and g(r)
module ctrlData
   use prec
   real(kind=dp),allocatable :: frcSPA(:,:,:,:), grSPA(:,:), explicitDist(:)
   integer :: crdLines
   logical :: explicit_R
end module ctrlData


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!  Main Program  !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program compute_avgForce
   use prec
   implicit none
   character(len=128) :: histFile, cfgFile, outFile
   real(kind=dp) :: omp_get_wtime, ti, tf, seconds
   integer :: hours, minutes

   ti = omp_get_wtime()

   ! make list of average direct force from 'collapsed' file.
   call parse_command_line(cfgFile)

   ! read config file
   call read_cfg(cfgFile, histFile, outFile)

   ! make list of average direct force from 'collapsed' file.
   call make_hist_table(histFile)

   ! Now that we have the relevant information spline the g and f arrays along r.
   call spline_hist_array

   ! read in LJ--LJ dist array from file
   call R_list

   ! setup for computing the average force integral.
   call setup_compute_avg_force

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
   write(*,'(a,i4,a,i2,a,f6.3,a)') "Total time elapsed:   ", hours, "h ", minutes, "m ", seconds, "s"

   flush(6)

end program compute_avgForce


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse commandline for relevant files.
subroutine parse_command_line(cfgFile)
   implicit none
   character(len=128) :: cfgFile
   character(len=16) :: arg
   integer :: i
   logical :: cfgFileFlag, cfgExist

   cfgFileFlag = .false.
   cfgExist = .false.
   i=1
   do
      call get_command_argument(i,arg)
      select case (arg)

      case ('-cfg')
         i = i+1
         call get_command_argument(i,cfgFile)
         cfgFileFlag=.true.
         INQUIRE(FILE=cfgFile, EXIST=cfgExist)
         write(*,*) 'Config File:            ', cfgFile
         write(*,*) 'Config File Exists:         ', cfgExist
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

   ! 'ERROR STOP' if either file doesn't exist
   if (cfgExist.eqv..false.) then
      write(*,*) 'cfg file does not exist'
      error stop
   end if

   flush(6)

end subroutine parse_command_line


! read python cfg file for g(r) parameters
subroutine read_cfg(cfgFile, histFile, outFile)
   use cfgData
   implicit none
   character(len=128) :: cfgFile, histFile, outFile
   character(len=256) :: line
   character(len=32) :: firstWord, sep
   integer :: ios
   logical :: outFileFlag, histFileFlag, histExist, RstepSizeFlag, xzStepSizeFlag, RmaxFlag, RminFlag, xzRangeFlag, &
      & thetaBinsFlag, phiBinsFlag, psiBinsFlag, c_explicit_RFlag, TFlag, cutFlag, radiusFlag, offsetFlag, nThreadsFlag, &
      & soluteChgFlag

   histFileFlag = .false.;  histExist = .false.
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
   cutFlag = .false.
   radiusFlag = .false.
   offsetFlag = .false.
   nThreadsFlag = .false.
   soluteChgFlag = .false.

   ios = 0

   open(20,file=cfgFile)
   do while(ios>=0)
      read(20,'(a)',IOSTAT=ios) line
      call split(line,'=',firstWord, sep)
      if (line .ne. "") then
         if (firstWord .eq. "hist_file") then
            read(line,'(a)') histFile
            write(*,*) "Histogram File:            ", histFile
            histFileFlag = .true.
            INQUIRE(FILE=histFile, EXIST=histExist) ! check if it exists
         else if (firstWord .eq. "out_file") then
            read(line,*) outFile
            write(*,*) "Output File:            ", outFile
            outFileFlag = .true.
         else if (firstWord .eq. "RStepSize") then
            read(line,*) RStepSize
            write(*,*) "PMF Step Size:            ", RStepSize
            RstepSizeFlag = .true.
         else if (firstWord .eq. "xzStepSize") then
            read(line,*) xzStepSize
            write(*,*) "Solvent Grid Step Size:      ", xzStepSize
            xzStepSizeFlag = .true.
         else if (firstWord .eq. "explicit_R") then
            read(line,*) c_explicit_R
            write(*,*) "Use Explicit R Values:         ", c_explicit_R
            c_explicit_RFlag = .true.
         else if (firstWord .eq. "R_max") then
            read(line,*) R_max
            write(*,*) "R Maximum Value:         ", R_max
            RmaxFlag = .true.
         else if (firstWord .eq. "R_min") then
            read(line,*) R_min
            write(*,*) "R Minimum Value:         ", R_min
            RminFlag = .true.
         else if (firstWord .eq. "xz_range") then
            read(line,*) xz_range
            write(*,*) "XZ - Range:            ", xz_range
            xzRangeFlag = .true.
         else if (firstWord .eq. "theta_bins") then
            read(line,*) cfgCosThBins
            write(*,*) "Theta Bins:         ", cfgCosThBins
            thetaBinsFlag= .true.
         else if (firstWord .eq. "phi_bins") then
            read(line,*) cfgPhiBins
            write(*,*) "Phi Bins:         ", cfgPhiBins
            phiBinsFlag= .true.
         else if (firstWord .eq. "psi_bins") then
            read(line,*) cfgPsiBins
            write(*,*) "Psi Bins:         ", cfgPsiBins
            psiBinsFlag= .true.
         else if (firstWord .eq. "temperature") then
            read(line,*) T
            write(*,*) "Temperature (K):         ", T
            TFlag= .true.
         else if (firstWord .eq. "bicubic_cutoff") then
            read(line,*) cut
            write(*,*) "Bicubic/Bilinear Cutoff:      ", cut
            cutFlag= .true.
         else if (firstWord .eq. "solute_radius") then
            read(line,*) radius
            write(*,*) "Solute radius:      ", radius
            radiusFlag= .true.
         else if (firstWord .eq. "offset") then
            read(line,*) offset
            write(*,*) "Solvent offset distance:      ", offset
            offsetFlag= .true.
         else if (firstWord .eq. "num_threads") then
            read(line,*) nThreads
            write(*,*) "Number of Parallel Threads:      ", nThreads
            nThreadsFlag= .true.
         else if (firstWord .eq. "solute_charge") then
            read(line,*) soluteChg(1)
            soluteChg(2) = -soluteChg(1)
            write(*,*) "Solute Charges:      ", soluteChg
            soluteChgFlag= .true.
         end if
      end if
   end do
   close(20)

   if (histFileFlag.eqv..false.) then
      write(*,*) "Config file must have a 'hist_file' value"
      stop
   end if
   if (histExist.eqv..false.) then
      write(*,*) "Config file must point to a 'hist_file' that exists: ", histFile, " doesn't exist."
      stop
   end if
   if (outFileFlag.eqv..false.) then
      write(*,*) "Config file must have a 'out_file' value"
      stop
   end if
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
   if (cutFlag.eqv..false.) then
      write(*,*) "Config file must have a 'bicubic_cuttoff' value"
      stop
   end if
   if (radiusFlag.eqv..false.) then
      write(*,*) "Config file must have a 'solute_radius' value"
      stop
   end if
   if (offsetFlag.eqv..false.) then
      write(*,*) "Config file must have a 'offset' value"
      stop
   end if
   if (nThreadsFlag.eqv..false.) then
      write(*,*) "Config file must have a 'num_threads' value"
      stop
   end if
   if (soluteChgFlag.eqv..false.) then
      write(*,*) "Config file must have a 'solute_charge' value"
      stop
   end if

   flush(6)

end subroutine read_cfg


! read force file and make a lookup table.
subroutine make_hist_table(histFile)
   use histData; use cfgData
   implicit none
   character(len=128) :: histFile
   character(len=512) :: line
   integer :: ios, ios2, i, j, k, nHistLines
   real(kind=dp),allocatable :: histTmp(:,:)

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
   write(*,*) "nHistLines", nHistLines

   allocate( histTmp(19,nHistLines) )

   ! populate hist arrays
   ios = 0; i = 1
   open(20,file=histFile)
   ! read file ignoring comment lines at the beginning
   do while(ios>=0)
      read(20,'(a)',IOSTAT=ios) line
      if ((line(1:1) .ne. "#") .and. (ios .ge. 0)) then
         !                     r         cos(Th)         phi/3
         read(line,*) histTmp(1,i), histTmp(2,i), histTmp(3,i), &
            !        g(r)+         g(r)-
            & histTmp(4,i), histTmp(5,i), &
               !     <fLJ.r>+       <fLJ.s>+    <fLJ.t>+
               & histTmp(6,i), histTmp(7,i), histTmp(8,i), &
                  !     <fLJ.r>-       <fLJ.s>-       <fLJ.t>-
                  & histTmp(9,i), histTmp(10,i), histTmp(11,i), &
                     !     <fC.r>+       <fC.s>+    <fC.t>+
                     & histTmp(12,i), histTmp(13,i), histTmp(14,i), &
                        !     <fC.r>-       <fC.s>-       <fC.t>-
                        & histTmp(15,i), histTmp(16,i), histTmp(17,i), &
                           !       gc(r)+       gc(r)-
                           & histTmp(18,i), histTmp(19,i)
         i = i + 1
      end if
   end do
   close(20)

   ! Unique value determination
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
         if (( histTmp(2,i) .gt. (histTmp(2,1)-1d-6) ) .and. ( histTmp(2,i) .lt. (histTmp(2,1)+1d-6) ) .and. ( ios .eq. 0 ) .and. &
            & ( ios2 .eq. 1 )) then
            ! note: this statement will trigger when i = histCosThBins+1 because it finds the first repeated element
            histCosThBins = (i - 1)/histPhiBins
            ios = 1
         end if
         if (( histTmp(3,i) .gt. (histTmp(3,1)-1d-6) ) .and. ( histTmp(3,i) .lt. (histTmp(3,1)+1d-6) ) .and. ( ios2 .eq. 0 )) then
            ! note: this statement will trigger when i = histPhiBins+1 because it finds the first repeated element
            histPhiBins = i - 1
            ios2 = 1
         end if
      end if
   end do
   write(*,*) "Histogram Distance Bins:   ", histDistBins
   write(*,*) "Histogram Cosine Theta Bins:   ", histCosThBins
   write(*,*) "Histogram Phi Bins:      ", histPhiBins

   allocate( histDist(histDistBins), histCosTh(histCosThBins), histPhi(histPhiBins), g(2,histDistBins,histCosThBins,histPhiBins), &
      & fLJr(2,histDistBins,histCosThBins,histPhiBins), fLJs(2,histDistBins,histCosThBins,histPhiBins), &
      & fLJt(2,histDistBins,histCosThBins,histPhiBins), fCr(2,histDistBins,histCosThBins,histPhiBins), &
      & fCs(2,histDistBins,histCosThBins,histPhiBins), fCt(2,histDistBins,histCosThBins,histPhiBins), &
      & gc(2,histDistBins,histCosThBins,histPhiBins) )
   
   ! populate arrays that will be used in the rest of the calculation from temp array
   do i = 1, histDistBins      ! the values written out from python script are at half-bin distances
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
            g(1,i,j,k)  = histTmp( 4, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! g(r,cos,phi)+ currently g
            g(2,i,j,k)  = histTmp( 5, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! g(r,cos,phi)- currently g
            fLJr(1,i,j,k) = histTmp( 6, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <fLJ.r>(r,cos,phi) +
            fLJs(1,i,j,k) = histTmp( 7, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <fLJ.s>(r,cos,phi) +
            fLJt(1,i,j,k) = histTmp( 8, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <fLJ.t>(r,cos,phi) +
            fLJr(2,i,j,k) = histTmp( 9, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <fLJ.r>(r,cos,phi) -
            fLJs(2,i,j,k) = histTmp(10, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <fLJ.s>(r,cos,phi) -
            fLJt(2,i,j,k) = histTmp(11, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <fLJ.t>(r,cos,phi) -
            fCr(1,i,j,k) = histTmp(12, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <fC.r>(r,cos,phi) +
            fCs(1,i,j,k) = histTmp(13, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <fC.s>(r,cos,phi) +
            fCt(1,i,j,k) = histTmp(14, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <fC.t>(r,cos,phi) +
            fCr(2,i,j,k) = histTmp(15, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <fC.r>(r,cos,phi) -
            fCs(2,i,j,k) = histTmp(16, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <fC.s>(r,cos,phi) -
            fCt(2,i,j,k) = histTmp(17, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <fC.t>(r,cos,phi) -
            gc(1,i,j,k) = histTmp(18, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! gc(r,cos,phi) +
            gc(2,i,j,k) = histTmp(19, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! gc(r,cos,phi) +
         end do
      end do
   end do

   histDistStepSize = histDist(2) - histDist(1)
   write(*,*) "Histogram Distance Step Size:      ", histDistStepSize
   histCosThStepSize = histCosTh(2) - histCosTh(1)
   write(*,*) "Histogram Cosine Theta Step Size:   ", histCosThStepSize
   histPhiStepSize = histPhi(2) - histPhi(1)
   write(*,*) "Histogram Phi Step Size:      ", histPhiStepSize 

   flush(6)
   !print*, 'look4', g(1,30,1,1), g(2,30,1,1)

end subroutine make_hist_table


! calculate reduced mean field emf1D(r) from g(r,cosTh,phi)
subroutine reduced_mean_field(ia)
   use histData; use functions
   implicit none
   integer :: ia, ir, ith, iphi
   real(kind=dp) :: boltz, boltz_sum
   real(kind=dp),allocatable :: u02D(:,:)
!   real(kind=dp) :: x, y1, y2 !debug

   write(*,*) 'Calculating reduced mean field...'; flush(6)

   allocate( u02D(histDistBins,histCosThBins) )
   u02D = 0_dp
   u02D = minval(-g(ia,:,:,:),dim=3)

   avgCosTh(ia,:) = 0_dp
   do ir = 1, histDistBins
      boltz_sum = 0_dp
      do ith = 1, histCosThBins
         do iphi = 1, histPhiBins
            boltz = exp(g(ia,ir,ith,iphi) )!- u02D(ir,ith)) ! i think this should be a + sign but the data looks wrong when it is.
            avgCosTh(ia,ir) = avgCosTh(ia,ir) + (boltz * histCosTh(ith))
            boltz_sum = boltz_sum + boltz
         end do
      end do
      if (abs(boltz_sum).lt.1.0d-6) then
         avgCosTh(ia,ir) = 0_dp
      else
         avgCosTh(ia,ir) = avgCosTh(ia,ir) / boltz_sum
      end if
   end do

   do ir = 1, histDistBins
      emf1D(ia,ir) = (3*avgCosTh(ia,ir) - avgCosTh(ia,ir)*(6*avgCosTh(ia,ir)**2 + avgCosTh(ia,ir)**4 - 2*avgCosTh(ia,ir)**6)/5.) &
         & / (1-avgCosTh(ia,ir)**2) ! inverse Langevin function
   end do

   if (ia.eq.2) then !debug
      open(99,file='emf.out',status='replace')
      write(99,*) '#  1. Distance'
      write(99,*) '#  2. emf +'
      write(99,*) '#  3. <p> +'
      write(99,*) '#  4. emf -'
      write(99,*) '#  5. <p> -'
      !do ir = 1, histDistBins*100
         !x = ir/100.0/10.0
         !call splint(histDist,emf1D(1,:),emf1D2(1,:),histDistBins,x, y1)
         !call splint(histDist,emf1D(2,:),emf1D2(2,:),histDistBins,x, y2)
         !write(99,*) x, y1, y2
      !end do
      do ir = 1, histDistBins
         write(99,*) histDist(ir), emf1D(1,ir), avgCosTh(1,ir), emf1D(2,ir), avgCosTh(2,ir)
      end do
      close(99)
   end if
end subroutine reduced_mean_field


! spline the r dimension of each theta phi stack and then average over phi for 2D
subroutine spline_hist_array
   use constants; use functions; use histData; use cfgData; use idealSolv
   implicit none
   integer :: ia, ir, ith, iphi, imin, igo, igo1, igo2
   real(kind=dp) :: norm_factor, boltz, boltz_sum
   real(kind=dp),allocatable :: idealHist(:,:,:,:,:), idealHist2D(:,:,:,:), u02D(:,:), idealHist1D(:,:,:), u01D(:)
   integer,allocatable :: ispline(:,:), ispline2D(:)
   integer :: ispline1D
   !real(kind=dp)   :: xx, yy   !debug
   !integer  :: rTmp !debug

   write(*,*) 'Editing input histogram arrays with ideal arrays in 3D...'

   ! Calculate a 4D array idealHist(lj+/lj-,g/f,r,th,phi)
   allocate( idealHist(2,13,histDistBins,histCosThBins,histPhiBins), ispline(histCosThBins,histPhiBins), ispline2D(histCosThBins) )
   idealHist = 0_dp; ispline = 0_dp
   call ideal_CL3(histDistBins,histDistStepSize,histCosThBins,cosTh_min,cosTh_max,histPhiBins,phi_min,phi_max,radius,offset,T, &
      & soluteChg, idealHist)
   !print*, 'look1', idealHist(1,1,30,1,1), idealHist(2,1,30,1,1) !debug

   ! Spline the log(g) and force arrays using ideal values for the slopes at small r. This populates the second derivative arrays.
   ! This requires ideal values that have been averaged over phi.
   allocate( idealHist2D(2,9,histDistBins,histCosThBins) )
   idealHist2D = 0_dp
   call ideal_3D_to_2D(idealHist,histDistBins,histCosThBins,histPhiBins, idealHist2D)

   allocate( idealHist1D(2,5,histDistBins) )
   idealHist1D = 0_dp
   call ideal_2D_to_1D(idealHist2D,histDistBins,histCosThBins, idealHist1D)


   ! Allocate explicit data arrays for 3D -> 2D transformation
   allocate( g2D(2,histDistBins,histCosThBins), fLJr2D(2,histDistBins,histCosThBins), fLJs2D(2,histDistBins,histCosThBins), &
      & fCr2D(2,histDistBins,histCosThBins), fCs2D(2,histDistBins,histCosThBins), u02D(histDistBins,histCosThBins) )

   g2D = 0_dp; fLJr2D = 0_dp; fLJs2D = 0_dp; fCr2D = 0_dp; fCs2D = 0_dp; u02D = 0_dp

   ! Allocate explicit data arrays for 2D -> 1D transformation
   allocate( g1D(2,histDistBins), fLJr1D(2,histDistBins), fCr1D(2,histDistBins), u01D(histDistBins) )
   g1D = 0_dp; fLJr1D = 0_dp; fCr1D = 0_dp; u01D = 0_dp

   allocate( g1D2(2,histDistBins), fLJr1D2(2,histDistBins), fCr1D2(2,histDistBins), emf1D(2,histDistBins), &
      & emf1D2(2,histDistBins), avgCosTh(2,histDistBins) )
   g1D2 = 0_dp; fLJr1D2 = 0_dp; fCr1D2 = 0_dp; emf1D = 0_dp; emf1D2 = 0_dp


   ! Edit the input hist arrays to more smoothly transition to -/+ infinity with the help of idealHist.
   do ia = 1, 2 ! which solute
      do ith = 1, histCosThBins
         do iphi = 1, histPhiBins
            imin = 0
            ! Normalization factor or each theta phi array.
            norm_factor = gc(ia,histDistBins,ith,iphi)/(g(ia,histDistBins,ith,iphi)*4*pi*histDist(histDistBins)**2)
            ! Find the first non-zero g(r) bin for each theta/phi array and set 'imin' to that 'ir' index
find:       do ir = 1, histDistBins
               if (g(ia,ir,ith,iphi).gt.1d-6) then
                  imin = ir
                  exit find
               end if
            end do find
             
            ! Note: Add the ideal values to bins with no sampling. And half counts to bins that probably should have had sampling.
            igo = 0
            do ir = histDistBins, 1, -1
               if (ir.ge.imin) then
                  if (g(ia,ir,ith,iphi).gt.1d-6) then
                     g(ia,ir,ith,iphi) = log(g(ia,ir,ith,iphi)) ! g is ln(g) now
                  else  ! note: This is a zero bin where there probably should have been something. So put a half count in.
                     g(ia,ir,ith,iphi) = log(real(0.5,dp)/(4*pi*histDist(ir)**2)/norm_factor)
                     fLJr(ia,ir,ith,iphi) = idealHist(ia, 2,ir,ith,iphi)
                     fLJs(ia,ir,ith,iphi) = idealHist(ia, 3,ir,ith,iphi)
                     fLJt(ia,ir,ith,iphi) = idealHist(ia, 4,ir,ith,iphi)
                     fCr(ia,ir,ith,iphi)  = idealHist(ia, 8,ir,ith,iphi)
                     fCs(ia,ir,ith,iphi)  = idealHist(ia, 9,ir,ith,iphi)
                     fCt(ia,ir,ith,iphi)  = idealHist(ia,10,ir,ith,iphi)
                  end if
               else if (ir.lt.imin) then ! Analytic Continuation: .lt.imin ==> in the region of no sampling. Set the FE (log(g)) to
                  ! the direct energy shifted by a constant energy term , the difference between the last sampled indirect and
                  ! direct energies.
                  ! ln(g(r<r0)) = -u_dir(r)/T - ( u_pmf(r0)/T - u_dir(r0)/T )
                  ! ln(g(r<r0)) = -u_dir(r)/T + ln(g(r0)) - u_dir(r0)/T
                  g(ia,ir,ith,iphi) = ( -idealHist(ia,1,ir,ith,iphi) + idealHist(ia,1,imin,ith,iphi) ) + g(ia,imin,ith,iphi) ! ln(g)
                  !print*, 'look2', idealHist(ia,1,30,1,1) !debug
                  if ((g(ia,ir,ith,iphi).lt.cut).and.(igo.eq.0)) then ! the largest r to go past the cutoff
                     ispline(ith,iphi) = ir
                     igo = 1
                  end if
                  fLJr(ia,ir,ith,iphi) = idealHist(ia, 2,ir,ith,iphi)
                  fLJs(ia,ir,ith,iphi) = idealHist(ia, 3,ir,ith,iphi)
                  fLJt(ia,ir,ith,iphi) = idealHist(ia, 4,ir,ith,iphi)
                  fCr(ia,ir,ith,iphi)  = idealHist(ia, 8,ir,ith,iphi)
                  fCs(ia,ir,ith,iphi)  = idealHist(ia, 9,ir,ith,iphi)
                  fCt(ia,ir,ith,iphi)  = idealHist(ia,10,ir,ith,iphi)
               end if
            end do ! ir
         end do ! iphi
      end do ! ith

      ! Set all ln(g) values past the largest r bin to reach the cuttoff to the cuttoff value.
      igo = 0
      do iphi = 1, histPhiBins
         do ith = 1, histCosThBins
            do ir = histDistBins, 1, -1
               if ((g(ia,ir,ith,iphi).lt.-abs(cut)).and.(igo.eq.0)) then
                  g(ia,1:ir,ith,iphi) = -abs(cut)
                  igo = 1
               end if
            end do
            igo = 0
         end do
      end do

      call reduced_mean_field(ia) ! this way the emf1D is made with the already analytically continued g(r,th,phi) which is ln(g)

      ! Set all LJ forces past (to the left of) the first (largest r) bin to reach the cutoff to the cutoff value.
      igo = 0; igo1 = 0; igo2 = 0
      do iphi = 1, histPhiBins
         do ith = 1, histCosThBins
            do ir = histDistBins, 1, -1
               if ((fLJr(ia,ir,ith,iphi).gt.abs(cut)).and.(igo.eq.0)) then
                  fLJr(ia,1:ir,ith,iphi) = abs(cut)
                  igo = 1
               end if
               if ((fLJs(ia,ir,ith,iphi).gt.abs(cut)).and.(igo1.eq.0)) then
                  fLJs(ia,1:ir,ith,iphi) = abs(cut)
                  igo1 = 1
               end if
               if ((fLJt(ia,ir,ith,iphi).gt.abs(cut)).and.(igo2.eq.0)) then
                  fLJt(ia,1:ir,ith,iphi) = abs(cut)
                  igo2 = 1
               end if
            end do
            igo = 0; igo1 = 0; igo2 = 0
         end do
      end do

      ! Set all C forces past (to the left of) the first (largest r) bin to reach the cutoff to the cutoff value.
      igo = 0; igo1 = 0; igo2 = 0
      do iphi = 1, histPhiBins
         do ith = 1, histCosThBins
            do ir = histDistBins, 1, -1
               if ((abs(fCr(ia,ir,ith,iphi)).gt.abs(cut)).and.(igo.eq.0)) then
                  fCr(ia,1:ir,ith,iphi) = sign(real(1,dp),idealHist(ia,8,ir,ith,iphi)) * abs(cut)
                  igo = 1
               end if
               if ((abs(fCs(ia,ir,ith,iphi)).gt.abs(cut)).and.(igo1.eq.0)) then
                  fCs(ia,1:ir,ith,iphi) = sign(real(1,dp),idealHist(ia,9,ir,ith,iphi)) * abs(cut)
                  igo1 = 1
               end if
               if ((abs(fCt(ia,ir,ith,iphi)).gt.abs(cut)).and.(igo2.eq.0)) then
                  fCt(ia,1:ir,ith,iphi) = sign(real(1,dp),idealHist(ia,10,ir,ith,iphi)) * abs(cut)
                  igo2 = 1
               end if
            end do
            igo = 0; igo1 = 0; igo2 = 0
         end do
      end do
!      do ir = 1, histDistBins !debug
!         print*, 'g3D: ', ir, g(ia,ir,1,1) !debug
!      end do !debug

      ! Average the 3D input histograms into 2D
      write(*,*) 'Averaging input histograms from 3D into 2D for solute: ',ia

      ! Find the minimum value of u(phi; r,th) ==> u02D(r,th)
      ! dim=3 in this case means the phi dimension. Replace the array in phi at each r,th with the minimum value of the array,
      ! making an array u02D(r,th).
      u02D = minval(-g(ia,:,:,:),dim=3)
!      do ir = 1, histDistBins !debug
!         print*, 'u0_2D: ',ir, u02D(ir,1) !debug
!      end do !debug

      do ith = 1, histCosThBins
         do ir = 1, histDistBins
            boltz_sum = 0_dp
            do iphi = 1, histPhiBins
               boltz = exp(g(ia,ir,ith,iphi) + u02D(ir,ith))
               !print*, 'boltz: ', g(ia,ir,ith,iphi), u02D(ir,ith), boltz !debug
               g2D(ia,ir,ith) = g2D(ia,ir,ith) + exp(g(ia,ir,ith,iphi)) ! g is g
               !print*, 'g2Dloop: ', ith,ir,iphi, g2D(ia,ir,ith) !debug
               fLJr2D(ia,ir,ith) = fLJr2D(ia,ir,ith) + (boltz * fLJr(ia,ir,ith,iphi)) ! fLJ.r
               fLJs2D(ia,ir,ith) = fLJs2D(ia,ir,ith) + (boltz * fLJs(ia,ir,ith,iphi)) ! fLJ.s
               fCr2D(ia,ir,ith) = fCr2D(ia,ir,ith) + (boltz * fCr(ia,ir,ith,iphi)) ! fC.r
               fCs2D(ia,ir,ith) = fCs2D(ia,ir,ith) + (boltz * fCs(ia,ir,ith,iphi)) ! fC.s
               boltz_sum = boltz_sum + boltz ! denominator for averaging over phi
            end do
            !print*, 'boltz_sum: ', boltz_sum !debug
            g2D(ia,ir,ith) = log(g2D(ia,ir,ith) / real(histPhiBins,dp)) ! finish average over phi by dividing and converting to log(g)
            if (g2D(ia,ir,ith).lt.-abs(cut)) g2D(ia,ir,ith) = -abs(cut)
            fLJr2D(ia,ir,ith) = fLJr2D(ia,ir,ith) / boltz_sum
            fLJs2D(ia,ir,ith) = fLJs2D(ia,ir,ith) / boltz_sum
            fCr2D(ia,ir,ith) = fCr2D(ia,ir,ith) / boltz_sum
            fCs2D(ia,ir,ith) = fCs2D(ia,ir,ith) / boltz_sum
         end do
         !print*, ith, fLJr2D(ia,1:40,ith) !debug
      end do
!      do ir = 1, histDistBins !debug
!         print*, 'g2D: ', ir, g2D(ia,ir,1) !debug
!      end do !debug

      ! Average the 2D input histograms into 1D
      write(*,*) 'Averaging input histograms from 2D into 1D for solute: ',ia

      u01D(:) = minval(-g2D(ia,:,:),dim=2)
!      do ir = 1, histDistBins !debug
!         print*, 'u0_1D: ',ir, u01D(ir) !debug
!      end do !debug

      do ir = 1, histDistBins
         boltz_sum = 0_dp
         do ith = 1, histCosThBins
            boltz = exp(g2D(ia,ir,ith) + u01D(ir))
            !print*, 'boltz: ', g2D(ia,ir,ith), u01D(ir), boltz !debug
            g1D(ia,ir) = g1D(ia,ir) + exp(g2D(ia,ir,ith)) ! g1D is g
            fLJr1D(ia,ir) = fLJr1D(ia,ir) + (boltz * fLJr2D(ia,ir,ith)) ! fLJ.r
            fCr1D(ia,ir) = fCr1D(ia,ir) + (boltz * fCr2D(ia,ir,ith)) ! fC.r
            boltz_sum = boltz_sum + boltz ! denominator for averaging over theta
         end do
         !print*, 'boltz_sum: ', boltz_sum !debug
         g1D(ia,ir) = log(g1D(ia,ir) / real(histCosThBins,dp)) ! finish average over cosTh by dividing and converting to log(g)
         fLJr1D(ia,ir) = fLJr1D(ia,ir) / boltz_sum
         fCr1D(ia,ir) = fCr1D(ia,ir) / boltz_sum
      end do
      !print*, 'fLJ1D: ', fLJr1D(ia,1:40) !debug
      ! After the averaging is done enforce the cutoff
      do ir = 1, histDistBins
         if (g1D(ia,ir).lt.cut) then
            g1D(ia,ir) = cut
            ispline1D = ir
         end if
         if (fLJr1D(ia,ir).gt.abs(cut)) then
            fLJr1D(ia,ir) = -cut
         end if
         if (abs(fCr1D(ia,ir)).gt.abs(cut)) then
            fCr1D(ia,ir) = sign(real(1,dp),idealHist1D(ia,4,ir))*cut
         end if
      end do
      ! spline g1D for the solute
      call spline(histDist,g1D(ia,:),ispline1D,histDistBins,(idealHist1D(ia,2,ispline1D)+idealHist1D(ia,4,ispline1D)),real(0,dp), &
         & g1D2(ia,:))
      call spline(histDist,fLJr1D(ia,:),ispline1D,histDistBins,idealHist1D(ia,3,ispline1D),real(0,dp), fLJr1D2(ia,:))
      call spline(histDist,fCr1D(ia,:),ispline1D,histDistBins,idealHist1D(ia,5,ispline1D),real(0,dp), fCr1D2(ia,:))
      call spline(histDist,emf1D(ia,:),ispline1D,histDistBins,real(0,dp),real(0,dp), emf1D2(ia,:))
   end do ! ia

   ! note: write out the effective input histogram after averaging/alterations.
   write(*,*) 'Writing input histogram after averaging/alterations to "input_hist.out" ...'
   open(91,file='input_hist.out',status='replace')
   write(91,*) '#  1. Distance'
   write(91,*) '#  2. g+'
   write(91,*) '#  3. fLJ.r+'
   write(91,*) '#  4. fC.r+'
   write(91,*) '#  5. emf.r+'
   write(91,*) '#  6. g-'
   write(91,*) '#  7. fLJ.r-'
   write(91,*) '#  8. fC.r-'
   write(91,*) '#  9. emf.r-'
   do ir = 1, histDistBins
      write(91,*) histDist(ir), g1D(1,ir), fLJr1D(1,ir), fCr1D(1,ir), emf1D(1,ir), g1D(2,ir), fLJr1D(2,ir), fCr1D(2,ir), emf1D(2,ir)
   end do
   close(91)

   ! note: write out the ideal histogram after averaging/alterations.
   write(*,*) 'Writing ideal histogram after averaging/alterations to "ideal_hist.out" ...'
   open(92,file='ideal_hist.out',status='replace')
   write(92,*) '# 1. Distance'
   write(92,*) '# 2. g+'
   write(92,*) '# 3. fLJ.r+'
   write(92,*) '# 4. fC.r+'
   write(92,*) '# 5. g-'
   write(92,*) '# 6. fLJ.r-'
   write(92,*) '# 7. fC.r-'
   do ir = 1, histDistBins
      write(92,*) histDist(ir), -idealHist1D(1,1,ir), idealHist1D(1,2,ir), idealHist1D(1,4,ir), -idealHist1D(2,1,ir), &
         & idealHist1D(2,2,ir), idealHist1D(2,4,ir)
   end do
   close(92)

   open(93,file='3dhist.out',status='replace')
   write(93,*) '#  1. Distance'
   write(93,*) '#  2. cosTh'
   write(93,*) '#  3. phi'
   write(93,*) '#  4. g+'
   write(93,*) '#  5. g-'
   do ir = 1, histDistBins
      do ith = 1, histCosThBins
         do iphi = 1, histPhiBins
            write(93,*) histDist(ir), histCosTh(ith), histPhi(iphi), g(1,ir,ith,iphi), g(2,ir,ith,iphi)
         end do
      end do
   end do
   close(93)
   !print*, 'look4', g(1,30,1,1), g(2,30,1,1)

   !debug
   write(*,*) 'Writing ideal histogram after averaging/alterations to "ideal_3D.out" ...'
   open(94,file='ideal_3D.out',status='replace')
   write(94,*) '# 1. Distance'
   write(94,*) '# 2. cosTheta'
   write(94,*) '# 3. Phi'
   write(94,*) '# 4. g+'
   write(94,*) '# 5. fLJ.r+'
   write(94,*) '# 6. fC.r+'
   write(94,*) '# 7. g-'
   write(94,*) '# 8. fLJ.r-'
   write(94,*) '# 9. fC.r-'
   do ir = 1, histDistBins
      do ith = 1, histCosThBins
         do iphi = 1, histPhiBins
            write(94,*) histDist(ir), histCosTh(ith), histPhi(iphi), -idealHist(1,1,ir,ith,iphi), idealHist(1,2,ir,ith,iphi), &
               & idealHist(1,8,ir,ith,iphi), -idealHist(2,1,ir,ith,iphi), idealHist(2,2,ir,ith,iphi), idealHist(2,8,ir,ith,iphi)
         end do
      end do
   end do
   close(94)

   !debug   difference between ideal and measured 
   !rTmp = int(real(4,dp)/histDistStepSize)
   !print*, rTmp, histDistStepSize
   !do ith=1,histCosThBins
      !write(55,*) histCosTh(ith), g2D(rTmp,ith), idealhist2D(1,rTmp,ith)
   !end do

   !debug
   !do ir=1,histdistbins
      !write(45,*) histDist(ir), g1D(ir), g1D2(ir), idealHist1D(1,ir), fLJr1D(ir), idealHist1D(2,ir)
   !end do
   !do ir=1,100*histdistbins
      !xx = ir*(histDistStepSize/100_dp)
      !call splint(histDist,g1D,g1D2,histDistBins,xx, yy)
      !write(55,*) xx, yy
   !end do
   !write(65,*) histDist(ispline1D), g1D(ispline1D)

end subroutine spline_hist_array


! read LJ--LJ displacements from file
subroutine R_list
   use cfgData; use ctrlData
   implicit none
   integer :: ios, i
   character(len=16) :: junk
   character(len=128) :: line

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


! setup for the average force integral
subroutine setup_compute_avg_force
   use cfgData; use angleData; use ctrlData
   implicit none
   integer :: i
   real(kind=dp) :: psiLF

   write(*,*) "Setting up for average force iteration..."

   if (explicit_R .eqv. .true.) then
      cfgRBins = crdLines
      write(*,*) "Number of R Bins:      ", cfgRBins
   else if (explicit_R .eqv. .false.) then
      cfgRBins = int( (R_max - R_min)/RStepSize + 1 )
      if (cfgRBins .eq. 0) then
         cfgRBins = 1
      end if
      write(*,*) "Number of R Bins:      ", cfgRBins
   end if
   xBins = int( (2 * xz_range)/xzStepSize )
   write(*,*) "Number of X Bins:      ", xBins
   zBins = int( (xz_range)/xzStepSize )
   write(*,*) "Number of Z Bins:      ", zBins

   ! allocate array sizes for axes and average force
   allocate( R_axis(cfgRBins), fAvg(2,cfgRBins), x_axis(xBins), z_axis(zBins) )
   R_axis = 0_dp; fAvg = 0_dp; x_axis = 0_dp; z_axis = 0_dp

   ! allocate arrays for control arrays
   ! frcSPA(solutes, LJ/Coulomb, x, z)
   allocate( frcSPA(2, 2, xBins, zBins), grSPA(xBins, zBins) )

   ! Distance Axes
   do i = 1, cfgRBins
      if (explicit_R .eqv. .true.) then
         R_axis(i) = explicitDist(i)
      else if (explicit_r .eqv. .false.) then
         R_axis(i) = (i-1) * RStepSize + R_min
      end if
   end do
   do i = 1, xBins
      x_axis(i) = (i-1) * xzStepSize - xz_range + xzStepSize/2_dp
   end do
   do i = 1, zBins
      z_axis(i) = (i-1) * xzStepSize + xzStepSize/2_dp
   end do

   ! ANGLES
   allocate( cosThetaLF(cfgCosThBins), sinThetaLF(cfgCosThBins), sinPsiLF(cfgPsiBins), cosPsiLF(cfgPsiBins) )

   ! Theta
   ! tilt off of z
   cfgCosThStepSize = (cosTh_max - cosTh_min) / real(cfgCosThBins, dp)
   do i = 1, cfgCosThBins
      cosThetaLF(i) = (i-0.5_dp)*cfgCosThStepSize - cosTh_max
      sinThetaLF(i) = sqrt(abs(1_dp-cosThetaLF(i)**2))
   end do
   write(*,*) "Config Cos(Theta) Step Size:      ", cfgCosThStepSize

   ! Psi
   ! processison about z
   cfgPsiStepSize = (psi_max - psi_min) / real(cfgPsiBins, dp)
   do i = 1, cfgPsiBins
      psiLF = (i+0.5_dp)*cfgPsiStepSize
      sinPsiLF(i) = sin(psiLF)
      cosPsiLF(i) = cos(psiLF)
   end do
   write(*,*) "Config Psi Step Size:         ", cfgPsiStepSize

end subroutine setup_compute_avg_force


! do the average force integral
subroutine compute_avg_force
   use cfgData; use histData; use angleData; use ctrlData; use constants; use functions
   implicit none
   integer :: r, i, j, ip, omp_get_thread_num !, omp_get_num_threads
   real(kind=dp) :: gx1, gx2, flj, emfVec(3), emf, fx1, fx2, polMeanVec(3), emfVecMag

   write(*,*) "Computing average force..."; flush(6)

   allocate(longRange(cfgRBins)) !debug

   ! Calculate the average force integral for top half of bisecting plane of cylinder
   do r = 1, cfgRBins ! loop lj--lj distances 
!   do r = cfgRBins, cfgRBins ! loop lj--lj distances !debug
      frcSPA = 0_dp; grSPA = 0_dp
      !$omp PARALLEL DEFAULT( none ) &
      !$omp PRIVATE( ip, i, j, gx1, gx2, fx1, fx2, flj, emfVec, emf, polMeanVec, emfVecMag ) &
      !$omp SHARED( nThreads, r, xBins, zBins, cut, R_axis, x_axis, z_axis, histDist, histDistBins, g1D, g1D2, fLJr1D, fLJr1D2, &
      !$omp&   fCr1D, fCr1D2, emf1D, emf1D2, frcSPA, grSPA, soluteChg, dielectric, density, dipole_moment ) &
      !$omp NUM_THREADS( nThreads )
      if ((omp_get_thread_num().eq.0).and.(r.eq.1)) then
         write(*,*) 'Parallel CPUs:         ', nThreads
         !write(*,*) 'Parallel CPUs:         ', omp_get_num_threads()
         flush(6)
      end if
      !$omp DO SCHEDULE( guided )
      do ip = 1, (xBins*zBins)
         ! Convert single index 'ip' to the x and z indicies 'i' and 'j' respectively.
         i = int((ip-1)/zBins)+1 ! x integer
         j = mod(ip-1,zBins)+1   ! z integer

         rSolv1(1) = -R_axis(r)/2_dp - x_axis(i)
         rSolv1(2) = 0_dp
         rSolv1(3) = -z_axis(j)
         rSolvn(1) = euclid_norm(rSolv1)

         rSolv2(1) = R_axis(r)/2_dp - x_axis(i)
         rSolv2(2) = 0_dp
         rSolv2(3) = -z_axis(j)
         rSolvn(2) = euclid_norm(rSolv2)
         !if ((i.eq.350).and.(j.eq.200)) print*, 'x,z,r1,r2', x_axis(i), z_axis(j), rSolvn !debug 1
         !if ((i.eq.350).and.(j.eq.200)) print*, 'r1,r2', rSolv1, rSolv2 !debug 1

         if (rSolvn(1).gt.histDist(histDistBins)) then
            gx1 = 0_dp
         else
            call splint(histDist,g1D(1,:),g1D2(1,:),histDistBins,rSolvn(1), gx1)   ! solute 1 density
         end if

         if (rSolvn(2).gt.histDist(histDistBins)) then
            gx2 = 0_dp
         else
            call splint(histDist,g1D(2,:),g1D2(2,:),histDistBins,rSolvn(2), gx2)   ! solute 2 density
         end if
         gx1 = exp(gx1+gx2)

         if ((gx1.gt.1d-6).and.((rSolvn(1).le.histDist(histDistBins)).or.(rSolvn(2).le.histDist(histDistBins)))) then ! inside vol
            if (rSolvn(1).gt.histDist(histDistBins)) then
               ! outside solute 1 interaction volume use ideal field in CDD
               !emfVec = (coulomb_const_kcalmolAq2/dielectric*soluteChg(1)/rSolvn(1)**2) * (-rSolv1/rSolvn(1))
               emfVec = (3.0/4.0/pi/density/dipole_moment*soluteChg(1)*(1.-1./dielectric)/rSolvn(1)**2) * (-rSolv1/rSolvn(1))
               fx1 = 0_dp ! LJ force is zero
            else
               call splint(histDist,fLJr1D(1,:),fLJr1D2(1,:),histDistBins,rSolvn(1), flj) ! lj force
               call splint(histDist,emf1D(1,:),emf1D2(1,:),histDistBins,rSolvn(1), emf) ! field
               emfVec = emf * rSolv1/rSolvn(1)
               !print*, 'emf:', emf !debug
               fx1 = flj
               !if ((i.eq.350).and.(j.eq.200)) print*, 'emf1', emf !debug 1
            end if

            if (rSolvn(2).gt.histDist(histDistBins)) then ! solute 2
               ! outside solute 2 interaction volume use ideal field in CDD
               !emfVec = emfVec + (coulomb_const_kcalmolAq2/dielectric*soluteChg(2)/rSolvn(2)**2) * (-rSolv2/rSolvn(2))
               emfVec = (3.0/4.0/pi/density/dipole_moment*soluteChg(2)*(1.-1./dielectric)/rSolvn(2)**2) * (-rSolv2/rSolvn(2))
               fx2 = 0_dp ! LJ force is zero
            else
               call splint(histDist,fLJr1D(2,:),fLJr1D2(2,:),histDistBins,rSolvn(2), flj)
               call splint(histDist,emf1D(2,:),emf1D2(2,:),histDistBins,rSolvn(2), emf)
               emfVec = emfVec + emf * rSolv2/rSolvn(2)
               fx2 = flj
               !if ((i.eq.350).and.(j.eq.200)) print*, 'emf2', emf !debug 1
               !if ((i.eq.350).and.(j.eq.200)) print*, 'emf_sum', emfVec !debug 1
            end if

            emfVecMag = euclid_norm(emfVec)
            polMeanVec = (tanh(emfVecMag)**(-1) - (emfVecMag)**(-1))*emfVec/emfVecMag ! Langevin fxn Emf --> <p> @ cell
            ! Once you have the polariztion via Langevin, calculate force
            emfVec = soluteChg(1)*dipole_moment*coulomb_const_kcalmolAq2/rSolvn(1)**3 * &
               & (3*dot_product(rSolv1/rSolvn(1),polMeanVec) * rSolv1/rSolvn(1) - polMeanVec) ! dipole force solute 1
            frcSPA(1,2,i,j) = frcSPA(1,2,i,j) - (gx1 * emfVec(1)) ! fd*g.R^hat coulomb solute 1
            !if ((i.eq.350).and.(j.eq.200)) print*, 'frc1', euclid_norm(emfVec) !debug 1

            emfVec = soluteChg(2)*dipole_moment*coulomb_const_kcalmolAq2/rSolvn(2)**3 * &
               & (3*dot_product(rSolv2/rSolvn(2),polMeanVec) * rSolv2/rSolvn(2) - polMeanVec) ! dipole force solute 2
            frcSPA(2,2,i,j) = frcSPA(2,2,i,j) + (gx1 * emfVec(1)) ! fd*g.R^hat coulomb solute 2
            !if ((i.eq.350).and.(j.eq.200)) print*, 'frc2', euclid_norm(emfVec) !debug 1
             
            frcSPA(1,1,i,j) = frcSPA(1,1,i,j) + (gx1 * fx1 * (-rSolv1(1)/rSolvn(1))) ! (f.r)*g.R^{hat} lj solute 1
            frcSPA(2,1,i,j) = frcSPA(2,1,i,j) + (gx1 * fx2 * ( rSolv2(1)/rSolvn(2))) ! (f.r)*g.R^{hat} lj solute 2
            grSPA(i,j) = grSPA(i,j) + gx1 ! gx1 is the SPA at this point
         end if
      end do !ip
      !$omp END DO
      !$omp END PARALLEL
       
      ! NOTE: Long range correction to mean force for Constant Density Dielectric medium at large distances
      emfVecMag = -soluteChg(1)*soluteChg(2)*coulomb_const_kcalmolAq2*(1-dielectric**(-1)) * &
         & ( (R_axis(r)*(8*histDist(histDistBins)-3*R_axis(r)))/(24*histDist(histDistBins)**4) &
         & - log(1+R_axis(r)/histDist(histDistBins))/(8*R_axis(r)**2) &
         & + histDist(histDistBins)/(8*R_axis(r)*(R_axis(r)+histDist(histDistBins))**2) &
         & - R_axis(r)**2/(32*histDist(histDistBins)**4) &
         & + 3/(16*histDist(histDistBins)**2) )
       
      ! Add each cell forces to average and normalize
      do i = 1, xBins
         do j = 1, zBins
            fAvg(1,r) = fAvg(1,r) + (frcSPA(1,1,i,j) + frcSPA(2,1,i,j)) * real(0.5,dp) * z_axis(j) ! lj
            fAvg(2,r) = fAvg(2,r) + (frcSPA(1,2,i,j) + frcSPA(2,2,i,j)) * real(0.5,dp) * z_axis(j) ! Coulomb
            frcSPA(1,1,i,j) = frcSPA(1,1,i,j)/grSPA(i,j) ! lj solute 1
            frcSPA(2,1,i,j) = frcSPA(2,1,i,j)/grSPA(i,j) ! lj solute 1
            frcSPA(1,2,i,j) = frcSPA(1,2,i,j)/grSPA(i,j) ! coulomb solute 2
            frcSPA(2,2,i,j) = frcSPA(2,2,i,j)/grSPA(i,j) ! coulomb solute 2
         end do !z again
      end do !x again
      call write_test_out(r) ! write grSPA and frcSPA arrays
       
      ! NOTE : After the fact multiply all elements by 2*pi*density/8/pi/pi (2*2pi*pi/3 (4pi**2)/3 steradians from orientations)
      !       Number density of chloroform per Angstrom**3 == 0.00750924
      fAvg(1,r) = fAvg(1,r)*2*pi*density*xzStepSize*xzStepSize
      fAvg(2,r) = fAvg(2,r)*2*pi*density*xzStepSize*xzStepSize !debug+ emfVecMag ! with long range correction added
      longRange(r) = emfVecMag !debug
   end do !r

   !debug
   open(95,file='longRange.out',status='replace')
   do r = 1, cfgRBins
      write(95,*) R_axis(r), longRange(r)
   end do
   close(95)

end subroutine compute_avg_force


! integrate the average force from 'compute_avg_force' to get the PMF.
subroutine integrate_force
   use cfgData; use ctrlData
   implicit none
   integer :: d, f

   allocate( u_dir(2,cfgRBins) )
   u_dir = 0_dp

   do f = 1, 2
      if (explicit_R .eqv. .false.) then
         do d = 1, cfgRBins
            if (d .eq. 1) then
               u_dir(f,cfgRBins) = fAvg(f,cfgRBins) * RStepSize
            else
               u_dir(f,cfgRBins-(d-1)) = u_dir(f,cfgRBins-(d-2)) + fAvg(f,cfgRBins-(d-1)) * RStepSize
            end if
         end do
      else if (explicit_R .eqv. .true.) then
         do d = 1, cfgRBins
            if (d .eq. 1) then
               u_dir(f,cfgRBins) = fAvg(f,cfgRBins) * (R_axis(cfgRBins)-R_axis(cfgRBins-1))
               !print*, (R_axis(cfgRBins)-R_axis(cfgRBins-1))
            else
               ! FIXME: is the delta R part of this correct?
               u_dir(f,cfgRBins-(d-1)) = u_dir(f,cfgRBins-(d-2)) + fAvg(f,cfgRBins-(d-1)) * &
                  (R_axis(cfgRBins-(d-1))-R_axis(cfgRBins-d))
               ! it looks like the first value is getting printed twice. Also, the values might be wrong. Should it be (d-1) and
               ! (d-0)? instead of -2 and -1?
               !print*, (R_axis(cfgRBins-(d-1))-R_axis(cfgRBins-d))
            end if
         end do
      end if
   end do
end subroutine integrate_force


! write force out and g(r) out to compare against explicit
subroutine write_test_out(r)
   use cfgData; use ctrlData
   implicit none
   integer :: r, i_f, i, j
   character(len=32) :: temp, filename
   character(len=8) :: frmt

   i_f = (r-1) * int(RStepSize*10)
   frmt = '(I3.3)' ! an integer of width 3 with zeroes on the left
   write(temp,frmt) i_f ! converting integer to string using 'internal file'
   filename='hist1D_output.'//trim(temp)//'.dat'


   open(35,file=filename,status='replace')
   write(6,*) "Writing test file:   ", filename
   write(35,*) "# 1.   X Distance"
   write(35,*) "# 2.   Z Distance"
   write(35,*) "# 3.   g(r)"
   write(35,*) "# 4.   Force.r + LJ"
   write(35,*) "# 5.   Force.r + Coulomb"
   write(35,*) "# 6.   Force.r - LJ"
   write(35,*) "# 7.   Force.r - Coulomb"
   write(35,*) "# "
   do j = 1, zBins
      do i = 1, xBins
         write(35,898) x_axis(i), z_axis(j), grSPA(i,j), frcSPA(1,1,i,j), frcSPA(1,2,i,j), frcSPA(2,1,i,j), frcSPA(2,2,i,j)
      end do
   end do
   close(35)

   flush(6)
   
898      format (2(1x,es14.7),5(1x,es14.7))
end subroutine write_test_out


! write output file
subroutine write_output(outFile)
   use cfgData
   implicit none
   character(len=128) :: outFile
   integer :: r

   open(35,file=outFile,status='replace')
   write(6,*) "Writing output file:   ", outFile
   write(35,*) "# 1.   R Distance"
   write(35,*) "# 2.   <f>_LJ"
   write(35,*) "# 3.   <f>_Coulomb"
   write(35,*) "# 4.   PMF LJ"
   write(35,*) "# 5.   PMF Coulomb"
   do r = 1, cfgRBins
      write(35,899) R_axis(r), fAvg(1,r), fAvg(2,r), u_dir(1,r), u_dir(2,r)
   end do
   close(35)

   flush(6)

899      format (5(1x,es14.7)) ! scientific format
end subroutine write_output
