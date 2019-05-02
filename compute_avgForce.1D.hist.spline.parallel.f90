! USAGE: ./this_file.x -cfg [CONFIGURATION FILE] -hist [3D HIST FILE]
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
   real(kind=dp),allocatable   :: histDist(:), histCosTh(:), histPhi(:), g(:,:,:), fr(:,:,:), fs(:,:,:), ft(:,:,:), g2(:,:), &
      & fr2(:,:), fs2(:,:), gc(:,:,:), gTmp1(:,:), gTmp2(:,:), frTmp1(:,:), fsTmp1(:,:)
   real(kind=dp),allocatable   :: g2D(:,:), fr2D(:,:), fs2D(:,:), g2D2(:,:), fr2D2(:,:), fs2D2(:,:)
   real(kind=dp),allocatable   :: g1D(:), fr1D(:), g1D2(:), fr1D2(:)
   real(kind=dp)   :: histDistStepSize, histCosThStepSize, histPhiStepSize
   integer   :: histDistBins, histCosThBins, histPhiBins

   !$omp THREADPRIVATE( gTmp1, gTmp2, frTmp1, fsTmp1 )
end module histData

! data from the config file.
module cfgData
   use prec
   use constants
   real(kind=dp),allocatable :: x_axis(:), z_axis(:), R_axis(:), fAvg(:), u_dir(:)
   real(kind=dp) :: RStepSize, xzStepSize, R_min, R_max, xz_range, cfgCosThStepSize, cfgPsiStepSize, T, cut, offset
   character(len=8) :: c_explicit_R
   integer :: cfgRBins, cfgCosThBins, cfgPhiBins, cfgPsiBins, radius
!
   integer   :: xBins, zBins
   real(kind=dp)   :: density = 0.00750924_dp ! numerical density of chloroforms per Angstrom**3
   real(kind=dp)   :: cosTh_max = 1_dp
   real(kind=dp)   :: cosTh_min = -1_dp
   real(kind=dp)   :: phi_max = pi/3_dp
   real(kind=dp)   :: phi_min = 0_dp
   real(kind=dp)   :: psi_max = 2_dp*pi
   real(kind=dp)   :: psi_min = 0_dp
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
   real(kind=dp),allocatable   :: frcSPA(:,:,:), grSPA(:,:), explicitDist(:)
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
   real(kind=dp)       :: omp_get_wtime, ti, tf, seconds
   integer            :: hours, minutes

   ti = omp_get_wtime()

   ! make list of average direct force from 'collapsed' file.
   call parse_command_line(histFile, cfgFile)

   ! read config file
   call read_cfg(cfgFile, outFile)

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
subroutine parse_command_line(histFile, cfgFile) !,outFile)
   implicit none
   character(len=64)   :: histFile, cfgFile !, outFile
   character(len=16)   :: arg
   integer            :: i
   logical            :: histFileFlag, cfgFileFlag, histExist, cfgExist

   histFileFlag = .false.
   cfgFileFlag = .false.
   histExist = .false.
   cfgExist = .false.
   i=1
   do
      call get_command_argument(i,arg)
      select case (arg)

      case ('-hist')
         i = i+1
         call get_command_argument(i,histFile)
         histFileFlag=.true.
         INQUIRE(FILE=histFile, EXIST=histExist)
         write(*,*) 'Histogram File:         ', histFile
         write(*,*) 'Histogram File Exists:         ', histExist
      case ('-cfg')
         i = i+1
         call get_command_argument(i,cfgFile)
         cfgFileFlag=.true.
         INQUIRE(FILE=cfgFile, EXIST=cfgExist)
         write(*,*) 'Config File:            ', cfgFile
         write(*,*) 'Config File Exists:         ', cfgExist
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

   ! 'ERROR STOP' if either file doesn't exist
   if ((histExist.eqv..false.).or.(cfgExist.eqv..false.)) then
      write(*,*) 'hist or cfg files do not exist'
      error stop
   end if

   flush(6)

end subroutine parse_command_line


! read python cfg file for g(r) parameters
subroutine read_cfg(cfgFile, outFile)
   use cfgData
   implicit none
   character(len=64)    :: cfgFile, outFile
   character(len=128)   :: line
   character(len=32)    :: firstWord, sep
   integer              :: ios
   logical              :: outFileFlag, RstepSizeFlag, xzStepSizeFlag, RmaxFlag, RminFlag, xzRangeFlag, thetaBinsFlag, phiBinsFlag,&
      & psiBinsFlag, c_explicit_RFlag, TFlag, cutFlag, radiusFlag, offsetFlag

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

   ios = 0

   open(20,file=cfgFile)
   do while(ios>=0)
      read(20,'(a)',IOSTAT=ios) line
      call split(line,'=',firstWord, sep)
      if (line .ne. "") then
         if (firstWord .eq. "out_file") then
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

   flush(6)

end subroutine read_cfg


! read force file and make a lookup table.
subroutine make_hist_table(histFile)
   use histData
   use cfgData
   implicit none
   character(len=64)            :: histFile
   character(len=64)            :: junk
   character(len=512)         :: line
   integer                     :: ios, ios2, i, j, k, nHistLines
   real(kind=dp),allocatable   :: histTmp(:,:)

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
         !                  dist         cos(Th)         phi/3         g(r)+         g(r)-      <f.r>+         <f.s>+      <f.t>+
         read(line,*) histTmp(1,i), histTmp(2,i), histTmp(3,i), histTmp(4,i), junk, histTmp(5,i), histTmp(6,i), histTmp(7,i), &
            !      gc(r)+    gc(r)-
            & histTmp(8,i), junk
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
   
   allocate( histDist(histDistBins), histCosTh(histCosThBins), histPhi(histPhiBins), g(histDistBins,histCosThBins,histPhiBins), &
      & fr(histDistBins,histCosThBins,histPhiBins), fs(histDistBins,histCosThBins,histPhiBins), &
      & ft(histDistBins,histCosThBins,histPhiBins), gc(histDistBins,histCosThBins,histPhiBins) )
   
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
            g(i,j,k) = histTmp(4, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)   ! g(r,cos,phi) currently g
            fr(i,j,k) = histTmp(5, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <f.r>(r,cos,phi)
            fs(i,j,k) = histTmp(6, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <f.s>(r,cos,phi)
            ft(i,j,k) = histTmp(7, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! <f.t>(r,cos,phi)
            gc(i,j,k) = histTmp(8, (i-1)*histCosThBins*histPhiBins + (j-1)*histPhiBins + k)  ! gc(r,cos,phi)
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

end subroutine make_hist_table


! spline the r dimension of each theta phi stack and then average over phi for 2D
subroutine spline_hist_array
   use constants
   use functions
   use histData
   use cfgData
   use idealSolv
   implicit none
   integer                     :: i, ir, ith, iphi, imin, igo, ir2, igo1, igo2
   real(kind=dp)               :: x, y, norm_factor, boltz, boltz_sum
   real(kind=dp),allocatable   :: idealHist(:,:,:,:), idealHist2D(:,:,:), u0(:,:), idealHist1D(:,:), u01D(:)
   integer,allocatable         :: ispline(:,:), ispline2D(:), ispline1D
   !real(kind=dp)   :: xx, yy   !debug
   integer  :: rTmp !debug

   write(*,*) 'Editing input histogram arrays with ideal arrays in 3D...'

   ! Calculate a 4D array idealHist(g/f,r,th,phi)
   allocate( idealHist(7,histDistBins,histCosThBins,histPhiBins), ispline(histCosThBins,histPhiBins), ispline2D(histCosThBins) )
   idealHist = 0_dp; ispline = 0_dp
   call ideal_CL3(histDistBins,histDistStepSize,histCosThBins,cosTh_min,cosTh_max,histPhiBins,phi_min,phi_max,radius,offset, &
      & T, idealHist)

   ! Edit the input hist arrays to more smoothly transition to -/+ infinity with the help of idealHist.
   do ith = 1, histCosThBins
      do iphi = 1, histPhiBins
         imin = 0
         ! Normalization factor or each theta phi array.
         norm_factor = gc(histDistBins,ith,iphi)/(g(histDistBins,ith,iphi)*4*pi*histDist(histDistBins)**2)
         ! Find the first non-zero g(r) bin for each theta/phi array and set 'imin' to that 'ir' index
find:      do ir = 1, histDistBins
            if (g(ir,ith,iphi).gt.1d-6) then
               imin = ir
               exit find
            end if
         end do find
          
         ! Note: Add the ideal values to bins with no sampling. And half counts to bins that probably should have had sampling.
         igo = 0
         do ir = histDistBins, 1, -1
            if (ir.ge.imin) then
               if (g(ir,ith,iphi).gt.1d-6) then
                  g(ir,ith,iphi) = log(g(ir,ith,iphi))
               else   ! note: This is a zero bin where there probably should have been something. So put a half count in.
                  g(ir,ith,iphi) = log(real(0.5,dp)/(4*pi*histDist(ir)**2)/norm_factor)
                  fr(ir,ith,iphi) = idealHist(2,ir,ith,iphi)
                  fs(ir,ith,iphi) = idealHist(3,ir,ith,iphi)
                  ft(ir,ith,iphi) = idealHist(4,ir,ith,iphi)
               end if
            else if (ir.lt.imin) then ! .lt.imin ==> in the region of no sampling. Set the FE (log(g)) to the direct energy shifted by
               ! a constant energy term , which is the last sampled indirect energy.
               ! ln(g(r<r0)) = -u_dir(r)/T - ( u_pmf(r0)/T - u_dir(r0)/T )
               ! ln(g(r<r0)) = -u_dir(r)/T + ln(g(r0)) - u_dir(r0)/T
               g(ir,ith,iphi) = ( -idealHist(1,ir,ith,iphi) + idealHist(1,imin,ith,iphi) ) + g(imin,ith,iphi)
               if ((g(ir,ith,iphi).lt.cut).and.(igo.eq.0)) then ! the largest r to go past the cutoff
                  ispline(ith,iphi) = ir
                  igo = 1
               end if
               fr(ir,ith,iphi) = idealHist(2,ir,ith,iphi)
               fs(ir,ith,iphi) = idealHist(3,ir,ith,iphi)
               ft(ir,ith,iphi) = idealHist(4,ir,ith,iphi)
            end if
         end do
      end do
   end do

   ! Set all forces past the first (largest r) bin to reach the cutoff to the cutoff value.
   igo = 0; igo1 = 0; igo2 = 0
   do iphi = 1, histPhiBins
      do ith = 1, histCosThBins
         do ir = histDistBins, 1, -1
            if ((fr(ir,ith,iphi).gt.-cut).and.(igo.eq.0)) then
               fr(1:ir,ith,iphi) = -cut
               igo = 1
            end if
            if ((fs(ir,ith,iphi).gt.-cut).and.(igo1.eq.0)) then
               fs(1:ir,ith,iphi) = -cut
               igo1 = 1
            end if
            if ((ft(ir,ith,iphi).gt.-cut).and.(igo2.eq.0)) then
               ft(1:ir,ith,iphi) = -cut
               igo2 = 1
            end if
         end do
         igo = 0; igo1 = 0; igo2 = 0
      end do
   end do

   ! Average the 3D input histograms into 2D
   write(*,*) 'Averaging input histograms from 3D into 2D...'

   allocate( g2D(histDistBins,histCosThBins), fr2D(histDistBins,histCosThBins), fs2D(histDistBins,histCosThBins), &
      & u0(histDistBins,histCosThBins) )
   g2D = 0_dp; fr2D = 0_dp; fs2D = 0_dp; u0 = 0_dp

   ! Find the minimum value of u(phi; r,th) ==> u0(r,th)
   ! dim=3 in this case means the phi dimension. Replace the array in phi at each r,th with the minimum value of the array, making
   ! an array u0(r,th).
   u0 = minval(-g(:,:,:),dim=3)

   do ith = 1, histCosThBins
      do ir = 1, histDistBins
         boltz_sum = 0_dp
         do iphi = 1, histPhiBins
            boltz = exp(g(ir,ith,iphi) + u0(ir,ith))
            g2D(ir,ith) = g2D(ir,ith) + exp(g(ir,ith,iphi)) ! g is g
            fr2D(ir,ith) = fr2D(ir,ith) + (boltz * fr(ir,ith,iphi)) ! f.r
            fs2D(ir,ith) = fs2D(ir,ith) + (boltz * fs(ir,ith,iphi)) ! f.s
            boltz_sum = boltz_sum + boltz ! denominator for averaging over phi
         end do
         g2D(ir,ith) = log(g2D(ir,ith) / real(histPhiBins,dp)) ! finish average over phi by dividing and converting to log(g)
         fr2D(ir,ith) = fr2D(ir,ith) / boltz_sum
         fs2D(ir,ith) = fs2D(ir,ith) / boltz_sum
      end do
   end do

   ! Average the 2D input histograms into 1D
   write(*,*) 'Averaging input histograms from 2D into 1D...'

   allocate( g1D(histDistBins), fr1D(histDistBins), u01D(histDistBins) )
   g1D = 0_dp; fr1D = 0_dp; u01D = 0_dp

   u01D(:) = minval(-g2D(:,:),dim=2)

   do ir = 1, histDistBins
      boltz_sum = 0_dp
      do ith = 1, histCosThBins
         boltz = exp(g2D(ir,ith) + u01D(ir))
         g1D(ir) = g1D(ir) + exp(g2D(ir,ith)) ! g1D is g
         fr1D(ir) = fr1D(ir) + (boltz * fr2D(ir,ith)) ! f.r
         boltz_sum = boltz_sum + boltz ! denominator for averaging over theta
      end do
      g1D(ir) = log(g1D(ir) / real(histCosThBins,dp)) ! finish average over cosTh by dividing and converting to log(g)
      fr1D(ir) = fr1D(ir) / boltz_sum
   end do
   ! After the averaging is done enforce the cutoff
   do ir = 1, histDistBins
      if (g1D(ir).lt.cut) then
         g1D(ir) = cut
         fr1D(ir) = -cut
         ispline1D = ir
      end if
   end do

   ! note: write out the effective input histogram after averaging/alterations.
   write(*,*) 'Writing input histogram after averaging/alterations to "input_hist.out" ...'
   open(91,file='input_hist.out',status='replace')
   write(91,*) '# 1. Distance'
   write(91,*) '# 2. g'
   write(91,*) '# 3. f.r'
   write(91,*) '#'
   write(91,*) '#'
   write(91,*) '#'
   write(91,*) '#'
   do ir = 1, histDistBins
      write(91,*) histDist(ir), g1D(ir), fr1D(ir)
   end do
   close(91)

   ! Spline the log(g) and force arrays using ideal values for the slopes at small r. This populates the second derivative arrays.
   ! This requires ideal values that have been averaged over phi.
   allocate( idealHist2D(5,histDistBins,histCosThBins) )
   call ideal_3D_to_2D(idealHist,histDistBins,histCosThBins,histPhiBins, idealHist2D)

   allocate( idealHist1D(3,histDistBins) )
   call ideal_2D_to_1D(idealHist2D,histDistBins,histCosThBins, idealHist1D)

   allocate( g1D2(histDistBins), fr1D2(histDistBins) )
   g1D2 = 0_dp; fr1D2 = 0_dp

   call spline(histDist,g1D,ispline1D,histDistBins,idealHist1D(2,ispline1D),real(0,dp), g1D2)
   call spline(histDist,fr1D,ispline1D,histDistBins,idealHist1D(3,ispline1D),real(0,dp), fr1D2)

   !debug   difference between ideal and measured 
   rTmp = int(real(4,dp)/histDistStepSize)
   print*, rTmp, histDistStepSize
   do ith=1,histCosThBins
      write(55,*) histCosTh(ith), g2D(rTmp,ith), idealhist2D(1,rTmp,ith)
   end do

   !debug
   do ir=1,histdistbins
      write(45,*) histDist(ir), g1D(ir), g1D2(ir), idealHist1D(1,ir), fr1D(ir), idealHist1D(2,ir)
   end do
   !do ir=1,100*histdistbins
      !xx = ir*(histDistStepSize/100_dp)
      !call splint(histDist,g1D,g1D2,histDistBins,xx, yy)
      !write(55,*) xx, yy
   !end do
   !write(65,*) histDist(ispline1D), g1D(ispline1D)

end subroutine spline_hist_array


! read LJ--LJ displacements from file
subroutine R_list
   use cfgData
   use ctrlData
   implicit none
   integer            :: ios, i
   character(len=16)   :: junk
   character(len=64)   :: line

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
   use cfgData
   use angleData
   use ctrlData
   implicit none
   integer         :: i
   real(kind=dp)   :: psiLF

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
   allocate( R_axis(cfgRBins), fAvg(cfgRBins), x_axis(xBins), z_axis(zBins) )
   R_axis = 0_dp; fAvg = 0_dp; x_axis = 0_dp; z_axis = 0_dp

   ! allocate arrays for control arrays
   allocate( frcSPA(1, xBins, zBins), grSPA(xBins, zBins) )

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
   use cfgData
   use histData
   use angleData
   use ctrlData
   use constants
   use functions
   implicit none
   integer         :: r, i, j, ip, tid, omp_get_thread_num, omp_get_num_threads
   real(kind=dp)   :: gx, gx2, fx

   write(*,*) "Computing average force..."
   flush(6)

   ! Calculate the average force integral for top half of bisecting plane of cylinder
   do r = 1, cfgRBins ! loop lj--lj distances
      frcSPA = 0_dp; grSPA = 0_dp
      !$omp PARALLEL DEFAULT( none ) &
      !$omp PRIVATE( ip, i, j, gx, gx2, fx ) &
      !$omp SHARED( r, xBins, zBins, cut, R_axis, x_axis, z_axis, histDist, histDistBins, g1D, g1D2, fr1D, fr1D2, frcSPA, grSPA )
      !!$omp   NUM_THREADS( 1 )
      if ((omp_get_thread_num().eq.0).and.(r.eq.1)) then
         write(*,*) 'Parallel CPUs:         ', omp_get_num_threads()
         flush(6)
      end if
      !$omp DO SCHEDULE( guided )
      do ip = 1, (xBins*zBins)
         ! Convert single index 'ip' to the x and z indicies 'i' and 'j' respectively.
         i = int((ip-1)/zBins)+1      ! x integer
         j = mod(ip-1,zBins)+1      ! z integer
          
         rSolv1(1) = -R_axis(r)/2_dp - x_axis(i)
         rSolv1(2) = 0_dp
         rSolv1(3) = -z_axis(j)
         rSolvn(1) = norm2(rSolv1)
         rSolv2(1) = R_axis(r)/2_dp - x_axis(i)
         rSolv2(2) = 0_dp
         rSolv2(3) = -z_axis(j)
         rSolvn(2) = norm2(rSolv2)


         call splint(histDist,g1D,g1D2,histDistBins,rSolvn(1), gx)   ! solute 1
         call splint(histDist,g1D,g1D2,histDistBins,rSolvn(2), gx2)   ! solute 2
         if (rSolvn(1).gt.histDist(histDistBins)) then
            gx = 0_dp
         end if
         if (rSolvn(2).gt.histDist(histDistBins)) then
            gx2 = 0_dp
         end if
         gx = exp(gx+gx2)
          
         if (gx .gt. 1d-6) then ! if gx == 0 then don't waste time with the rest of the calculation
            if (rSolvn(1).gt.histDist(histDistBins)) then
               fx = 0_dp
            else
               call splint(histDist,fr1D,fr1D2,histDistBins,rSolvn(1), fx)
            end if
             
            frcSPA(1,i,j) = frcSPA(1,i,j) + (gx * fx * (-rSolv1(1)/rSolvn(1)))   ! (f.r)*g.R^{hat}
            grSPA(i,j) = grSPA(i,j) + gx
         end if
      end do !ip
      !$omp END DO
      !$omp END PARALLEL
       
      ! Add each cell forces to average and normalize
      do i = 1, xBins
         do j = 1, zBins
            fAvg(r) = fAvg(r) + (frcSPA(1,i,j) * z_axis(j))
            frcSPA(1,i,j) = frcSPA(1,i,j)/grSPA(i,j)
         end do !z again
      end do !x again
      call write_test_out(r) ! write grSPA and frcSPA arrays

      ! NOTE : After the fact multiply all elements by 2*pi*density/8/pi/pi (2*2pi*pi/3 (4pi**2)/3 steradians from orientations)
      !       Number density of chloroform per Angstrom**3 == 0.00750924
      fAvg(r) = fAvg(r)*2*pi*density*xzStepSize*xzStepSize
   end do !r

end subroutine compute_avg_force


! integrate the average force from 'compute_avg_force' to get the PMF.
subroutine integrate_force
   use cfgData
   use ctrlData
   implicit none
   integer       :: d

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
subroutine write_test_out(r)
   use cfgData
   use ctrlData
   implicit none
   integer            :: r, i_f, i, j
   character(len=32)   :: temp, filename
   character(len=8)   :: frmt

   i_f = (r-1) * int(RStepSize*10)
   frmt = '(I3.3)' ! an integer of width 3 with zeroes on the left
   write(temp,frmt) i_f ! converting integer to string using 'internal file'
   filename='hist1D_output.'//trim(temp)//'.dat'


   open(35,file=filename,status='replace')
   write(6,*) "Writing test file:   ", filename
   write(35,*) "# 1.   X Distance"
   write(35,*) "# 2.   Z Distance"
   write(35,*) "# 3.   g(r)"
   write(35,*) "# 4.   Force.r"
   write(35,*) "# "
   write(35,*) "# "
   do j = 1, zBins
      do i = 1, xBins
         write(35,898) x_axis(i), z_axis(j), grSPA(i,j), frcSPA(1,i,j)
      end do
   end do
   close(35)

   flush(6)
   
898      format (4(1x,es14.7))

end subroutine write_test_out


! write output file
subroutine write_output(outFile)
   use cfgData
   implicit none
   character(len=64)    :: outFile
   integer          :: r

   open(35,file=outFile,status='replace')
   write(6,*) "Writing output file:   ", outFile
   write(35,*) "# 1.   R Distance"
   write(35,*) "# 2.   Avg Force"
   write(35,*) "# 3.   PMF"
   do r = 1, cfgRBins
      write(35,899) R_axis(r), fAvg(r), u_dir(r)
   end do
   close(35)

   flush(6)

899      format (3(1x,es14.7)) ! scientific format

end subroutine write_output
