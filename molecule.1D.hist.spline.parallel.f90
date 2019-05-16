!  USAGE: 
!     make molecule1D
!     ./molecule.1D.hist.spline.parallel.x -cfg [cfg_file]
!
!
!  Config File Format:
!
!     gr3File = ../../3_grPDI_3D/d0.00/PDI.d0.00.gr3
!     crdFile = ../../3_grPDI_3D/d0.00/PDI.d0.00.crd
!     gr1CfgFile = ../../9_PDI_hists/d0.00/group.cfg
!     path2gr1Files = ../../9_PDI_hists/d0.00/
!     outFile = molecule.1D.out
!
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!   Modules   !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cfgData
   character(len=64) :: gr3File, crdFile, gr1CfgFile, path2gr1Files, outFile
end module cfgData

module gr3Data
   use prec
   integer :: boxBins, xBins
   real(kind=dp),allocatable :: xyzGr3(:), gGr3(:)
end module gr3Data

module crdData
   use prec
   integer :: nAtoms
   integer,allocatable :: typeCrd(:)
   real(kind=dp),allocatable :: xyzCrd(:,:)
end module crdData

module gr1Data
   use prec
   character(len=64),allocatable :: filePrefix(:)
   character(len=64) :: outFileSuffix
   real(kind=dp),allocatable :: rGr1(:), f(:,:), f2(:,:)
   integer :: gr1Bins
end module gr1Data

module frcData
   use prec
   real(kind=dp),allocatable :: fVec(:,:)
end module frcData


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!  Main Program  !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program molecule_avgForce
   use prec; use cfgData
   implicit none
   character(len=64) :: cfgFile
   real(kind=dp) :: omp_get_wtime
   real(kind=dp) :: ti, tf, seconds
   integer :: hours, minutes

   ti = omp_get_wtime()

   ! read command line
   call parse_command_line(cfgFile)

   ! read config file
   call read_cfg(cfgFile)

   ! import gr3 data
   call make_gr3_table

   ! import solute atom positions
   call make_crd_table

   ! read gr1CfgFile
   call read_gr1_config
   call gr1_forces

   ! calculate force vectors for each solute atom
   call calculate_force_vectors

   ! write output to file
   call write_output

   ! Write time taken to finish calculation.
   tf = omp_get_wtime()

   hours = (tf-ti)/3600
   minutes = mod((tf-ti),3600d0)/60
   seconds = mod(mod((tf-ti),3600d0),60d0)

   write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
   write(*,'(a,i4,a,i2,a,f6.3,a)') "Total time elapsed:   ", hours, "h ", minutes, "m ", seconds, "s"
   flush(6)
end program molecule_avgForce


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse commandline for relevant files.
subroutine parse_command_line(cfgFile)
   implicit none
   character(len=64) :: cfgFile
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
         write(*,*) 'Config File ==> ', cfgFile
         write(*,*) 'Config File Exists ==> ', cfgExist
      case default
         write(*,*) 'Unrecognized command-line option: ', arg
         write(*,*) 'Usage: molecule_avgForce.x -cfg [cfg file]'
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read python cfg file for g(r) parameters
subroutine read_cfg(cfgFile)
   use cfgData
   implicit none
   character(len=64) :: cfgFile
   character(len=128) :: line
   character(len=32) :: firstWord, sep
   integer :: ios
   logical :: gr3FileFlag, crdFileFlag, gr1CfgFileFlag, path2gr1FilesFlag, outFileFlag
   logical :: gr3Exist, crdExist, gr1CfgExist

   gr3FileFlag = .false.; gr3Exist = .false.
   crdFileFlag = .false.; crdExist = .false.
   gr1CfgFileFlag = .false.; gr1CfgExist = .false.
   path2gr1FilesFlag = .false.
   outFileFlag = .false.

   ios = 0

   open(20,file=cfgFile)
   do while(ios>=0)
      read(20,'(a)',IOSTAT=ios) line
      call split(line,'=',firstWord, sep)
      if (line .ne. "") then
         if (firstWord .eq. "gr3File") then
            read(line,'(a)') gr3File
            write(*,*) "Explicit Density File ==> ", trim(gr3File)
            gr3FileFlag = .true.
            INQUIRE(FILE=gr3File, EXIST=gr3Exist) ! check if it exists
         else if (firstWord .eq. "crdFile") then
            read(line,'(a)') crdFile
            write(*,*) "Explicit Coordinate File ==> ", trim(crdFile)
            crdFileFlag = .true.
            INQUIRE(FILE=crdFile, EXIST=crdExist) ! check if it exists
         else if (firstWord .eq. "gr1CfgFile") then
            read(line,'(a)') gr1CfgFile
            write(*,*) "Grouped Atom Forces Config File ==> ", trim(gr1CfgFile)
            gr1CfgFileFlag = .true.
            INQUIRE(FILE=gr1CfgFile, EXIST=gr1CfgExist) ! check if it exists
         else if (firstWord .eq. "path2gr1Files") then
            read(line,'(a)') path2gr1Files
            write(*,*) "Grouped Atom Forces Location ==> ", trim(path2gr1Files)
            path2gr1FilesFlag = .true.
         else if (firstWord .eq. "outFile") then
            read(line,'(a)') outFile
            write(*,*) "Output File ==> ", trim(outFile)
            outFileFlag = .true.
         end if
      end if
   end do
   close(20)

   if (gr3FileFlag.eqv..false.) then
      write(*,*) "Config file must have a 'gr3File' value"
      stop
   else if (gr3Exist.eqv..false.) then
      write(*,*) "Config file must point to a 'gr3File' that exists: ", gr3File, " doesn't exist."
      stop
   end if

   if (crdFileFlag.eqv..false.) then
      write(*,*) "Config file must have a 'crdFile' value"
      stop
   else if (crdExist.eqv..false.) then
      write(*,*) "Config file must point to a 'crdFile' that exists: ", crdFile, " doesn't exist."
      stop
   end if

   if (gr1CfgFileFlag.eqv..false.) then
      write(*,*) "Config file must have a 'gr1CfgFile' value"
      stop
   else if (gr1CfgExist.eqv..false.) then
      write(*,*) "Config file must point to a 'gr1CfgFile' that exists: ", gr1CfgFile, " doesn't exist."
      stop
   end if

   if (outFileFlag.eqv..false.) then
      write(*,*) "Config file must have a 'out_file' value"
      stop
   end if

   flush(6)

end subroutine read_cfg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read gr3 density file
subroutine make_gr3_table
   use gr3Data; use cfgData
   implicit none
   character(len=512) :: line
   character(len=32) :: junk
   integer :: ios, i

   write(*,*) 'Reading Explicit Density File, this may take some time...'; flush(6)
   ! read number of lines in gr3File and allocate that many points in temporary histogram list, histTmp.
   ios = 0; boxBins = -1
   open(20,file=gr3File)
   do while(ios>=0)
      read(20,'(a)',IOSTAT=ios) line
      if (line(1:1) .ne. "#") then
         boxBins = boxBins + 1
      end if
   end do
   close(20)
   xBins = boxBins**(1./3.)
   write(*,*) "Total Bins ==> ", boxBins
   write(*,*) "Cube side length ==> ", xBins, ", (box is assumed to be cubic.)"; flush(6)

   ! populate gr3 arrays
   allocate( xyzGr3(xBins), gGr3(boxBins) )
   ios = 0; i = 1
   open(20,file=gr3File)
   ! read file ignoring comment lines at the beginning
   do while(ios>=0)
      read(20,'(a)',IOSTAT=ios) line
      if ((line(1:1) .ne. "#") .and. (ios .ge. 0)) then
         if (i<=xBins) then
            read(line,*) junk, junk, xyzGr3(i), gGr3(i)
         else
            read(line,*) junk, junk, junk, gGr3(i)
         end if
         i = i + 1
      end if
   end do
   close(20)
end subroutine make_gr3_table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read solute atom coordinates file
subroutine make_crd_table
   use cfgData; use crdData
   implicit none
   character(len=512) :: line
   character(len=64) :: junk
   integer :: a, ios

   ! populate crd array
   ios = 0; a = 1
   open(20,file=crdFile)
   read(20,*) nAtoms, junk ! first line has nAtoms and nTypes. Different format than rest of file
   allocate( xyzCrd(nAtoms,3), typeCrd(nAtoms) )
   do while(ios>=0)
      read(20,'(a)',IOSTAT=ios) line
      if ((line(1:1) .ne. "#") .and. (ios .ge. 0)) then
         !             typeIndex        xCrd         yCrd         zCrd
         !read(line,*) typeCrd(a), xyzCrd(a,1), xyzCrd(a,2), xyzCrd(a,3) ! find the solute atom type from crdFile
         read(line,*) junk, xyzCrd(a,1), xyzCrd(a,2), xyzCrd(a,3) ! find the solute atom type from the gr1CfgFile
         a = a + 1
      end if
   end do
   close(20)
end subroutine make_crd_table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read grouped atom gr1 config file
subroutine read_gr1_config
   use cfgData; use gr1Data; use crdData
   implicit none
   character(len=256) :: line
   character(len=64) :: firstWord, sep
   integer :: ios, nLines, i, atomNumber

   ! get number of lines to allocate arrays
   ios = 0; nLines = -1
   open(20,file=gr1CfgFile)
   do while(ios>=0)
      read(20,'(a)',iostat=ios) line
      if (line(1:1) .ne. '#') then
         if (line .ne. "") then
            nLines=nLines+1
         end if
      end if
   end do
   close(20)
   nLines = nLines - 4  ! four lines for atom prefix & suffix
   allocate( filePrefix(nLines) )

   ! read gr1CfgFile into arrays
   write(*,*) ''; write(*,*) 'Atom Groups:'; write(*,*) '~~~~~~~~~~~~'
   ios = 0; i = 1
   open(20,file=gr1CfgFile)
   do while(ios>=0)
      read(20,'(a)',iostat=ios) line
      call split(line,'=',firstWord, sep)
      !print*, line, ' -- ', firstWord; flush(6)   !debug
      if (line(1:1) .ne. "#") then
         if (line .ne. "") then
            if (firstWord=='file_before_atom_number') then ! not needed
            else if (firstWord=='file_after_atom_number_lj') then ! not needed
            else if (firstWord=='file_after_atom_number_coulomb') then ! not needed
            else if (firstWord=='grouped_file_suffix') then
               read(line,'(a)') outFileSuffix
            else
               read(firstWord,'(a)') filePrefix(i)
               write(*,*) trim(filePrefix(i)), '  ==>  ', trim(line); flush(6)
               do ! loop over all values in list
                  !print*, 'hello ', firstWord, line; flush(6)   !debug
                  call split(line,',',firstWord, sep)
                  !print*, 'hello1 ', firstWord, line; flush(6)   !debug
                  if (line .ne. "") then  ! not the last value
                     read(firstWord,*) atomNumber   ! convert string into integer
                     !print*, 'hey ', atomNumber, size(typeCrd); flush(6)  !debug
                     typeCrd(atomNumber) = i
                  else  ! last value
                     read(firstWord,*) atomNumber   ! convert string into integer
                     typeCrd(atomNumber) = i
                     exit  ! last value so exit the infinite loop
                  end if
               end do
               i = i + 1   ! iterate index into filePrefix
            end if
         end if
      end if
   end do
   close(20)
   !do i = 1, nAtoms  !debug
      !print*, i, typeCrd(i); flush(6)
   !end do
end subroutine read_gr1_config
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make gr1 force arrays
subroutine gr1_forces
   use cfgData; use gr1Data; use crdData; use functions
   implicit none
   character(len=128) :: groupFile, line
   integer :: ios, i, a, imin(size(filePrefix))
   real(kind=dp),allocatable :: fGr1(:,:,:)

   ! get number of lines to allocate arrays
   groupFile = trim(path2gr1Files) // trim(filePrefix(1)) // trim(outFileSuffix)
   gr1Bins = -1; ios = 0
   open(20,file=groupFile)
   do while(ios>=0)
      read(20,'(a)',iostat=ios) line
      if (line(1:1) .ne. "#") then
         if (line .ne. "") then
            gr1Bins=gr1Bins+1
         end if
      end if
   end do
   close(20)
   allocate( rGr1(gr1Bins), fGr1(size(filePrefix),2,gr1Bins), f(size(filePrefix),gr1Bins), f2(size(filePrefix),gr1Bins) )
   rGr1 = 0_dp; fGr1 = 0_dp; f = 0_dp; f2 = 0_dp; imin = 0

   do a = 1, size(filePrefix)
      groupFile = trim(path2gr1Files) // trim(filePrefix(a)) // trim(outFileSuffix)
      i = 1; ios = 0
      open(20,file=groupFile)
      do while((ios>=0).and.(i<=gr1Bins))
         read(20,'(a)',iostat=ios) line
         if (line(1:1) .ne. '#') then
            !             radius           lj      coulomb
            read(line,*) rGr1(i), fGr1(a,1,i), fGr1(a,2,i)
            i = i + 1
         end if
      end do
      close(20)
      ! spline using just the points with data
      ios = 0
      do i = 1, gr1Bins
         f(a,i) = fGr1(a,1,i) + fGr1(a,2,i)  ! total force
         if ((abs(fGr1(a,1,i))>1d-6).and.(ios==0)) then ! grab the first point with non-zero value
            imin(a) = i
            ios = 1
         end if
      end do
      do i = 1, imin(a)-1  !note: linear continuation of force
         f(a,i) = f(a,imin(a)) + (imin(a)-i)*( f(a,imin(a)) - f(a,imin(a)+1) )
      end do
      call spline(rGr1,f(a,:),1,gr1Bins,real(0,dp),real(0,dp), f2(a,:))
      call spline(rGr1,f(a,:),1,gr1Bins,real(0,dp),real(0,dp), f2(a,:))
   end do
end subroutine gr1_forces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! parallel loop through grid points of box and calculate weighted force vector for each solute atom.
subroutine calculate_force_vectors
   use crdData; use gr1Data; use gr3Data; use functions; use frcData
   implicit none
   integer :: a, ip, ixyz(3), j, tid, omp_get_thread_num, omp_get_num_threads, nThreads
   real(kind=dp),allocatable :: fTmpVec(:,:)
   real(kind=dp) :: rVec(3), dist, force

   write(*,*) ''; write(*,*) 'Calculating Forces...'
   !$omp PARALLEL DEFAULT( none ) PRIVATE( tid ) SHARED( nThreads )
   tid = omp_get_thread_num()
   if (tid==0) then
      nThreads = omp_get_num_threads()
   end if
   !$omp END PARALLEL
   allocate( fVec(nAtoms,3), fTmpVec(nThreads,3) )
   fVec = 0_dp

   do a = 1, nAtoms  ! loop over solute atoms
      fTmpVec = 0_dp
      !$omp PARALLEL DEFAULT( none ) &
      !$omp PRIVATE( tid, ip, ixyz, j, rVec, dist, force ) &
      !$omp SHARED( a, boxBins, xBins, xyzGr3, xyzCrd, typeCrd, gGr3, rGr1, f, f2, gr1Bins, fTmpVec ) 
      !!$omp num_threads( 1 )
      tid = omp_get_thread_num()
      if ((tid==0).and.(a==1)) then
         write(*,*) 'Parallel CPUs ==> ', omp_get_num_threads(); flush(6)
      end if
      if (tid==0) then
         write(*,*) 'Solute atoms complete ==> ', a; flush(6)
      end if
      !$omp DO SCHEDULE( guided )
      do ip = 1, boxBins
         ixyz(1) = int((ip-1)/xBins/xBins)+1       ! x index
         ixyz(2) = mod(int((ip-1)/xBins),xBins)+1  ! y index
         ixyz(3) = mod((ip-1),xBins)+1             ! z index
         dist = 0_dp
         do j = 1, 3
            rVec(j) = xyzCrd(a,j) - xyzGr3(ixyz(j))
            dist = dist + rVec(j)**2
         end do
         dist = sqrt(dist)
         if (dist>rGr1(gr1Bins)) then  ! keep splint from extrapolating
            force = 0_dp
         else
            call splint(rGr1,f(typeCrd(a),:),f2(typeCrd(a),:),gr1Bins,dist, force)
         end if
         fTmpVec(tid+1,:) = fTmpVec(tid+1,:) + gGr3(ip)*force*rVec/dist
      end do
      !$omp END DO
      !$omp END PARALLEL
      do j = 1, 3 ! atomic add
         fVec(a,j) = sum(fTmpVec(:,j))
      end do
   end do
end subroutine calculate_force_vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_output
   use cfgData; use crdData; use frcData
   implicit none
   integer :: a

   open(35,file=outFile,status='replace')
   write(6,*) 'Writing output file ==> ', outFile; flush(6)
   write(35,*) '# 1.    Atom Number'
   write(35,*) '# 2.    X Position'
   write(35,*) '# 3.    Y Position'
   write(35,*) '# 4.    Z Position'
   write(35,*) '# 5.    X Force'
   write(35,*) '# 6.    Y Force'
   write(35,*) '# 7.    Z Force'
   do a = 1, nAtoms
      write(35,899) a, xyzCrd(a,1), xyzCrd(a,2), xyzCrd(a,3), fVec(a,1), fVec(a,2), fVec(a,3)
   end do
   close(35)
899   format (i3,6(1x,es14.7))
end subroutine write_output
