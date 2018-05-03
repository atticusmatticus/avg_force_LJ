! USAGE: ./this_file.x -first [first file number] -last [last file number] -out [out file]
!
! create list of LJ--LJ displacements from explicit .crd files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!   Modules   !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! precision module
module prec
	integer, parameter	:: dp = kind(1.0d0)
endmodule prec

! CRD file name prefix
module filePrefix
	character(len=32)	:: prefix = 'crd_files/LJ2_d0.00_us.'
endmodule filePrefix

! command line data
module cmdData
	integer	:: ifirst, ilast, istep
endmodule cmdData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!  Main Program  !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program compute_avgForce
	use prec
	implicit none
	character(len=64) 	:: cfirst, clast, cstep, outFile
	real(kind=dp) 		:: omp_get_wtime, ti, tf

	ti = omp_get_wtime()

	! parse command line for output file and file range.
	call parse_command_line(cfirst, clast, cstep, outFile)

	! 
	call iterate_files(outFile)
	
	! Print time taken to finish calculation.
	tf = omp_get_wtime()
	write(*,*) "Total time elapsed: ", tf-ti, "seconds"

endprogram compute_avgForce


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Subroutines !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! parse commandline for relevant files and information.
subroutine parse_command_line(cfirst, clast, cstep, outFile)
	use cmdData
	implicit none
	character(len=64) 	:: cfirst, clast, cstep, outFile
	character(len=16) 	:: arg
	integer 			:: i
	logical 			:: firstFlag, lastFlag, stepFlag, outFileFlag

	firstFlag = .false.
	lastFlag = .false.
	stepFlag = .false.
	outFileFlag = .false.
	i=1
	do
		call get_command_argument(i,arg)
		select case (arg)

		case ('-first')
			i = i+1
			call get_command_argument(i,cfirst)
			firstFlag=.true.
			read(cfirst,*) ifirst
			print*, "first: ", ifirst
		case ('-last')
			i = i+1
			call get_command_argument(i,clast)
			lastFlag=.true.
			read(clast,*) ilast
			print*, "last: ", ilast
		case ('-step')
			i = i+1
			call get_command_argument(i,cstep)
			stepFlag=.true.
			read(cstep,*) istep
			print*, "step size: ", istep
		case ('-out')
			i = i+1
			call get_command_argument(i,outFile)
			outFileFlag=.true.
			print*, "out File: ", outFile
		case default
			print*, 'Unrecognized command-line option: ', arg
			print*, 'Usage: crd_list.x -first [first number]'
			print*, 'Usage: crd_list.x -last [last number]'
			print*, 'Usage: crd_list.x -step [step size]'
			print*, 'Usage: crd_list.x -out [out file]'
			stop

		end select
		i = i+1
		if (i.ge.command_argument_count()) exit
	enddo

	if (firstFlag.eqv..false.) then
		write(*,*) "Must provide a first number using command line argument -first [first number]"
		stop
	endif

	if (lastFlag.eqv..false.) then
		write(*,*) "Must provide a last number using command line argument -last [last number]"
		stop
	endif

	if (stepFlag.eqv..false.) then
		write(*,*) "Must provide a step size using command line argument -step [step size]"
		stop
	endif

	if (outFileFlag.eqv..false.) then
		write(*,*) "Must provide a out file using command line argument -out [out file name]"
		stop
	endif

endsubroutine parse_command_line

! iterate through .crd files
subroutine iterate_files(outFile)
	use prec
	use filePrefix
	use cmdData
	implicit none
	integer				:: i
	character(len=64)	:: temp, filename, outFile
	character(len=6)	:: frmt
	real(kind=dp)		:: dist

	open(25,file=outFile,status='replace')
	do i = ifirst, ilast, istep
		if ((i .ge. 0) .and. (i .le. 9)) then
			frmt = '(I1)' ! an integer with width 1, no leading zeroes
		else if ((i .ge. 10) .and. (i .le. 99)) then
			frmt = '(I2)' ! an integer with width 2, no leading zeroes
		else if ((i .ge. 100) .and. (i .le. 999)) then
			frmt = '(I3)' ! an integer with width 3, no leading zeroes
		endif
		write(temp,frmt) i ! converting integer to string using 'internal file'
		filename = trim(prefix)//trim(temp)//'.crd'
		print*, 'Reading file: ', filename

		call read_crd(filename, dist)

		write(25,999) i, dist
	enddo
	close(25)

999 format (1x,i3,1x,f16.12)

endsubroutine iterate_files

! read the .crd files
subroutine read_crd(filename, dist)
	use prec
	implicit none
	integer				:: i
	character(len=64)	:: filename
	real(kind=dp)		:: x, y, z, xT, yT, zT, dist
	character(len=16)	:: junk

	xT = 0_dp; yT = 0_dp; zT = 0_dp
	
	open(20,file=filename,status='old')
	do i = 1, 3
		if (i .eq. 1) then
			read(20,*) junk, junk
		else
			read(20,*) junk, x, y, z
			xT = xT + abs(x)
			yT = yT + abs(y)
			zT = zT + abs(z)
		endif
	enddo
	close(20)
	
	dist = dsqrt(xT*xT + yT*yT + zT*zT)

endsubroutine read_crd
