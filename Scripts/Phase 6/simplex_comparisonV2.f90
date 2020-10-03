program simplex_comparison
	implicit real*8 (a-h, o-z)
	integer iconf, iat, ientry, icorrect, nconfig
	dimension ref_vec(400), fp_vec(1600, 400), results(400), ref_matrix(101, 400)
	integer, dimension(101, 3) :: iconfigs
	character(1) symb_test, symb_truth
	acorrect = 0.0
	nconfig = 792
	call readVector(fp_vec)
	
	do i = 1, 101
		call readConfigurations(iconfigs)
		iconf = iconfigs(i, 3)
		iat = iconfigs(i, 2)
		call getRefVector(iconf, iat, ref_vec)
		ref_matrix(i, :) = ref_vec
	enddo

	do k = 1, nconfig		
		bmin = 100.0
		do i = 1, 101
			ref_vec = ref_matrix(i, :)
			iat = iconfigs(i, 2)
			iconf = iconfigs(i, 3)
			fpd = 0.0
			do j = 1, 100
				fpd = fpd + (ref_vec(j) - fp_vec(k ,j))**2
			enddo
			fpd = sqrt(fpd)
			if (fpd < bmin) then
				bmin = fpd
				call getAtsymb_test(iconf, iat, symb_test)
				!write(*,*) symb_test, iconfigs(i, 1), iat, iconf, bmin
			end if
		enddo
		!write(*,*) "------------------------"
		
		call getAtsymb_truth(k, symb_truth)
		if (symb_truth == symb_test) then
			acorrect = acorrect + 1.0
			write(*,*) "Correct Guess: ", symb_test
		end if
	enddo
	aRatio = acorrect/nconfig
	write(*,*) "Ratio of correct guesses:", aRatio
	

end program




subroutine getRefVector(iconf, iat, fp_vec)
  implicit real*8 (a-h, o-z)
  integer iconf, iat, ientry
  dimension fp_vec(100)

 
  ientry = (iconf-1)*8+iat

  open(unit=20, file="FP.dat")
  do i= 1, ientry-1
  	read(20, *)
  enddo
  read(20, *) fp_vec
  close(20)
end subroutine getRefVector


subroutine readConfigurations(iconfigs)
	implicit real*8 (a-h, o-z)
	integer, dimension(101, 3) :: iconfigs

	open(21, file="configs")
	do i = 1, 101
		read(21, *) iconfigs(i, :)
	enddo
	close(21)
end subroutine readConfigurations


subroutine getAtsymb_test(iconf, iat, symb_test)
	implicit real*8 (a-h, o-z)
	integer iconf, iat, natom
	character(1) symb_test
	natom = 8*(iconf-1)+iat
	if (modulo(natom, 8) < 4 .and. modulo(natom, 8) /= 0) then
		symb_test = "C"
	else
		symb_test = "N"
	end if
end subroutine getAtsymb_test


subroutine readVector(fp_vec)
	implicit real*8 (a-h, o-z)
	integer ivec
	dimension fp_vec(1600, 400)
	open(22, file="FP.dat")
	do i = 1, 1600
		read(22,*) fp_vec(i, :)
	enddo
end subroutine readVector


subroutine getAtsymb_truth(iatom, symb_truth)
	implicit real*8 (a-h, o-z)
	character(1) symb_truth
	if (modulo(iatom, 16) < 7 .and. modulo(iatom, 16) /= 0) then
		symb_truth = "C"
	else
		symb_truth = "N"
	end if
end subroutine getAtsymb_truth


