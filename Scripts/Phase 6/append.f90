program appenddata
implicit real*8 (a-h, o-z)
parameter(nmax = 2376, nstart = 793)
integer :: conf, iat
open(unit=20, file='fort.25')

conf = 100
iat = 1


do k = nstart, nmax
	 write(20,*) k-792, iat, conf
	 iat = iat + 1
	 if (iat.eq.17) then
	 	iat = 1
	 	conf = conf + 1
	 endif
enddo

	
close(20)

end program appenddata