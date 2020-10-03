program Fingerprint
  ! Main program to test silicon potential subroutines
  implicit real*8 (a-h,o-z)
  parameter(nat=8,natx_sphere=100,ns=1,np=1,nconf=99)
  dimension rxyz(3,nat),rcov(nat),alat(3,3)
  dimension fpall((ns+3*np)*natx_sphere,nat,nconf)
  character(len=2), dimension(nat,nconf) :: symb
  character(len=100) :: filetype
  integer, dimension(nat) :: symb_table
  integer :: progress
  symb_table = (/0, 0, 0, 1, 1, 1, 1, 1 /)
  progress = nconf/10

  !filetype='ascii'
  filetype='ascii'
  if (filetype.eq.'ascii') then
     open(unit=10,file='POSALL.ascii')
  else if (filetype.eq.'xyz') then
     open(unit=10,file='POSALL.xyz')
  else 
      stop 'unrecognized file format'
  endif

do iconf=1,nconf
  if (modulo(iconf, progress).eq.0) then
    write(*,'(F10.1,A)') (float(iconf)/float(nconf))*100,"%"
  endif
	call readfiles(filetype,nat,natp,rxyz,alat,rcov,symb(1,iconf))
	if (natp.ne.nat) stop 'MAIN: natx'
	if (filetype.eq.'ascii') then
		! periodic structures
		call back2cell(nat,rxyz,alat)
	else
		! non-periodic structures
		do j=1,3
			do i=1,3
				alat(i,j)=0.d0
			enddo
		enddo
	endif


	call fingerprint_eval(nat, natx_sphere, ns,np, alat, rxyz, rcov, fpall(1,1,iconf), symb_table)
enddo
close(10)


open(unit=20, file='FP.dat')

do iconf=1, nconf
	do iat=1,nat
		open(unit=10,file='POSALL.ascii')
		
		write(20,*) fpall(:,iat,iconf)
			
		
		!write(*,*) "end of vector"
	enddo
enddo
	
close(10)



end program




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   S U B R O U T I N E S   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fingerprint_eval(nat, natx_sphere, ns,np, alat, rxyz, rcov, fp, symb_table)
      implicit real*8 (a-h,o-z)
      parameter(nwork=100)
      dimension workalat(nwork)
      dimension rxyz_sphere(3,natx_sphere),rcov_sphere(natx_sphere),indat(natx_sphere)
      dimension fp((ns+3*np)*natx_sphere,nat),amplitude(natx_sphere),deramplitude(natx_sphere)
      dimension rxyz(3,nat),rcov(nat)!,eval((ns+3*np)*natx_sphere)
      dimension alat(3, 3),alatalat(3,3),eigalat(3)
      dimension alpha(natx_sphere), cs(10),cp(10)
      integer, dimension(nat) :: symb_table
      real*8,allocatable ::  ovrlp(:,:),ovrla(:,:),eval(:)


! parameters for cutoff function: width_cutoff is the width of the Gauusian
! approximated by a polynomial with exponent nex_cutoff
     

    if (alat(1,1)*alat(2,2)*alat(3,3).eq.0.d0 ) then 
        ixyzmax=0
    else

      do i=1,3
      do j=1,3
          alatalat(i,j)=alat(1,i)*alat(1,j)+alat(2,i)*alat(2,j)+alat(3,i)*alat(3,j)
      enddo
      enddo
      call dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
      ixyzmax= int(sqrt(1.d0/eigalat(1))*radius_cutoff) + 1
    endif


       do i=1,10
         cs(i)=sqrt(2.d0)**(i-1)
         cp(i)=sqrt(2.d0)**(i-1)
       enddo

! loop over all center atoms
      natsmax=0
      natsmin=1000000
    do lat=1,nat
!---------------------------------------------------------------------
      if (symb_table(lat).eq.0) then
        width_cutoff = 4.5d0
      else
      width_cutoff = 3.0d0
      endif
!---------------------------------------------------------------------
      nex_cutoff=2
      radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff  ! radius where the polynomial is zero
      

       call atoms_sphere(width_cutoff,nex_cutoff,lat,llat,ixyzmax,nat,natx_sphere,nat_sphere,alat,rxyz,rxyz_sphere, & 
                        rcov,rcov_sphere,indat,amplitude,deramplitude)
       natsmax=max(natsmax,nat_sphere)
       natsmin=min(natsmin,nat_sphere)
!       call xyz2eval(nat_sphere,rxyz_sphere,rcov_sphere,ns,np,amplitude,fp(1,lat))

       norb= nat_sphere*(ns+np*3)
       allocate(ovrlp(norb,norb),ovrla(norb,norb),eval(norb))

        do iat=1,nat_sphere
           alpha(iat)=.5d0/(1.0d0*rcov_sphere(iat))**2
        enddo

        call crtovrlp(nat_sphere,rxyz_sphere,alpha,cs,cp,ns,np,ovrlp)
        call multamp(nat_sphere,ovrlp,amplitude,norb,ns,np,ovrla)
         trace=0.d0
         do i=1,norb
         trace=trace+ovrla(i,i)
         enddo

        call diagonalizematrix(norb, ovrla, eval)


! eigenvalues in decreasing order
         evals=0.d0
         do i=1,norb
         evals=evals+eval(norb-i+1)
         fp(i,lat)=eval(norb-i+1)
         enddo
         norbx=natx_sphere*(ns+np*3)
         do i=norb+1,norbx
         fp(i,lat)=0.d0
         enddo

        deallocate(ovrlp,ovrla,eval)


    enddo

  end subroutine 


subroutine diagonalizeMatrix(n, aa, eval)
  implicit real*8 (a-h,o-z)
  
  ! Calling arguments
  real*8,dimension(n,n) :: aa
  real*8,dimension(n) :: eval
  
  ! Local variables
  real*8 ,dimension(:),allocatable:: work
  
  lwork=100*n
  allocate(work(lwork))
  call dsyev('v','l', n, aa, n, eval, work, lwork, info)
  if(info/=0) stop ' ERROR in dsyev'
  deallocate(work)

end subroutine diagonalizeMatrix


