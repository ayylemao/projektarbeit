program CheckForces
  ! Main program to test silicon potential subroutines
  implicit real*8 (a-h,o-z)
  parameter(nat=8,natx_sphere=100,ns=1,np=1,nconf=99)
  dimension rxyz(3,nat),rcov(nat),alat(3,3)
  dimension fpall((ns+3*np)*natx_sphere,nat,nconf)
  character(len=2), dimension(nat,nconf) :: symb
  character(len=100) :: filetype


  !filetype='ascii'
  filetype='ascii'
  if (filetype.eq.'ascii') then
     open(unit=10,file='POSALL.ascii')
  else if (filetype.eq.'xyz') then
     open(unit=10,file='POSALL.xyz')
  else 
      stop 'unrecognized file format'
  endif
  ncount = 0
  open(unit=24, file='FP.dat')
  do iconf=1,nconf
  write(*,*) '------------------ ',iconf,' --------------------'
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


     call fingerprint_eval(nat, natx_sphere, ns,np, alat, rxyz, rcov, fpall(1,1,iconf))
    do iat = 1, nat
      write(24,"(1000(1x,e24.17))") (fpall(j,iat,iconf), j=1,(ns+3*np)*natx_sphere)
      ncount = ncount + 1
      write(25,*) ncount,iat,iconf
    enddo
    enddo
  close(10)


     open(unit=20,file='CC.dat')
     open(unit=21,file='CN.dat')
     open(unit=22,file='NN.dat')
  do iconf=1,nconf
  do jconf=iconf+1,nconf

  do iat=1,nat
  do jat=1,nat

     open(unit=10,file='POSALL.ascii')

    fpd=0.d0
    do l=1,(ns+3*np)*natx_sphere
    fpd=fpd+(fpall(l,iat,iconf)-fpall(l,jat,jconf))**2
    enddo
    fpd=sqrt(fpd)

   if (symb(iat,iconf) .eq. 'C' .and. symb(jat,jconf).eq.'C' ) then
       write(20,'(i6,i6,i3,i3,1x,e14.7)') iconf,jconf,iat,jat,fpd
   else if (symb(iat,iconf) .eq. 'N' .and. symb(jat,jconf).eq.'N' ) then
       write(22,'(i6,i6,i3,i3,1x,e14.7)') iconf,jconf,iat,jat,fpd
   else
       write(21,'(i6,i6,i3,i3,1x,e14.7)') iconf,jconf,iat,jat,fpd
   endif


  enddo
  enddo
  enddo
  enddo

!end comment
end program 




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   S U B R O U T I N E S   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fingerprint_eval(nat, natx_sphere, ns,np, alat, rxyz, rcov, fp)
      implicit real*8 (a-h,o-z)
      parameter(nwork=100)
      dimension workalat(nwork)
      dimension rxyz_sphere(3,natx_sphere),rcov_sphere(natx_sphere),indat(natx_sphere)
      dimension fp((ns+3*np)*natx_sphere,nat),amplitude(natx_sphere),deramplitude(natx_sphere)
      dimension rxyz(3,nat),rcov(nat)!,eval((ns+3*np)*natx_sphere)
      dimension alat(3, 3),alatalat(3,3),eigalat(3)
      dimension alpha(natx_sphere), cs(10),cp(10)
      real*8,allocatable ::  ovrlp(:,:),ovrla(:,:),eval(:)


! parameters for cutoff function: width_cutoff is the width of the Gauusian
! approximated by a polynomial with exponent nex_cutoff
      width_cutoff=4.d0
      nex_cutoff=2
! The following line has to be idential to the corresponding line in xyz2devaldr
      radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff  ! radius where the polynomial is zero

    if (alat(1,1)*alat(2,2)*alat(3,3).eq.0.d0 ) then 
        ixyzmax=0
    else

      do i=1,3
      do j=1,3
          alatalat(i,j)=alat(1,i)*alat(1,j)+alat(2,i)*alat(2,j)+alat(3,i)*alat(3,j)
      enddo
      enddo
      call dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
      !  write(*,*) 'alat EVals',eigalat
      !  write(*,*) 'ixyzmax',int(sqrt(1.d0/eigalat(1))*radius_cutoff)
      ixyzmax= int(sqrt(1.d0/eigalat(1))*radius_cutoff) + 1
    endif
    write(*,*) 'ixyzmax ',ixyzmax


       do i=1,10
         cs(i)=sqrt(2.d0)**(i-1)
         cp(i)=sqrt(2.d0)**(i-1)
       enddo

! loop over all center atoms
      natsmax=0
      natsmin=1000000
    do lat=1,nat
       call atoms_sphere(width_cutoff,nex_cutoff,lat,llat,ixyzmax,nat,natx_sphere,nat_sphere,alat,rxyz,rxyz_sphere, & 
                        rcov,rcov_sphere,indat,amplitude,deramplitude)
                   ! write(*,*) lat,rcov(lat),nat_sphere
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
    write(*,*) 'min, max number of atoms in sphere ',natsmin,natsmax

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


