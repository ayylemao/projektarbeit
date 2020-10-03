program Fingerprint
  ! Main program to test silicon potential subroutines
  implicit real*8 (a-h,o-z)
  parameter(nat=8,natx_sphere=100,ns=1,np=1,nconf=99)
  dimension rxyz(3,nat),rcov(nat),alat(3,3), rcov_fuck(nat)
  dimension fpall((ns+3*np)*natx_sphere,nat,nconf)
  dimension ref_vector((ns+3*np)*natx_sphere)
  real*8 N_rcov
  real*8 C_rcov
  real*8, dimension(2) :: radii
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

open(unit=22, file='fpvector.dat')
do l=1,(ns+3*np)*natx_sphere
  read(22,*) ref_vector(l)
enddo
close(22)

call readfiles(filetype,nat,natp,rxyz,alat,rcov,symb(1,1))
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

C_rcov = rcov(1)
N_rcov = rcov(5)
radii(1) = C_rcov
radii(2) = N_rcov

iconf = 1
do ia = 1, 2
  rcov(1) = radii(ia)
  do ib = 1,2
    rcov(2) = radii(ib)
    do ic = 1,2
      rcov(3) = radii(ic)
      do id = 1,2
        rcov(4) = radii(id)
        do ie = 1,2
          rcov(5) = radii(ie)
          do ig = 1,2
            rcov(6) = radii(ig)
            do ih = 1,2
              rcov(7) = radii(ih)
              do ij = 1,2
                rcov(8) = radii(ij)
                !write(*,*) rcov
                fpd = 0.0d0
                call fingerprint_eval(nat, natx_sphere, ns,np, alat, rxyz, rcov, fpall(1,1,iconf))
                do l=1,(ns+3*np)*natx_sphere
                  fpd=fpd+(fpall(l,1,iconf)-ref_vector(l))**2
                enddo
                if (fpd == 0.0) then
                  write(*,*) "FP distance", fpd
                  write(*,*) ia, ib, ic, id, ie, ig, ih, ij
                  call EXIT(0)
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
enddo






! do iconf=1,1
! 	do i = 1, 5
!     rcov(i) = 5.0
!     call fingerprint_eval(nat, natx_sphere, ns,np, alat, rxyz, rcov, fpall(1,1,iconf))
!     do iat=1,1
!       !open(unit=10,file='POSALL.ascii')
      
!     enddo
!   enddo
! enddo
! close(10)





	

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
    !write(*,*) 'ixyzmax ',ixyzmax


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
    !write(*,*) 'min, max number of atoms in sphere ',natsmin,natsmax

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


