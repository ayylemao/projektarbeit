  subroutine readfiles(filetype,natx,nat,rxyz,alat,rcov,symb)
  implicit real*8 (a-h,o-z)
  dimension rxyz(3, natx),rcov(natx),alat(3, 3)
  character(len=2), dimension(natx) :: symb
  character(len=20) :: ustring
   character(len=100) :: filetype
   logical :: debug

   debug=.false.
   !debug=.false.


  
 if (trim(filetype).eq."ascii") then



        read(10,*) nat!,t1,t2,ustring,energy,ustring,enthalpy,ustring,volume
        if (nat.gt.natx) stop 'natx'
        ustring ='angstroem'
        read(10,*) dxx,dyx,dyy
        read(10,*) dzx,dzy,dzz

        alat(1,1)=dxx
        alat(2,1)=0.d0
        alat(3,1)=0.d0

        alat(1,2)=dyx
        alat(2,2)=dyy 
        alat(3,2)=0.d0

        alat(1,3)=dzx
        alat(2,3)=dzy 
        alat(3,3)=dzz

       if (trim(ustring).eq.'angstroemd0' .or. trim(ustring).eq.'angstroem') then
               convert=1.d0/0.52917720859d0
       else
               convert=1.d0
       endif

        do j=1,3
        do i=1,3
        alat(i,j)=alat(i,j)*convert
        end do
        end do

        do i = 1, nat
           read(10, *) ( rxyz(j, i), j = 1, 3) ,symb(i)
           do l=1,3
           rxyz(l,i)=rxyz(l,i)*convert
           enddo
        end do

!    Assign the covalent radii
        do i = 1, nat
           call sym2rcov(symb(i), rcov(i))
        end do

 end if  !  ASCII
 


 if (trim(filetype).eq."xyz") then

        read(10,*) nat ,ustring, energy
        !read(10,*) nat,energy
        if (nat.gt.natx) stop 'natx'
        read(10,*) 
        !if (debug) write(*,*)  natp,ustring,energy
        !if (natp.ne.nat) stop "Number of atoms not equals the defined parameter!!!!" 
        !read(10,*) !Ignore the "free" line EDIT HERE FOR PERIODIC BOUNDARY

       !ustring= 'atomic' 

       if (trim(ustring).eq.'angstroemd0' .or. trim(ustring).eq.'angstroem') then
               convert=1.d0/0.52917720859d0
       else
               convert=1.d0
               !write(*,*) 'The units are assumed as Bohr!!!'
       endif
               !write(*,*) 'The units are assumed as angstroem'
               !convert=1.d0/0.52917720859d0
        do i = 1, nat
           read(10, *) symb(i),( rxyz(j, i), j = 1, 3) 
           do l=1,3
           rxyz(l,i)=(rxyz(l,i))*convert
           enddo
        end do

        do j=1,3
        do i=1,3
           alat(i,j) = 0.d0
        enddo
        enddo

!    Assign the covalent radii
        do i = 1, nat
           call sym2rcov(symb(i), rcov(i))
        end do
        
 end if  !  XYZ

        end subroutine

