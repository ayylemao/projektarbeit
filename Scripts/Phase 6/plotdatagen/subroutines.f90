subroutine atoms_sphere(width_cutoff,nex_cutoff,lat,llat,ixyzmax,nat,natx_sphere,nat_sphere,alat,rxyz,rxyz_sphere, & 
                        rcov,rcov_sphere,indat,amplitude,deramplitude)
  implicit real*8 (a-h,o-z)
  dimension rxyz_sphere(3, natx_sphere),rcov_sphere(natx_sphere)
  dimension amplitude(natx_sphere),deramplitude(natx_sphere)
  dimension rxyz(3,nat),rcov(nat)
  dimension alat(3, 3),indat(natx_sphere)
  

  radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff
  radius_cutoff2=radius_cutoff**2
  factor_cutoff=1.d0/(2.d0*nex_cutoff*width_cutoff**2)
!!  write(*,*) 'width_cutoff,radius_cutoff',width_cutoff,radius_cutoff

        
     nat_sphere=0
     do jat = 1, nat
         do ix = -ixyzmax,ixyzmax
           do iy = -ixyzmax,ixyzmax
             do iz = -ixyzmax,ixyzmax
                xj = rxyz(1, jat) + ix*alat(1,1)+iy*alat(1,2)+iz*alat(1,3)
                yj = rxyz(2, jat) + ix*alat(2,1)+iy*alat(2,2)+iz*alat(2,3)
                zj = rxyz(3, jat) + ix*alat(3,1)+iy*alat(3,2)+iz*alat(3,3)
                dist2 = (xj-rxyz(1, lat))**2+(yj-rxyz(2, lat))**2+(zj-rxyz(3, lat))**2
                !write(*,*) xj,rxyz(1, lat),yj,rxyz(2, lat),zj,rxyz(3,lat)


                if (dist2.le.radius_cutoff2) then
                    nat_sphere=nat_sphere+1
                    if (nat_sphere.gt.natx_sphere) stop 'enlarge natx_sphere'
                    !amplitude(nat_sphere)=(1.d0 - dist2*factor_cutoff)**nex_cutoff
                    temp = (1.d0 - dist2*factor_cutoff)**(nex_cutoff-1)
                    amplitude(nat_sphere) = temp*(1.d0 - dist2*factor_cutoff)
                    deramplitude(nat_sphere) = -2.d0*factor_cutoff*nex_cutoff*temp

                    rxyz_sphere(1,nat_sphere)=xj
                    rxyz_sphere(2,nat_sphere)=yj
                    rxyz_sphere(3,nat_sphere)=zj
                    rcov_sphere(nat_sphere)=rcov(jat)
                    indat(nat_sphere)=jat
                    if (dist2.eq.0.d0) llat=nat_sphere
                endif
             enddo
           enddo
         enddo
     enddo
end subroutine atoms_sphere


subroutine crtovrlp(nat,rxyz,alpha,cs,cp,ns,np,ovrlp, CN_symbs, const)
  implicit real*8 (a-h,o-z)
  real*8  rxyz(3,nat)
  real*8 ovrlp(nat*(ns+3*np),nat*(ns+3*np))
  real*8  alpha(nat), cs(10),cp(10)
  real*8 C_s, C_p, N_s, N_p, const
  integer, dimension(natx_sphere) :: CN_symbs
  integer :: natx_sphere
  natx_sphere = 100
  C_s = 0.500866
  C_p = 0.199186
  N_s = 0.676151
  N_p = 0.266297

  if(ns>10 .or. np > 10) stop 'ns > 10   .or.  np > 10  !'


 ! 1- setup the overlap matrix 

  !  <s|s>
  do jat=1,nat
    do js=1,ns
      jorb=(jat-1)*(ns+3*np)+js
      if (CN_symbs(jat).eq.0) then
            eps = C_s
          else
            eps = N_s
      endif
      aj=const*alpha(jat)/(cs(js)*eps)
      xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

      do iat=1,nat
        do is=1,ns
          !!iorb=iat+(is-1)*nat
          iorb=(iat-1)*(ns+3*np)+is
          if (CN_symbs(iat).eq.0) then
            eps = C_s
          else
            eps = N_s
          endif
            ai= const*alpha(iat)/(cs(is)*eps)

          xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

          xij=xi-xj; yij=yi-yj; zij=zi-zj
          r2=xij**2 + yij**2 + zij**2 
          t1=ai*aj
          t2=ai+aj

          ! normalized GTOs:
          sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2) 
          ovrlp(iorb,jorb)=sij

        enddo
      enddo
    enddo  
  enddo  


  !  <pi|sj>
  do jat=1,nat
    do js=1,ns
      if (CN_symbs(jat).eq.0) then
            eps = C_s
          else
            eps = N_s
      endif
      jorb=(jat-1)*(ns+3*np)+js
      aj=const*alpha(jat)/(cs(js)*eps)
      xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

      do iat=1,nat
        do ip=1,np
          if (CN_symbs(iat).eq.0) then
            eps = C_p
          else
            eps = N_p
          endif
          !!iorb=1+(iat-1)*3+ns*nat + (ip-1)*3*nat
          iorb=(iat-1)*(ns+3*np)+ns+ip
          ai= const*alpha(iat)/(cp(ip)*eps)
          xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

          xij=xi-xj; yij=yi-yj; zij=zi-zj
          r2=xij**2 + yij**2 + zij**2 

          t1=ai*aj
          t2=ai+aj

          ! normalized GTOs:
          sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2) 
          t3=-2.d0*sqrt(ai)*aj/t2
          ovrlp(iorb  ,jorb  )= t3 * xij *sij
          ovrlp(iorb+1,jorb  )= t3 * yij *sij
          ovrlp(iorb+2,jorb  )= t3 * zij *sij

        enddo
      enddo
    enddo  
  enddo  


  !  <si|pj> 
  do jat=1,nat
    do jp=1,np
      if (CN_symbs(jat).eq.0) then
          eps = C_p
        else
          eps = N_p
      endif
      jorb=(jat-1)*(ns+3*np)+ns+jp
      aj=const*alpha(jat)/(cp(jp)*eps)
      xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

      do iat=1,nat
        do is=1,ns
          if (CN_symbs(iat).eq.0) then
              eps = C_s
            else
              eps = N_s
          endif
          !!iorb=iat+(is-1)*nat
          iorb=(iat-1)*(ns+3*np)+is
          ai= const*alpha(iat)/(cs(is)*eps)
          xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

          xij=xi-xj; yij=yi-yj; zij=zi-zj
          r2=xij**2 + yij**2 + zij**2 

          t1=ai*aj
          t2=ai+aj

          ! normalized GTOs:
          sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2) 
          t3=+2.d0*sqrt(aj)*ai/t2
          ovrlp(iorb,jorb  )= t3 * xij *sij
          ovrlp(iorb,jorb+1)= t3 * yij *sij
          ovrlp(iorb,jorb+2)= t3 * zij *sij

        enddo
      enddo
    enddo  
  enddo  


  !  <p|p>
  do jat=1,nat
    do jp=1,np
      if (CN_symbs(jat).eq.0) then
          eps = C_p
        else
          eps = N_p
      endif
      jorb=(jat-1)*(ns+3*np)+ns+jp
      aj=const*alpha(jat)/(cp(jp)*eps)
      xj=rxyz(1,jat) ; yj=rxyz(2,jat); zj=rxyz(3,jat)

      do iat=1,nat
        do ip=1,np
          if (CN_symbs(iat).eq.0) then
            eps = C_p
          else
            eps = N_p
          endif
          iorb=(iat-1)*(ns+3*np)+ns+ip
          ai= const*alpha(iat)/(cp(ip)*eps)
          xi=rxyz(1,iat) ; yi=rxyz(2,iat); zi=rxyz(3,iat)

          xij=xi-xj; yij=yi-yj; zij=zi-zj
          r2=xij**2 + yij**2 + zij**2 
          t1=ai*aj
          t2=ai+aj

          ! normalized GTOs:
          sij=sqrt(2.d0*sqrt(t1)/t2)**3 * exp (-t1/t2*r2) 
          t4= 2.d0*sqrt(t1)/t2 
          t5=-2.d0*t1/t2 

          ovrlp(iorb  ,jorb  )= t4 *(1.d0 + t5* xij* xij)  * sij
          ovrlp(iorb+1,jorb  )= t4 *(       t5* yij* xij)  * sij
          ovrlp(iorb+2,jorb  )= t4 *(       t5* zij* xij)  * sij
          ovrlp(iorb  ,jorb+1)= t4 *(       t5* xij* yij)  * sij
          ovrlp(iorb+1,jorb+1)= t4 *(1.d0+  t5* yij* yij)  * sij
          ovrlp(iorb+2,jorb+1)= t4 *(       t5* zij* yij)  * sij
          ovrlp(iorb  ,jorb+2)= t4 *(       t5* xij* zij)  * sij
          ovrlp(iorb+1,jorb+2)= t4 *(       t5* yij* zij)  * sij
          ovrlp(iorb+2,jorb+2)= t4 *(1.d0+  t5* zij* zij)  * sij

        enddo
      enddo
    enddo  
  enddo  

end subroutine crtovrlp

subroutine multampspd(nat,ovrlp,amplitude,norb,ns,np,nd,ovrla)
  implicit real*8 (a-h,o-z)
  real*8 amplitude(nat), ovrlp(norb,norb),ovrla(norb,norb)

  do jat = 1,nat
    do j = 1,(ns+3*np+5*nd)
      jorb = (jat -1)*(ns+3*np+5*nd) + j
      do iat = 1,nat
        do i = 1,(ns+3*np+5*nd)
          iorb = (iat-1)*(ns+3*np+5*nd) +i
          ovrla(iorb,jorb) = ovrlp(iorb,jorb)*amplitude(iat)*amplitude(jat)
        end do    
      end do  
    end do
  end do 
  
end subroutine multampspd

subroutine multampoff(nat,ovrlp,amplitude,norb,ns,np,ovrla)
  implicit real*8 (a-h,o-z)
  real*8 amplitude(nat), ovrlp(norb,norb),ovrla(norb,norb)

  do jat = 1,nat
    do j = 1,(ns+3*np)
      jorb = (jat -1)*(ns+3*np) + j
      do iat = 1,nat
        do i = 1,(ns+3*np)
          iorb = (iat-1)*(ns+3*np) +i
          if (iat.eq.jat) then
          ovrla(iorb,jorb) = ovrlp(iorb,jorb)
          else
          ovrla(iorb,jorb) = ovrlp(iorb,jorb)*amplitude(iat)*amplitude(jat)
          endif
        end do    
      end do  
    end do
  end do 
  
end subroutine multampoff


subroutine multamp(nat,ovrlp,amplitude,norb,ns,np,ovrla)
  implicit real*8 (a-h,o-z)
  real*8 amplitude(nat), ovrlp(norb,norb),ovrla(norb,norb)

  do jat = 1,nat
    do j = 1,(ns+3*np)
      jorb = (jat -1)*(ns+3*np) + j
      do iat = 1,nat
        do i = 1,(ns+3*np)
          iorb = (iat-1)*(ns+3*np) +i
          ovrla(iorb,jorb) = ovrlp(iorb,jorb)*amplitude(iat)*amplitude(jat)
        end do    
      end do  
    end do
  end do 
  
end subroutine multamp




     subroutine frac2cart(nat, alat, xyzred, rxyz)
     implicit real*8 (a-h,o-z)
     dimension alat(3,3), xyzred(3,nat), rxyz(3,nat)
     
     do iat=1,nat
        do i = 1, 3
           t = 0.d0
           do j = 1, 3
              t = t + xyzred(j,iat) * alat(i, j)
           end do
           rxyz(i,iat) = t
        end do
     enddo

    !  do j=1,3
    !  do i=1,3
    !  alat(i,j)=alat(i,j)
    !  enddo
    !  enddo
     
     end subroutine frac2cart




subroutine cart2frac(nat,alat,rxyz,rxyzred)
  !This subrouine will convert the redernal coordinates into cartesian coordinates
  implicit real*8 (a-h,o-z)
  dimension rxyzred(3,nat),rxyz(3,nat),alat(3,3),alatinv(3,3)

    div=(alat(1,1)*alat(2,2)*alat(3,3)-alat(1,1)*alat(2,3)*alat(3,2)- & 
         alat(1,2)*alat(2,1)*alat(3,3)+alat(1,2)*alat(2,3)*alat(3,1)+ & 
         alat(1,3)*alat(2,1)*alat(3,2)-alat(1,3)*alat(2,2)*alat(3,1))
    div=1.d0/div
      alatinv(1,1) = (alat(2,2)*alat(3,3)-alat(2,3)*alat(3,2))*div
      alatinv(1,2) =-(alat(1,2)*alat(3,3)-alat(1,3)*alat(3,2))*div
      alatinv(1,3) = (alat(1,2)*alat(2,3)-alat(1,3)*alat(2,2))*div
      alatinv(2,1) =-(alat(2,1)*alat(3,3)-alat(2,3)*alat(3,1))*div
      alatinv(2,2) = (alat(1,1)*alat(3,3)-alat(1,3)*alat(3,1))*div
      alatinv(2,3) =-(alat(1,1)*alat(2,3)-alat(1,3)*alat(2,1))*div
      alatinv(3,1) = (alat(2,1)*alat(3,2)-alat(2,2)*alat(3,1))*div
      alatinv(3,2) =-(alat(1,1)*alat(3,2)-alat(1,2)*alat(3,1))*div
      alatinv(3,3) = (alat(1,1)*alat(2,2)-alat(1,2)*alat(2,1))*div

      do iat=1,nat
      rxyzred(1,iat)=alatinv(1,1)*rxyz(1,iat)+alatinv(1,2)*rxyz(2,iat)+alatinv(1,3)*rxyz(3,iat)
      rxyzred(2,iat)=alatinv(2,1)*rxyz(1,iat)+alatinv(2,2)*rxyz(2,iat)+alatinv(2,3)*rxyz(3,iat)
      rxyzred(3,iat)=alatinv(3,1)*rxyz(1,iat)+alatinv(3,2)*rxyz(2,iat)+alatinv(3,3)*rxyz(3,iat)
      enddo

    !  do j=1,3
    !  do i=1,3
    !  alat(i,j)=alat(i,j)
    !  enddo
    !  enddo
 end subroutine






     subroutine back2cell(nat,rxyz,alat)
     implicit real*8 (a-h,o-z)
     dimension rxyz(3,nat),xyzred(3,nat),alat(3,3)

        call cart2frac(nat,alat,rxyz,xyzred)
        do iat=1,nat
        do l=1,3
        xyzred(l,iat)=modulo(xyzred(l,iat),1.d0)
        enddo
        enddo
        call frac2cart(nat,alat,xyzred,rxyz)

     end subroutine   


