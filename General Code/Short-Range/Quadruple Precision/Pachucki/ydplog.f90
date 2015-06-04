      subroutine dplog(u,n,gg)
!----------------------------------------------------------------------
!
!                                                     m
!   Purpose           - calculates derivatives (-d/du) [ln(1+u)/u]
!
!----------------------------------------------------------------------
      implicit          real (16) (a-h,o-z)
      parameter         (lenl=200,lenl2=3*lenl)
      dimension         gg(0:lenl2),ff(0:lenl2)
      real (16), parameter :: zero = 0.0q0
      real (16), parameter :: one  = 1.0q0
      real (16), parameter :: threshold  = 0.05q0 

      epsm = epsilon(one)
      au = abs(u)
!
      if (au > threshold) then
         uinv = one/u
         u1inv = one/(one+u)
         ff(0) = log(one+u)
         ff(1) = -u1inv
         gg(0) = ff(0)*uinv
         gg(1) = uinv*(gg(0)+ff(1))
         do i=2,n
            ff(i) = (i-1)*ff(i-1)*u1inv
            gg(i) = uinv*(i*gg(i-1)+ff(i))
         end do
      else if (au == zero) then
         u1inv = one
         ff(0) = zero
         ff(1) = -u1inv
         do i=2,n+1
            ff(i) = (i-1)*ff(i-1)
         end do
!
         do i=n+1,1,-1
            gg(i-1) = -ff(i)/i
         end do
      else
         nadd = int(log(epsm)/log(au))+5
         if (n+nadd > lenl2) then
            write (*,*) 'Error!!! Insufficient array length in DPLOG!'
            stop
         end if
!
         u1inv = one/(one+u)
         ff(0) = log(one+u)
         ff(1) = -u1inv
         do i=2,n+nadd
            ff(i) = (i-1)*ff(i-1)*u1inv
         end do
!
         gg(n+nadd) = one
         do i=n+nadd,1,-1
            gg(i-1) = (u*gg(i)-ff(i))/i
         end do
      end if
!
      return
      end

