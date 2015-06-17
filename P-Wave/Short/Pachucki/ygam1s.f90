      subroutine ygam1s(alf,bet,gam,yg,n)
!----------------------------------------------------------------------
!   changed on February 11, 2008 by MP, KP.
!
!   Purpose           - calculates Gamma       (alpha,beta,gamma)
!                                       -1,-1,n
!
!                                  Gamma       (alpha,beta,gamma)
!                                       -2,-1,n
!
!   Remarks:         1: see R.A.Sack, C.C.J.Roothaan, and W.Kolos,
!                         J. Math. Phys., v.8 (1967) 1093, for notation
!
!----------------------------------------------------------------------
      implicit          real (16) (a-h,o-z)
      parameter         (lenl=200,lenl2=3*lenl)
      real (16), parameter :: zero = 0.0q0
      real (16), parameter :: one  = 1.0q0
      real (16), parameter :: threshold  = 0.05q0
      real (16)         yg(0:lenl,2)
      real (16), dimension(0:lenl2) :: cc,dd,ee,ff,gg,hh,ii1,ii2
!
      epsm = epsilon(one)
      agam = abs(gam)
      if (agam >= threshold) then
         nl = n
      else if (agam == zero) then
         nl = n+1
      else
         nadd = int(log(epsm)/log(agam))+5
         nl = n+nadd
         if (nl > lenl2) then
            write (*,*) 'Error!!! Insufficient array length in YGAM1!'
            stop
         end if
      end if
!
      ainv = one/(bet+gam)
      binv = one/(alf+gam)
      ginv = one/(alf+bet)
!
!                       b/[(b+c)(b-c)]*ln[(a+b)/(a+c)]
!
      if (bet > gam) then
         u = (bet-gam)/(alf+gam)
         ab = alf+bet
!
         call dplog(u,nl-1,gg)
!
         cc(0) = gg(0)
         cc(1) = gg(1)*binv
         temp = binv
         do i=2,nl-1
            temp = temp*binv
            cc(i) = gg(i)*temp
         end do
         ee(0) = bet*cc(0)
         do j=1,nl-1
            do i=0,nl-j-1
               cc(i) = binv*((i+j-1)*cc(i)-ab*cc(i+1))
            end do
            ee(j) = bet*cc(0)
         end do
         ee(0) = ainv*ee(0)
         do j=1,nl-1
            ee(j) = ainv*(j*ee(j-1)+ee(j))
         end do
         ee(0) = binv*ee(0)
         do j=1,nl-1
            ee(j) = binv*(j*ee(j-1)+ee(j))
         end do
      else
         u = (gam-bet)/(alf+bet)
         ag = alf+gam
!
         call dplog(u,nl-1,gg)
!
         cc(0) = gg(0)
         cc(1) = gg(1)*ginv
         ee(0) = bet*cc(0)
         ee(1) = bet*cc(1)
         temp = ginv
         do i=2,nl-1
            temp = temp*ginv
            cc(i) = gg(i)*temp
            ee(i) = bet*cc(i)
         end do
         ee(0) = ainv*ee(0)
         do j=1,nl-1
            ee(j) = ainv*(j*ee(j-1)+ee(j))
         end do
         do j=0,nl-1
            ee(j) = ginv*ee(j)
         end do
      end if
!
!                       a/[(a+c)(a-c)]*ln[(a+b)/(b+c)]
!                         and (a+b)/[(b+c)(a-c)]*ln[(a+b)/(b+c)]
!
      if (alf > gam) then
         u = (alf-gam)/(bet+gam)
         ab = alf+bet
!
         call dplog(u,nl-1,gg)
!
         cc(0) = gg(0)
         cc(1) = gg(1)*ainv
         temp = ainv
         do i=2,nl-1
            temp = temp*ainv
            cc(i) = gg(i)*temp
         end do
         dd(0) = alf*cc(0)
         hh(0) = ab*cc(0)
         do j=1,nl-1
            do i=0,nl-j-1
               cc(i) = ainv*((i+j-1)*cc(i)-ab*cc(i+1))
            end do
            dd(j) = alf*cc(0)
            hh(j) = ab*cc(0)
         end do
         dd(0) = binv*dd(0)
         hh(0) = ainv*hh(0)
         do j=1,nl-1
            dd(j) = binv*(j*dd(j-1)+dd(j))
            hh(j) = ainv*(j*hh(j-1)+hh(j))
         end do
         dd(0) = ainv*dd(0)
         hh(0) = ainv*hh(0)
         do j=1,nl-1
            dd(j) = ainv*(j*dd(j-1)+dd(j))
            hh(j) = ainv*(j*hh(j-1)+hh(j))
         end do
      else
         u = (gam-alf)/(alf+bet)
         bg = bet+gam
!
         call dplog(u,nl-1,gg)
!
         cc(0) = gg(0)
         cc(1) = gg(1)*ginv
         dd(0) = alf*cc(0)
         dd(1) = alf*cc(1)
         temp = ginv
         do i=2,nl-1
            temp = temp*ginv
            cc(i) = gg(i)*temp
            dd(i) = alf*cc(i)
         end do
         dd(0) = binv*dd(0)
         hh(0) = ainv*cc(0)
         do j=1,nl-1
            dd(j) = binv*(j*dd(j-1)+dd(j))
            hh(j) = ainv*(j*hh(j-1)+cc(j))
         end do
         do j=0,nl-1
            dd(j) = ginv*dd(j)
         end do
      end if
!
!                       Gamma       (alpha,beta,gamma)
!                            -1,-1,n
!
      if (agam >= threshold) then
         gl = log((alf+gam)/(bet+gam))
         dilog1 = dilog((alf+bet)/(alf+gam))
         dilog2 = dilog((alf+bet)/(bet+gam))
         yg(0,1) = (gl**2+2*(dilog1+dilog2)+pi**2/3)/(4*gam)
         gaminv = one/gam
         do i=1,n
            yg(i,1) = gaminv*(i*yg(i-1,1)-dd(i-1)-ee(i-1))
         end do
      else
         cc(nl) = one
         do i=nl,1,-1
            cc(i-1) = (gam*cc(i)+dd(i-1)+ee(i-1))/i
         end do
         do i=0,n
            yg(i,1) = cc(i)
         end do
      end if
!
!                       Gamma       (alpha,beta,gamma)
!                            -2,-1,n
!
      clogbg = log(bet+gam)
      ff(0) = clogbg**2/2
      ff(1) = -clogbg*ainv
      temp = ainv
      ff(2) = ainv*(ff(1)+temp)
      do i=3,n
         temp = temp*ainv*(i-2)
         ff(i) = ainv*((i-1)*ff(i-1)+temp)
      end do
!
      if (agam < threshold) then
         dilog2 = dilog((alf+bet)/(bet+gam))
      end if
      ag = alf+gam
      gg(0) = ag*yg(0,1)
      yg(0,2) = ff(0)+dilog2-gg(0)
      do i=1,n
         gg(i)= -i*yg(i-1,1)+ag*yg(i,1)
         yg(i,2) = ff(i)-hh(i-1)-gg(i)
      end do


!     added 06.01.2008 by MP                                                   
!     derivatives over -c of -log(b+c)+pi^2/6+1                                
!                                                                              

      ii1(0) = -log(bet+gam) + Pi**2/6 + 1
      ii1(1)=ainv

      do i=2,n
       ii1(i)=(i-1)*ii1(i-1)*ainv
      enddo

      yg(0,2) = yg(0,2)+ii1(0)

      do i=1,n
       yg(i,2) = yg(i,2)+ii1(i)
      enddo

      return
      end
