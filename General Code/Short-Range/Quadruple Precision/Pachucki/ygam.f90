      subroutine ygam(alf,bet,gam,yg,m,n,nintg)
!----------------------------------------------------------------------
!
!   Version           - July 10, 2000, changed on February 11, 2008 MP,KP.
!
!   Purpose           - calculates Gamma      (alpha,beta,gamma)
!                                       -1,m,n
!
!                                  Gamma      (alpha,beta,gamma)
!                                       -2,m,n
!
!   Arguments         - parameter NINTG (1-2) defines a number of integrals
!                         to be evaluated starting from Gamma
!                                                            -1,m,n
!
!   Remarks:          - see R.A.Sack, C.C.J.Roothaan, and W.Kolos,
!                         J. Math. Phys., v.8 (1967) 1093, for notation
!
!----------------------------------------------------------------------
      implicit          real (16) (a-h,o-z)
      parameter         (lenl=200)
      real (16), parameter ::  one  = 1.0q0
      real (16), parameter ::  half = 0.5q0
      dimension         yg(0:lenl,0:lenl,2)
      dimension         cc(0:lenl,0:lenl),dd(0:lenl,0:lenl),ee(0:2*lenl)
      dimension         ygt(0:lenl,0:lenl,10)
      dimension         gg(0:3*lenl)

      if (bet >= gam) then
         u = (bet-gam)/(alf+gam)
         bt = bet
         gm = gam
         ainv = one/(bet+gam)
         ag = alf+gam
         binv = one/ag
         ab = alf+bet
         ginv = one/ab
         mm = m
         nn = n
         mn = m+n
         ind = 0
      else
         u = (gam-bet)/(alf+bet)
         bt = gam
         gm = bet
         ainv = one/(bet+gam)
         ag = alf+bet
         binv = one/ag
         ab = alf+gam
         ginv = one/ab
         mm = n
         nn = m
         mn = m+n
         ind = 1
      end if
!
!                     dd(i,j) - derivatives of [ln(a+b)-ln(a+g)]/(a-g)
!
      call dplog(u,mn,gg)
      cc(0,0) = gg(0)
      cc(1,0) = gg(1)*binv
      temp = binv
      do i=2,mn
         temp = temp*binv
         cc(i,0) = gg(i)*temp
      end do
      do j=1,nn
         do i=0,mn-j
            cc(i,j) = binv*((i+j-1)*cc(i,j-1)-ab*cc(i+1,j-1))
         end do
      end do
      dd(0,0) = binv*cc(0,0)
      do i=1,mm
         dd(i,0) = binv*cc(i,0)
      end do
      do j=1,nn
         dd(0,j) = binv*(j*dd(0,j-1)+cc(0,j))
         do i=1,mm
            dd(i,j) = binv*(j*dd(i,j-1)+cc(i,j))
         end do
      end do
!
!                       Gamma
!                            -1
!
      yg(0,0,1) = ainv*dd(0,0)
      do i=1,mm
         yg(i,0,1) = ainv*(i*yg(i-1,0,1)+dd(i,0))
      end do
      do j=1,nn
         yg(0,j,1) = ainv*(j*yg(0,j-1,1)+dd(0,j))
         do i=1,mm
            yg(i,j,1) = ainv*(i*yg(i-1,j,1)+j*yg(i,j-1,1)+dd(i,j))
         end do
      end do
!
!                       Gamma
!                            -2
!
      if (nintg < 2) goto 10
      clogbg = log(ab*ag)
      ygt(0,0,1) = ainv*clogbg
      yg(0,0,2) = -alf*yg(0,0,1)-half*(dd(0,0)+ygt(0,0,1))
      dlog = -ginv
      do i=1,mm
         ygt(i,0,1) = ainv*(i*ygt(i-1,0,1)+dlog)
         yg(i,0,2) = -alf*yg(i,0,1)-half*(dd(i,0)+ygt(i,0,1))
         dlog = dlog*i*ginv
      end do
      dlog = -binv
      do j=1,nn
         ygt(0,j,1) = ainv*(j*ygt(0,j-1,1)+dlog)
         yg(0,j,2) = -alf*yg(0,j,1)-half*(dd(0,j)+ygt(0,j,1))
         dlog = dlog*j*binv
         do i=1,mm
            ygt(i,j,1) = ainv*(i*ygt(i-1,j,1)+j*ygt(i,j-1,1))
            yg(i,j,2) = -alf*yg(i,j,1)-half*(dd(i,j)+ygt(i,j,1))
         end do
      end do

!     added 06.01.2008 by MP                                                   
!     add derivatives over -b,-c of 1/(b+c)                                    
!                                                                              
      temp = ainv

!     ee(i) = i!/(b+c)^(i+1)                                                   

      ee(0)=ainv
      do i=1,nn+mm
       ee(i)=i*ee(i-1)*ainv
      enddo

      do i=0,mm
       do j=0,nn
        yg(i,j,2) = yg(i,j,2) + ee(i+j)
       enddo
      enddo

!
!                    backward transposition if necessary
!
   10 if (ind == 1) then
         do k=1,nintg
            do i=0,min(m,n)
               do j=0,i-1
                  temp = yg(i,j,k)
                  yg(i,j,k) = yg(j,i,k)
                  yg(j,i,k) = temp
               end do
            end do
            do j=m+1,n
               do i=0,m
                  yg(i,j,k) = yg(j,i,k)
               end do
            end do
            do i=n+1,m
               do j=0,n
                  yg(i,j,k) = yg(j,i,k)
               end do
            end do
         end do
      end if
!
      return
      end

