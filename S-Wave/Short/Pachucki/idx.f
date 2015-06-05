c    idx2(n1,n2) -> N; ni >=0

      integer function idx2(n1,n2)
      IMPLICIT NONE
      integer n1,n2,nx,sum

      if ((n1<0).or.(n2<0))then
       idx2=0
       return
      endif

       nx  = n1
       sum = nx+1
       nx  = nx+n2
       sum = sum + (nx*(1 + nx))/2

       idx2 = sum

       return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c    idx3(n1,n2,n3) -> N; ni >=0

      integer function idx3(n1,n2,n3)
      IMPLICIT NONE
      integer n1,n2,n3,nx,sum

      if ((n1<0).or.(n2<0).or.(n3<0))then
       idx3=0
       return
      endif

       nx  = n1
       sum = nx+1
       nx  = nx+n2
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n3
       sum = sum + (nx*(1 + nx)*(2 + nx))/6

       idx3 = sum

       return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c    idx5(n1,n2,n3,n4,n5) -> N; ni >=0

      integer function idx5(n1,n2,n3,n4,n5)
      IMPLICIT NONE
      integer n1,n2,n3,n4,n5,nx,sum

      if ((n1<0).or.(n2<0).or.(n3<0).or.(n4<0).or.(n5<0))then
       idx5=0
       return
      endif

       nx  = n1
       sum = nx+1
       nx  = nx+n2
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n3
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n4
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n5
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120

       idx5 = sum

       return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c    idx6(n1,n2,n3,n4,n5,n6) -> N; ni >=0

      integer function idx6(n1,n2,n3,n4,n5,n6)
      USE PRMTS
      IMPLICIT NONE
      integer n1,n2,n3,n4,n5,n6,nx,sum

      if ((n1>-1).and.(n2>-1).and.(n3>-1).and.(n4>-1).and.(n5>-1).and.(n6>-1)) then
       nx  = n1
       sum = nx+1
       nx  = nx+n2
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n3
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n4
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n5
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       nx  = nx+n6
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx)*(5 + nx))/720
       idx6= sum
       return
      endif

      IF (SWITCH.EQ.0) THEN
       idx6 = 0
       RETURN
      ENDIF

      if((n1.eq.-1).and.(n2>-1).and.(n3>-1).and.(n4>-1).and.(n5>-1).and.(n6>-1)) then
       nx  = n2
       sum = nx+1
       nx  = nx+n3
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n4
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n5
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n6
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       idx6= MAXFN+sum
       return
      endif

      if((n1>-1).and.(n2.eq.-1).and.(n3>-1).and.(n4>-1).and.(n5>-1).and.(n6>-1)) then
       nx  = n1
       sum = nx+1
       nx  = nx+n3
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n4
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n5
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n6
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       idx6= MAXFN+MAXFM+sum
       return
      endif

      if((n1>-1).and.(n2>-1).and.(n3.eq.-1).and.(n4>-1).and.(n5>-1).and.(n6>-1)) then
       nx  = n1
       sum = nx+1
       nx  = nx+n2
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n4
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n5
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n6
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       idx6= MAXFN+2*MAXFM+sum
       return
      endif

      if((n1>-1).and.(n2>-1).and.(n3>-1).and.(n4.eq.-1).and.(n5>-1).and.(n6>-1)) then
       nx  = n1
       sum = nx+1
       nx  = nx+n2
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n3
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n5
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n6
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       idx6= MAXFN+3*MAXFM+sum
       return
      endif

      if((n1>-1).and.(n2>-1).and.(n3>-1).and.(n4>-1).and.(n5.eq.-1).and.(n6>-1)) then
       nx  = n1
       sum = nx+1
       nx  = nx+n2
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n3
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n4
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n6
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       idx6= MAXFN+4*MAXFM+sum
       return
      endif

      if((n1>-1).and.(n2>-1).and.(n3>-1).and.(n4>-1).and.(n5>-1).and.(n6.eq.-1)) then
       nx  = n1
       sum = nx+1
       nx  = nx+n2
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n3
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n4
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n5
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       idx6= MAXFN+5*MAXFM+sum
       return
      endif

      IF (SWITCH.EQ.1) THEN
       idx6 = 0
       RETURN
      ENDIF

      if((n1.eq.-2).and.(n2>-1).and.(n3>-1).and.(n4>-1).and.(n5>-1).and.(n6>-1)) then
       nx  = n2
       sum = nx+1
       nx  = nx+n3
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n4
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n5
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n6
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       idx6= MAXFN+6*MAXFM+sum
       return
      endif

      if((n1>-1).and.(n2.eq.-2).and.(n3>-1).and.(n4>-1).and.(n5>-1).and.(n6>-1)) then
       nx  = n1
       sum = nx+1
       nx  = nx+n3
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n4
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n5
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n6
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       idx6= MAXFN+7*MAXFM+sum
       return
      endif

      if((n1>-1).and.(n2>-1).and.(n3.eq.-2).and.(n4>-1).and.(n5>-1).and.(n6>-1)) then
       nx  = n1
       sum = nx+1
       nx  = nx+n2
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n4
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n5
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n6
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       idx6= MAXFN+8*MAXFM+sum
       return
      endif

      if((n1>-1).and.(n2>-1).and.(n3>-1).and.(n4.eq.-2).and.(n5>-1).and.(n6>-1)) then
       nx  = n1
       sum = nx+1
       nx  = nx+n2
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n3
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n5
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n6
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       idx6= MAXFN+9*MAXFM+sum
       return
      endif

      if((n1>-1).and.(n2>-1).and.(n3>-1).and.(n4>-1).and.(n5.eq.-2).and.(n6>-1)) then
       nx  = n1
       sum = nx+1
       nx  = nx+n2
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n3
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n4
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n6
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       idx6= MAXFN+10*MAXFM+sum
       return
      endif

      if((n1>-1).and.(n2>-1).and.(n3>-1).and.(n4>-1).and.(n5>-1).and.(n6.eq.-2)) then
       nx  = n1
       sum = nx+1
       nx  = nx+n2
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n3
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       nx  = nx+n4
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx))/24
       nx  = nx+n5
       sum = sum + (nx*(1 + nx)*(2 + nx)*(3 + nx)*(4 + nx))/120
       idx6= MAXFN+11*MAXFM+sum
       return
      endif

      idx6 = 0
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccc

!      real (16) function delta(n)
!      IMPLICIT NONE
!       integer n
!       if (n==0) then
!        delta=1.0q0
!       else
!        delta=0.0q0
!       endif
!       return
!      end

cccccccccccccccccccccccccccccccccccccccccccccccccccc
