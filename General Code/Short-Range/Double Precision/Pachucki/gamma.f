       MODULE PRMTS
       INTEGER MAXF,MAXG,MAXFM,MAXFN,SHELL,SWITCH
       PARAMETER(SHELL=30)  ! This gives us up to omega = 10.
!      SHELL=SHELL+1 if grad. min. or hfs
!      SHELL=SHELL+2*L, L - angular momentum
!      max(SHELL)=2*6+6 = 18 in quad
!      max(SHELL)=2*9+6 = 24 in sixr
!      max(SHELL)=2*12+6= 30 in octr

       PARAMETER(MAXG=(7+SHELL)*(8+SHELL)*(9+SHELL)/6-1)
       PARAMETER(MAXFN=((1+SHELL)*(2+SHELL)*(3+SHELL)*(4+SHELL)*(5+SHELL)*(6+SHELL))/720)
!      PARAMETER(MAXFM=((1+SHELL)*(2+SHELL)*(3+SHELL)*(4+SHELL)*(5+SHELL))/120)
       PARAMETER(MAXF = MAXFN) ! if SWITCH=0
!      PARAMETER(MAXF = MAXFN+6*MAXFM)  ! if SWITCH=1
!      PARAMETER(MAXF = MAXFN+12*MAXFM) ! if SWITCH=2
       END MODULE PRMTS


c gamma(l,m,n); l,m,n >=-2; l+m+n<=sm
      subroutine p_g(tg,a,b,c,sm)
       USE PRMTS
       IMPLICIT NONE

       integer sm
     
       ! macierz wspolczynnikow gamma w artykule
       real*16 tg(0:MAXG)

       ! macierze B oraz A w artykule
       real*16 tb(0:MAXG)
       real*16  yg1(0:200,0:200,2)
       real*16  yg2(0:200,2)
       real*16 a,b,c,p
      
       ! aktualna powloka

       integer sc  

       !  identyfikator funkcji idx zwracajacej indeks w tablicy tg,tb

       integer idx
       integer l,m,n
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	!write (*,*) "Test 4"

       call p_tb(sm,tb,a,b,c)
      
       do sc=0,sm

        if (sc==0) then 
         tg(idx(0,0,0)) = 1/(a+b)/(b+c)/(c+a)

         goto 10
        endif


        do l=0, sc      
         do m=0, sc-l             

          n=sc-(l+m)          

          tg(idx(l,m,n))=1/(a+b)*tb(idx(l,m,n))

          if (l>0) then 
           tg(idx(l,m,n))=tg(idx(l,m,n))
     -     +1/(a+b)*l*tg(idx(l-1,m,n))
          endif

          if (m>0) then 
           tg(idx(l,m,n))=tg(idx(l,m,n))
     -     +1/(a+b)*m*tg(idx(l,m-1,n))
          endif  

          enddo
         enddo     
 10      continue
        enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! to powinno byc z kodu Korobova 

        call ygam(a,b,c,yg1,sm+2,sm+2,2)
        do m=0,sm+1
         do n=0,sm+1-m
          tg(idx(-1,m,n))=yg1(m,n,1)
         enddo
        enddo
        do m=0,sm+2
         do n=0,sm+2-m
          tg(idx(-2,m,n))=yg1(m,n,2)
         enddo
        enddo

        call ygam(b,a,c,yg1,sm+2,sm+2,2)
        do m=0,sm+1
         do n=0,sm+1-m
          tg(idx(m,-1,n))=yg1(m,n,1)
         enddo
        enddo
        do m=0,sm+2
         do n=0,sm+2-m
          tg(idx(m,-2,n))=yg1(m,n,2)
         enddo
        enddo

        call ygam(c,b,a,yg1,sm+2,sm+2,2)
        do m=0,sm+1
         do n=0,sm+1-m
          tg(idx(n,m,-1))=yg1(m,n,1)
         enddo
        enddo     
        do m=0,sm+2
         do n=0,sm+2-m
          tg(idx(n,m,-2))=yg1(m,n,2)
         enddo
        enddo     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        p=b;b=c;c=p

        call ygam1s(a,b,c,yg2,sm+3)
        do n=0,sm+2
         tg(idx(-1,-1,n))= yg2(n,1)
        enddo
        do n=0,sm+3
         tg(idx(-2,-1,n))= yg2(n,2)
        enddo
        call ygam1s(b,a,c,yg2,sm+3)
        do n=0,sm+3
         tg(idx(-1,-2,n))= yg2(n,2)
        enddo

        call ygam1s(c,b,a,yg2,sm+3)
        do n=0,sm+2
         tg(idx(n,-1,-1))= yg2(n,1) 
        enddo
        do n=0,sm+3
         tg(idx(n,-1,-2))= yg2(n,2) 
        enddo

        call ygam1s(b,c,a,yg2,sm+3)
        do n=0,sm+3
         tg(idx(n,-2,-1))= yg2(n,2) 
        enddo

        call ygam1s(a,c,b,yg2,sm+3)
        do n=0,sm+2
         tg(idx(-1,n,-1))= yg2(n,1)
        enddo
        do n=0,sm+3
         tg(idx(-2,n,-1))= yg2(n,2)
        enddo

        call ygam1s(c,a,b,yg2,sm+3)
        do n=0,sm+3
         tg(idx(-1,n,-2))= yg2(n,2)
        enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

 20     continue
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! wyliczanie elementow beta

      subroutine p_tb(sm,tb,a,b,c)
       use prmts
       IMPLICIT NONE
       real*16 a,b,c
       real*16 tb(0:MAXG)

       integer sm,sc
       integer l,m,n,idx
       real*16 ta(0:sm+3)

       do sc=0,sm

        if (sc==0) then 
         tb(idx(0,0,0))=1/((a+c)*(b+c))         
         ta(0) = 1/(b+c)
         goto 20
        endif 

  
        ta(sc) = ta(sc-1)*sc/(b+c)

        do l=0,sc
         do m=0,sc-l  
          n=sc-(l+m)
          if (n>0) then          
           tb(idx(l,m,n))= 1/(a+c)*n*tb(idx(l,m,n-1))
          else
           tb(idx(l,m,n))=0 
          endif
 
          if (l==0) then
           if (n==0) then
            tb(idx(l,m,n))=tb(idx(l,m,n))
     -      +1/(b+c)*m*tb(idx(l,m-1,n))
           else          
            tb(idx(l,m,n))=tb(idx(l,m,n))
     -      +1/(a+c)*ta(m+n)
           endif
          else
           tb(idx(l,m,n))=tb(idx(l,m,n))+1/(a+c)*l*tb(idx(l-1,m,n))
          endif

c          write(*,*) 'tb(lmn):',l,' ',m,' ',n,' ',tb(idx(l,m,n))

         enddo
        enddo
 20     continue
       enddo 


       return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    idx(n1,n2,n3) -> N; ni >=-2

      integer function idx(n1,n2,n3)
      USE PRMTS
      IMPLICIT NONE
      integer n1,n2,n3,nx,sum

      if ((n1<-2).or.(n2<-2).or.(n3<-2))then
       idx=0
       return
      endif

       nx  = n1+2
       sum = nx+1
       nx  = nx+n2+2
       sum = sum + (nx*(1 + nx))/2
       nx  = nx+n3+2
       sum = sum + (nx*(1 + nx)*(2 + nx))/6
       idx = sum

       if (idx>MAXG) then
        write(*,*) 'ZAKRES MAXG:',MAXG,idx,n1,n2,n3
        STOP 'too small MAXG'
       endif

       return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                         
