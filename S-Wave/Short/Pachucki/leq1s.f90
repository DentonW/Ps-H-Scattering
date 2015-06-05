      subroutine leq1s (a,n,b,m,ib,ijob,wa,ier)
!-------------------------------------------------------------
!
!   latest revision     - april 19, 1993
!
!   purpose             - linear equation solution -
!                           indefinite matrix - symmetric
!                           storage mode - space economizer
!                           solution
!
!   usage               - call leq1s (a,n,b,m,ib,ijob,wa,ier)
!
!   arguments    a      - input/output vector of dimension
!                           n*(n+1)/2. see parameter ijob.
!                n      - order of matrix a and the number of
!                           rows in b. (input)
!                b      - input/output matrix of dimension n
!                           by m. on input, b contains the m
!                           right-hand sides of the equation
!                           ax = b. on output, the solution
!                           matrix x replaces b. if ijob = 1,
!                           b is not used.
!                m      - number of right hand sides (columns
!                           in b). (input)
!                ib     - column dimension of b exactly as
!                           specified in the dimension
!                           statement in the calling program.
!                           (input)
!                ijob   - input option parameter. ijob = i
!                           implies when
!                           i = 0, factor the matrix a and
!                             solve the equation ax = b. on
!                             input, a contains the coefficient
!                             matrix of the equation ax = b,
!                             where a is assumed to be an
!                             n by n symmetric matrix. a is
!                             stored in symmetric storage mode
!                             and therefore has dimension
!                             n*(n+1)/2. on output, a is
!                             replaced by its factorized form.
!                           i = 1, factor the matrix a.
!                             a contains the same input/output
!                             information as if ijob = 0. b is
!                             not used.
!                           i = 2, solve the equation ax = b.
!                             this option implies that leq1s
!                             has already been called using
!                             ijob = 0 or 1 so that the matrix
!                             a has already been factored.
!                             in this case, output matrix a
!                             must have been saved for reuse
!                             in the call to leq1s.
!                wa     - work area of length n+1,
!                           on output wa(n+1) contains the
!                           number of eigenvalues below zero
!                ier    - error parameter. (output)
!                         warning error (with fix)
!                           ier = 65 indicates that matrix a
!                             is algorithmically singular.
!                         terminal error
!                           ier = 129 indicates that matrix a
!                             is exactly singular on one step
!                             of procedure.
!
!-------------------------------------------------------------
      implicit real (16) (a-h,o-z)
      dimension          a(*),b(ib,*),wa(*)
      real (16)     :: zero=0.0q0
      real (16)     :: alpha=0.6403882q0
      real (16)     :: sixtn=16.0q0
      save               zero,alpha,sixtn
!
!                                  first executable statement
!                                    initialize ier
!
      ier = 0
!
!   ****  first stage  ****        decompose a into the
!                                    product m*d*m-transpose
!                                    where m is unit lower
!                                    triangular and d is block
!                                    diagonal with blocks of
!                                    order 1 or 2.
!                                    m and d are written over a.
!
      if (n.le.0) goto 9005
      if (ijob.eq.2) goto 80
!
!                                  calculate equilibration
!                                    factors
!
      rn = sixtn*n
      wa(n+1) = 0
      ll = 1
      do 10 i=1,n
         wa(i) = zero
         l = ll
         do 5 j=1,n
            temp = abs(a(l))
            if (temp.gt.wa(i)) wa(i) = temp
            if (j.lt.i) then
               l = l+1
            else
               l = l+j
            end if
    5    continue
         wa(i) = wa(i)*rn
         ll = ll+i
   10 continue
      i = n
      if (i.eq.1) goto 75
!
!                                  ****  begin of the loop  ****
!
   15 im1 = i-1
      im2 = i-2
      ix = i*im1/2
      ixi = ix+i
!
!                                  calculate maximum off
!                                    diagonal element in row i
!
      aii = abs(a(ixi))
      save = zero
      do 20 l=1,im1
         temp = abs(a(ix+l))
         if (temp.gt.save) then
            save = temp
            j = l
         end if
   20 continue
      ichang = i
      if (aii.ge.alpha*save) goto 35
!
!                                  calculate maximum off
!                                    diagonal element in row j
!
      sigma = save
      jx = j*(j-1)/2
      jxj = jx+j
      jxt = jx+1
      ajj = abs(a(jxj))
      do 25 l=1,im1
         temp = abs(a(jxt))
         if (temp.gt.sigma.and.l.ne.j) sigma = temp
         if (l.lt.j) then
            jxt = jxt+1
         else
            jxt = jxt+l
         end if
   25 continue
!
!                                  choose a strategy for the
!                                    pivoting
!
      if (aii*sigma.ge.alpha*save*save) goto 35
!
      if (ajj.lt.alpha*sigma) goto 50
!
!                                  interchange rows i and j
!
      jxt = jx+1
      do 30 l=1,im1
         temp = a(ix+l)
         a(ix+l) = a(jxt)
         a(jxt) = temp
         if (l.lt.j) then
            jxt = jxt+1
         else
            jxt = jxt+l
         end if
   30 continue
      temp = a(jxj)
      a(jxj) = a(ixi)
      a(ixi) = a(ix+j)
      a(ix+j) = temp
      temp = wa(i)
      wa(i) = wa(j)
      wa(j) = temp
      aii = ajj
      ichang = j
!
!                                  we use a 1 by 1 pivot
!
   35 if (wa(i)+aii.le.wa(i)) ier = 65
      if (aii.eq.zero) goto 9005
      wa(i) = ichang
!
      aii = a(ixi)
      if (aii.lt.zero) wa(n+1) = wa(n+1)+1
      ixt = ix+1
      do 45 j=im1,1,-1
         save = -a(ix+j)/aii
         ixl = ixt-1
         ixt = ixt-j
         ish = ix-ixt+1
         mm = mod(j,5)
         if (mm.ne.0) then
            do 39 k=ixt,ixt+mm-1
               a(k) = a(k)+save*a(ish+k)
   39       continue
         end if
         ixf = ixt+mm
         do 40 k=ixf,ixl,5
            a(k) = a(k)+save*a(ish+k)
            a(k+1) = a(k+1)+save*a(ish+k+1)
            a(k+2) = a(k+2)+save*a(ish+k+2)
            a(k+3) = a(k+3)+save*a(ish+k+3)
            a(k+4) = a(k+4)+save*a(ish+k+4)
   40    continue
         a(ix+j) = save
   45 continue
      i = im1
      goto 75
!
!                                  we use a 2 by 2 pivot
!
   50 ix1 = ix-im1
      wa(n+1) = wa(n+1)+1
      if (j.eq.im1) goto 60
!
!                                  interchange rows i-1 and j
!
      jxt = jx+1
      do 55 l=1,im2
         temp = a(ix1+l)
         a(ix1+l) = a(jxt)
         a(jxt) = temp
         if (l.lt.j) then
            jxt = jxt+1
         else
            jxt = jxt+l
         end if
   55 continue
      temp = a(jxj)
      a(jxj) = a(ix)
      a(ix) = a(ix1+j)
      a(ix1+j) = temp
      temp = a(ix+j)
      a(ix+j) = a(ix+im1)
      a(ix+im1) = temp
      temp = wa(im1)
      wa(im1) = wa(j)
      wa(j) = temp
!
!                                  det must be negative
!
   60 det = a(ixi)*a(ix)-a(ixi-1)**2
!
      temp = max(wa(i),wa(im1))
      if (temp+abs(det).le.temp) ier = 65
      if (det.ge.zero) goto 9005
      wa(i) = j
      wa(im1) = det
!
      aim1i = a(ixi-1)/det
      aii = a(ixi)/det
      aim1 = a(ix)/det
      ixt = ix1
      do 70 j=im2,1,-1
         save = aim1i*a(ix1+j)-aim1*a(ix+j)
         temp = aim1i*a(ix+j)-aii*a(ix1+j)
         ixt = ixt-j
         do 65 k=1,j
            a(ixt+k) = a(ixt+k)+a(ix+k)*save+a(ix1+k)*temp
   65    continue
         a(ix+j) = save
         a(ix1+j) = temp
   70 continue
      i = im2
!
!                                  *** end of the loop ***
!
   75 if (i.gt.1) goto 15
      if (i.eq.1) then
         if (wa(1)+abs(a(1)).le.wa(1)) ier = 65
         if (a(1).eq.zero) goto 9005
         if (a(1).lt.zero) wa(n+1) = wa(n+1)+1
         wa(1) = 1
      end if
!
!   ****  second stage  ****       solve m*d*mt*x = b where
!                                    mt = m-transpose
!
   80 if (ijob.eq.1) go to 9000
!
      do 120 jc=1,m
!
!                                  solve m*d*y = b and store
!                                    y in b
!
         i = n
   85    if (i.le.1) goto 100
         im1 = i-1
         ix = i*im1/2
         ixi = ix+i
         ichang = wa(i)
         if (ichang.le.0) then
            write (*,'(2x,a)') 'error in routine leq1s |',              &
     &                         'array wa is wrong |'
            return
         end if
         save = b(ichang,jc)
         if (wa(im1).gt.zero) then
!
!                                  we use a 1 by 1 pivot
!
            b(ichang,jc) = b(i,jc)
            b(i,jc) = save/a(ixi)
            do 90 j=1,im1
               b(j,jc) = b(j,jc)+a(ix+j)*save
   90       continue
            i = im1
         else
!
!                                  we use a 2 by 2 pivot
!
            temp = b(i,jc)
            b(ichang,jc) = b(im1,jc)
            det = wa(im1)
            ix1 = ix-im1
            im2 = i-2
            b(i,jc) = (temp*a(ix)-save*a(ix+im1))/det
            b(im1,jc) = (save*a(ixi)-temp*a(ix+im1))/det
            do 95 j=1,im2
               b(j,jc) = b(j,jc)+a(ix+j)*temp+a(ix1+j)*save
   95       continue
            i = im2
         end if
         goto 85
  100    if (i.eq.1) then
            b(1,jc) = b(1,jc)/a(1)
            i = 2
         else
            i = 3
         end if
!
!                                  solve mt*x = y and store
!                                    x in b
!
  105    if (i.gt.n) goto 120
            ii = i
            if (wa(i).le.zero) ii = i+1
            im1 = i-1
            ix = i*im1/2
            do 115 k=i,ii
               save = b(k,jc)
               do 110 j=1,im1
                  save = save+a(ix+j)*b(j,jc)
  110          continue
               b(k,jc) = save
               ix = ix+k
  115       continue
            ichang = wa(ii)
            save = b(i,jc)
            b(i,jc) = b(ichang,jc)
            b(ichang,jc) = save
            i = ii+1
         goto 105
  120 continue
!
!                                  exit
!
 9000 if (ier.eq.0) return
      if (ier.eq.65) write (*,9065)
 9065 format(/5x,'warning error  ier=65  in  *** leq1s ***'/            &
     &       7x,'matrix a is algorithmically singular.'/)
      return
 9005 ier = 129
      write (*,9129)
 9129 format(/5x,'terminal error  ier=129  in  *** leq1s ***'/          &
     &       7x,'matrix a is exactly singular.'/)
      return
      end
