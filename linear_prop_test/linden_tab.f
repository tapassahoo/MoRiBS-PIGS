      program linden
      implicit double precision (a-h,o-z)

      character argum*30

      call getarg(1,argum)
      read(argum,*)temprt
      call getarg(2,argum)
      read(argum,*)nslice
      call getarg(3,argum)
      read(argum,*)bconst
      call getarg(4,argum)
      read(argum,*)iodevn

      cost1 = 0.0d0
      phi1 = 0.0d0
      cost3 = 0.0d0
      phi3 = 0.0d0
      call propos(cost1, phi1, cost3, phi3, temprt, nslice,
     & bconst, iodevn)

      stop 
      end

      subroutine propos(cost1, phi1, cost3, phi3, temprt, nslice,
     & bconst, iodevn)
      implicit double precision (a-h,o-z)

      parameter(maxl=500)
      parameter(taunit=1.4387752224d+00,
     +          pi=3.14159265358979323846d+00,eps=1.d-16)
      dimension pl(0:maxl)
      parameter (ncost2 = 101)
      parameter (nphi2 = 361)
      dimension distrb(ncost2*nphi2)

c ... calculate tau
c.... unit of tau is cm and bconst is cm-1 
      tau=taunit/(temprt*nslice)  

c======================================================================c
c ... judge the maximam l quantum number
      do l=0,maxl
        if(exp(-tau*bconst*l*(l+1)).lt.eps)then
          lmax=l
          goto 20
        endif
      enddo
      write(6,'(''!!! Warning: maxl is reached'')')
      lmax=maxl
   20 continue

c======================================================================c
      sint1  = sqrt(1.0-cost1*cost1)
      sint3  = sqrt(1.0-cost3*cost3)
c======================================================================c
      open(7,file='linden.out',status='unknown')

c======================================================================c

      cstep2 = 2.0/dfloat(ncost2-1)
      pstep2 = 2.0d0*pi/dfloat(nphi2-1)

      sum=0.0d0
      do ic2 = 1, ncost2
        cost2 = (ic2-1)*cstep2-1.d0
        sint2 = sqrt(1.0d0-cost2*cost2)
        do ip2 = 1,nphi2 
          phi2 = (ip2-1)*pstep2
          ii = ip2 + (ic2 - 1)*nphi2
          cost12 = sint1*sint2*cos(phi2-phi1)+cost1*cost2
          cost23 = sint2*sint3*cos(phi3-phi2)+cost2*cost3
          call exarho(cost12,lmax,maxl,pl,rho12,tau,bconst,
     &    iodevn,nslice)
          call exarho(cost23,lmax,maxl,pl,rho23,tau,bconst,
     &    iodevn,nslice)
          distrb(ii) = rho12*rho23
          sum = sum+distrb(ii)
        enddo
      enddo
      sum1 = 0.0d0
      do ii = 1, ncost2*nphi2
        distrb(ii) = distrb(ii)/sum
        write(7,'(1p,7(1x,E15.8))')distrb(ii)
        sum1 = sum1 + distrb(ii)
      enddo
      close(7,status='keep')
      write(*,*) sum1
      return
      end

c======================================================================c
      subroutine exarho(cost,lmax,maxl,pl,rho,tau,bconst,iodevn,nslice)
      implicit double precision(a-h,o-z)
      parameter(pi=3.14159265358979323846d+00,boltz=0.69503476d0)
      dimension pl(0:maxl)

      call lgnd(lmax,cost,pl)

      rho=0.d0
      do l=0,lmax
        if((mod(l,2).eq.iodevn).or.iodevn.eq.-1) then
          tmp=dfloat(2*l+1)*pl(l)
          tmp=tmp*exp(-tau*bconst*dfloat(l*(l+1)))
          rho=rho+tmp
        endif
      enddo
      rho=rho/(4.0*Pi)
      return
      end
ccccccccccccccccccccccccc     Program 5.4     cccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Please Note:                                                         c
c                                                                      c
c (1) This computer program is part of the book, "An Introduction to   c
c     Computational Physics," written by Tao Pang and published and    c
c     copyrighted by Cambridge University Press in 1997.               c
c                                                                      c
c (2) No warranties, express or implied, are made for this program.    c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE LGND(LMAX,X,P)
      implicit double precision(a-h,o-z)
C
C Subroutine to generate Legendre polynomials P_L(X)
C for L = 0,1,...,LMAX with given X.
C
      DIMENSION P(0:LMAX)
      P(0) = 1.
      P(1) = X
      DO 100 L = 1, LMAX-1
        P(L+1) = ((2.0*L+1)*X*P(L)-L*P(L-1))/(L+1)
  100 CONTINUE
      RETURN
      END
