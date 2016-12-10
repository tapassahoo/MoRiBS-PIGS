      implicit real*8(a-h,o-z)
      parameter (nn=5000)
      dimension energy(nn)
      open(1,file='analyze_input')
      open(2,file='average_analyze_input')
      sum=0.0d0
      do i=1,nn
        read(1,*)iblock,aa,energy(i),bb,cc,dd,ee,ff,gg,aa1,bb1
        sum=sum+energy(i)
      enddo
      write(2,*)nn,sum/float(nn)
      stop
      end
