      implicit real*8(a-h,o-z)
      parameter (nn=2000)
      dimension pot(nn)
      dimension energy(nn)
      open(1,file='analyze_input')
      open(2,file='average_analyze_input')
      param1 = param2
      sum1=0.0d0
      sum1var=0.0d0
      sum2=0.0d0
      sum2var=0.0d0
      do i=1,nn
        read(1,*)iblock,ke,pot(i),energy(i),aa,bb,cc,dd,ee,ff
        sum1=sum1+pot(i)
        sum1var=sum1var+pot(i)*pot(i)
        sum2=sum2+energy(i)
        sum2var=sum2var+energy(i)*energy(i)
      enddo
      dev1=sqrt(sum1var-sum1*sum1)
      dev2=sqrt(sum2var-sum2*sum2)
      write(2,*)nn,param1,sum1/float(nn),sum2/float(nn),dev1,dev2
      stop
      end
