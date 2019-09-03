      program main
      implicit none
      integer the,t,jmax,n,s,m
      real temp
      parameter(n=2,s=3,m=8+s+n)
      character(s) stype
      character(n) beads
      character(m) filename
      character(8) fmt1,fmt2
      character(5) string
      
      
      stype='H2O'!change s if name is changed
      temp=1.0 !(temp in K)
      t=50   !beads(time slices) change n, if t changed
      jmax=40

      if (temp.lt.10.0)                           fmt1="(f5.3)"
      if ((temp.ge.10.0).and.(temp.lt.100.0))     fmt1="(f5.2)"
      if ((temp.ge.100.0).and.(temp.lt.1000.0))   fmt1="(f5.1)"
      if (temp.ge.1000.0)                         fmt1="(f5.0)"
      
      if ((t.lt.10)) fmt2="(I1)"
      if ((t.ge.10).and.(t.lt.100)) fmt2="(I2)"
      if ((t.ge.100).and.(t.lt.1000)) fmt2="(I3)"
      if ((t.ge.1000).and.(t.lt.10000)) fmt2="(I4)"
      if ((t.ge.10000).and.(t.lt.100000)) fmt2="(I5)"

       write(string,fmt1) temp
       write(beads,fmt2) t

       open(1,file='run')
       write(1,* )'!/bin/sh'
       write(1,* )'#PBS -j oe'
       write(1,*) '#PBS -m ae'
       write(1,* )'cd ${PBS_O_WORKDIR}'
       do the=0,180,1
       write(1,*)'./asymrho.x  ', string,'  ', beads,'  -1',the,the,'  0.6666525 0.2306476 0.1769383', jmax 
       enddo
       write(1,*)'./compile.x'
       filename=stype//'_T'//string//'t'//beads
       write(1,*) 'mv rho.den_eng ',filename,'.eng'
       write(1,*) 'mv rho.den_rho ',filename,'.rho'
       write(1,*) 'mv rho.den_esq ',filename,'.esq'
       close(1)
       end

