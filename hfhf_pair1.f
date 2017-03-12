      subroutine cluster1(Rpt,Eulang_1,Eulang_2,E_12)
      implicit double precision(a-h,o-z)
      parameter (ndim = 3)
      dimension Eulang_1(ndim),Eulang_2(ndim)
      dm = 1.86d0
      E_12_au = dm*dm*(Eulang_1(0)*Eulang_2(0) + Eulang_1(1)*Eulang_2(1)
     &     - Eulang_1(2)*Eulang_2(2))/Rpt

      E_12_au=v3_elec !in Hartree
      E_12=E_12_au*autoK !in Kelvin
      E_12_kcal=E_12_au*autokcal !in Kcal/mol  
      E_12_kJ=E_12_au*autokJ !in kJ/mol
      E_12_cm=E_12_au*autocm !in cm-1      
 
c      write(*,*) E_12
c      write(44,72) RO_1_sf,RH1_1_sf,RH2_1_sf      
c      write(44,73) RO_2_sf,RH1_2_sf,RH2_2_sf
c      write(44,71) com_1(1:3),com_2(1:3),R12(1:3),R
c      write(44,74) E_12,E_12_cm
c      write(44,*)  "            "
   71 format("R1=(",3F8.4,")   R2=(",3F8.4,")   R12=(",3F8.4,")  R=",
     1      1F8.4, "   Angstrom")
   72 format("1st mol  O=(",3F12.8,")  H=(",3F12.8,")  H=(",3F12.8,")")
   73 format("2nd mol  O=(",3F12.8,")  H=(",3F12.8,")  H=(",3F12.8,")")
   74 format("EKel = ",1F10.4,"  Ecm-1 = ",1F10.4)
      return
      end
