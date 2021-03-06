import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
from numpy import *
import math
import inputFile
import mypkg.pkgAsymrho.support_asymrho as asym
import mypkg.pkgSymrho.support_symrho as sym

def error_message(number):
        if(number<=1):
                print("Warning!!!")  
                print("Increase the number of blocks in your computation!!!")
                exit(0)

def errorpropagation(mean, data):
        ndim   = len(data)
        error = np.std(data,ddof=0)/sqrt(ndim)
        return error
        
def maxError_byBining(mean, data, workingNdim):
        error_message(workingNdim)
        error    = np.zeros(workingNdim)
        i        = 0
        error[0] = errorpropagation(mean, data)
        for i in range(1,workingNdim):
                ndim         = int(len(data)/2)
                data1        = np.zeros(ndim)

                for j in range(ndim):
                        data1[j] = 0.5*(data[2*j]+data[2*j+1])

                data         = data1
                error[i]     = errorpropagation(mean,data)
        return np.max(error)

def makeexecutionfile(src_dir,TypeCal,ENT_TYPE, source_dir_exe):
        execution_file_dir  = source_dir_exe
        os.chdir(execution_file_dir)
        #call(["cp", "Makefile-gnu-Intel", "Makefile"])
        call(["cp", "Makefile-Copy", "Makefile"])
        call(["make", "clean"])
        call(["make"])
        print("")
        print("")
        print(".............................Compilation done!")
        print("")
        print("")
        os.chdir(src_dir)

def compile_rotmat(source_dir_exe, input_dir):
        path_enter_linden = source_dir_exe+"linear_prop/"
        os.chdir(path_enter_linden)
        call(["make", "clean"])
        call(["make"])
        path_exit_linden  = input_dir
        os.chdir(path_exit_linden)

def compile_cagepot(source_dir_exe, input_dir):
        path_enter_cagepot = source_dir_exe+"tabulated_potential/"
        os.chdir(path_enter_cagepot)
        call(["make", "clean"])
        call(["make"])
        path_exit_cagepot  = source_dir_exe+input_dir
        os.chdir(path_exit_cagepot)

def jackknife(mean,data):
        ai            = [((np.sum(data) - data[j])/(len(data) - 1.0)) for j in range(len(data))]
        deviation     = ai - mean
        devsq         = np.multiply(deviation,deviation)
        error         = sqrt(np.sum(devsq)*(len(data)-1.0)/len(data))
        return error

def levels(number):
        for j in range (10):
                jj=pow(2,j)

                if jj <= (number-1):
                        level = j
                else:
                        break
        return level

def dropzeros(number):
    mynum          = decimal.Decimal(number).normalize()
    # e.g 22000 --> Decimal('2.2E+4')
    return mynum.__trunc__() if not mynum % 1 else float(mynum)

def GetBconst(molecule_rot):
        '''
        This function calculates rotational Bconstant for linear rotor
        '''
        '''
        autocminverse  = 2.1947463137e+5
        energyj0       = -36117.5942855
        energyj1       = -35999.1009407
        bconst         = 0.5*(energyj1-energyj0)     # in cm^-1
        '''
        if (molecule_rot == "HF"):
                #bconst    = 20.9561                     # in cm^-1  and it is  taken from http://webbook.nist.gov/cgi/inchi?ID=C7664393&Mask=1000#Diatomic
                #bconst    = 20.9661                     # in cm^-1  and it is  taken from http://webbook.nist.gov/cgi/inchi?ID=C7664393&Mask=1000#Diatomic
                bconst     = 20.561                      # in cm^-1  and it is  taken from J. Opt. Soc. Am. Vol. 57, issue 12, page 1464, year 1967
        if (molecule_rot == "H2"):
                bconst     = 60.853
        return bconst

def replace(string_old, string_new, file1, file2):
        '''
        This function replaces old string by new string
        '''
        f1             = open(file1, 'r')
        f2             = open(file2, 'w')
        for line in f1:
                f2.write(line.replace(string_old, string_new))
        f1.close()
        f2.close()

def beads(tau,beta):
        '''
        This function determins number of beads
        '''
        numbbeads1     =beta/tau+1
        numbbeads2     = int(round(numbbeads1,0))
        if (numbbeads2 % 2 == 0):
                numbbeads2 = numbbeads2 + 1
        return numbbeads2

def file_operations(TypeCal,final_dir_in_work,numbmolecules,numbbeads):

        if (TypeCal == "ENT"):
                fileList = ["output.rden", "output.xyz"]
                file_old = final_dir_in_work+"/results/output.rden_old"
        else:
                fileList = ["output.eng", "output.xyz"]
                file_old = final_dir_in_work+"/results/output.eng_old"

        flag = False
        if (os.path.isfile(final_dir_in_work+"/results/"+fileList[0]) == True):
                flag = True
                if (os.path.isfile(file_old) == True):
                        for filecat in fileList:
                                if (filecat != "output.xyz"):
                                        col_data_new = np.genfromtxt(final_dir_in_work+"/results/"+filecat)
                                        index = int(col_data_new[0,0])
                                        col_data_old = np.genfromtxt(final_dir_in_work+"/results/"+filecat+"_old")
                                        merged_data  = np.concatenate((col_data_old[:index-1], col_data_new), axis=0)

                                        if (filecat == "output.eng"):
                                                np.savetxt(final_dir_in_work+"/results/"+filecat+"_old", merged_data, fmt='%d    %.6e    %.6e    %.6e    %.6e')
                                        else:
                                                np.savetxt(final_dir_in_work+"/results/"+filecat+"_old", merged_data, fmt='%.6e', delimiter='    ')

                                if (filecat == "output.xyz"):
                                        if "H2O1" in open(final_dir_in_work+"/results/"+filecat).read():
                                                rmstr = int(numbmolecules*numbbeads+3)
                                                file_temp = final_dir_in_work+"/results/"+filecat+"_temp"
                                                cmd1="tail -n +"+str(rmstr)+" "+final_dir_in_work+"/results/"+filecat+">"+file_temp
                                                os.system(cmd1)
                                                col_data_new = np.genfromtxt(file_temp)
                                                call(["rm", file_temp])
                                        else:
                                                col_data_new = np.genfromtxt(final_dir_in_work+"/results/"+filecat)
                                        index = int(col_data_new[0,0])
                                        col_data_old = np.genfromtxt(final_dir_in_work+"/results/"+filecat+"_old")
                                        merged_data  = np.concatenate((col_data_old[:index-1], col_data_new), axis=0)
                                        np.savetxt(final_dir_in_work+"/results/"+filecat+"_old", merged_data, fmt='%.6e', delimiter='    ')

                                call(["rm", final_dir_in_work+"/results/"+filecat])
                else:
                        for filemv in fileList:
                                call(["mv", final_dir_in_work+"/results/"+filemv, final_dir_in_work+"/results/"+filemv+"_old"])
        return flag

def fmtAverageEnergy(TypeCal,status,variable):
        '''
        This function gives us the output 
        '''
        if (variable == "Rpt"):
                unit = "(Angstrom)"
        else:
                unit = "(1/K)"

        if (status == "analysis"):
                output     ="# Unit of all kind of energies is expessed in Kelvin."
                output  += "\n"
                output    +="# "
                output  += "\n"
                output    +="# "
                if (TypeCal == "PIMC"):
                        output    += '{blocks:^10}{beads:^10}{var:^10}{rot:^16}{pot:^16}{tot:^16}{rotsq:^16}{cv:^16}{er1:^12}{er2:^12}{er3:^12}{er4:^12}{er5:^12}'.format(blocks='nBlocks',beads='nBeads', var=variable+' invK', rot='<K>', pot='<V>', tot='<E>', rotsq='<Ksq>', cv='<Cv>', er1='Err-K', er2='Err-V', er3='Err-E', er4='Err-Ksq', er5='Err-Cv',)
                        output    +="\n"
                        output    += '{0:=<170}'.format('#')
                        output    +="\n"

                if (TypeCal == "PIGS"):
                        output    += '{blocks:^10}{beads:^10}{var:^10}{pot:^16}{tot:^16}{er1:^12}{er2:^12}'.format(blocks='nBlocks',beads='nBeads', var=variable+' invK', pot='<V>', tot='<E>', er1='Err-V', er2='Err-E')
                        output    +="\n"
                        output    += '{0:=<90}'.format('#')
                        output    +="\n"

                if (TypeCal == "ENT"):
                        output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable, 'Avg. rotational', 'Avg. (E - V)', 'Avg. Potential', 'Avg. Total', 'Error of Rotational', 'Error of (E - V)', 'Error of Potential', 'Error of Total')
                        output    +="\n"
                        output    += '{0:=<90}'.format('#')
                        output    +="\n"

                return output

def GetAverageEnergy(TypeCal,numbbeads,variable,final_dir_in_work,preskip,postskip, numbblocks):
        '''
        This function gives us the output 
        '''
        if (os.path.isdir(final_dir_in_work)):
                condition = True

                file_old = final_dir_in_work+"/results/output.eng_old"
                if (os.path.isfile(file_old) == True):
                        file_old_1 = final_dir_in_work+"/results/output.eng_old_1"
                        file_new = final_dir_in_work+"/results/output.eng"
                        if (os.path.isfile(file_old_1) == True):
                                col_data_new = genfromtxt(final_dir_in_work+"/results/output.eng_old_1")
                                index = int(col_data_new[0,0])
                                col_data_old = genfromtxt(final_dir_in_work+"/results/output.eng_old")
                                marged_data  = np.concatenate((col_data_old[:index-1], col_data_new), axis=0)

                                col_data_new = genfromtxt(final_dir_in_work+"/results/output.eng")
                                index = int(col_data_new[0,0])
                                marged_data  = np.concatenate((marged_data[:index-1], col_data_new), axis=0)
                                aa = col_data_new[:,0]
                                final_data_set = marged_data[preskip:(int(aa[-1])-postskip),:]
                        elif ((os.path.isfile(file_new) == True) and (os.path.isfile(file_old_1) == False)):
                                print(final_dir_in_work + " -- Restarted data")
                                col_data_new = genfromtxt(final_dir_in_work+"/results/output.eng")
                                index = int(col_data_new[0,0])
                                col_data_old = genfromtxt(final_dir_in_work+"/results/output.eng_old")
                                marged_data  = np.concatenate((col_data_old[:index-1], col_data_new), axis=0)
                                aa = col_data_new[:,0]
                                final_data_set = marged_data[preskip:(int(aa[-1])-postskip),:]
                        elif ((os.path.isfile(file_new) == False) and (os.path.isfile(file_old_1) == False)):
                                final_data_set = genfromtxt(final_dir_in_work+"/results/output.eng_old", skip_header=preskip, skip_footer=postskip)
                else:
                        final_data_set = genfromtxt(final_dir_in_work+"/results/output.eng", skip_header=preskip, skip_footer=postskip)

        if (TypeCal == "PIMC"):
                col_block = final_data_set[:,0] 
                col_rot   = final_data_set[:,1] 
                col_pot   = final_data_set[:,2] 
                col_tot   = final_data_set[:,3] 
                col_rotsq = final_data_set[:,4] 
                col_Cv1   = final_data_set[:,5] 
                col_Cv2   = final_data_set[:,6] 

                ncol_block = len(col_block)
                if (int(len(col_block)) != numbblocks-(preskip+postskip)):
                        print(len(col_block))
                        print(final_dir_in_work)
        
                workingNdim   = int(math.log(len(col_tot))/math.log(2))
                trunc         = int(len(col_tot)-2**workingNdim)
        
                col_rot       = col_rot[trunc:]
                col_pot       = col_pot[trunc:]
                col_tot       = col_tot[trunc:]
                col_rotsq     = col_rotsq[trunc:]
                col_Cv1       = col_Cv1[trunc:]
                col_Cv2       = col_Cv2[trunc:]

                mean_rot      = np.mean(col_rot)
                mean_pot      = np.mean(col_pot)
                mean_tot      = np.mean(col_tot)
                mean_rotsq    = np.mean(col_rotsq)
                mean_Cv1      = np.mean(col_Cv1)
                mean_Cv2      = np.mean(col_Cv2)
                mcbeta        = variable*numbbeads
                mean_Cv       = mcbeta*mcbeta*(mean_Cv1-mean_Cv2-mean_tot*mean_tot)

                error_rot     = maxError_byBining(mean_rot, col_rot, workingNdim-6)
                error_pot     = maxError_byBining(mean_pot, col_pot, workingNdim-6)
                error_tot     = maxError_byBining(mean_tot, col_tot, workingNdim-6)
                error_rotsq   = maxError_byBining(mean_rotsq, col_rotsq, workingNdim-6)
                error_Cv1     = maxError_byBining(mean_Cv1, col_Cv1, workingNdim-6)
                error_Cv2     = maxError_byBining(mean_Cv2, col_Cv2, workingNdim-6)
                error_Cv      = mcbeta*mcbeta*(math.sqrt(error_Cv1*error_Cv1+error_Cv2*error_Cv2+4.0*mean_tot*mean_tot*error_tot*error_tot))

                output  = '{blocks:^12d}{beads:^10d}{var:^10.6f}{rot:^16.6f}{pot:^16.6f}{tot:^16.6f}{rotsq:^16.6f}{Cv:^16.6f}{er1:^12.6f}{er2:^12.6f}{er3:^12.6f}{er4:^12.6f}{er5:^12.6f}'.format(blocks=ncol_block,beads=numbbeads, var=variable, rot=mean_rot, pot=mean_pot, tot=mean_tot, rotsq=mean_rotsq, Cv=mean_Cv, er1=error_rot, er2=error_pot, er3=error_tot, er4=error_rotsq, er5=error_Cv)
                output  += "\n"

        if (TypeCal == "PIGS"):
                col_block = final_data_set[:,0] 
                col_pot   = final_data_set[:,1]
                col_tot   = final_data_set[:,2]

                ncol_block = len(col_block)
                if (int(len(col_block)) != numbblocks-(preskip+postskip)):
                        print(len(col_block))
                        print(final_dir_in_work)
        
                workingNdim   = int(math.log(len(col_tot))/math.log(2))
                trunc         = int(len(col_tot)-2**workingNdim)
        
                col_pot       = col_pot[trunc:]
                col_tot       = col_tot[trunc:]

                mean_pot      = np.mean(col_pot)
                mean_tot      = np.mean(col_tot)

                error_pot     = maxError_byBining(mean_pot, col_pot, workingNdim-6)
                error_tot     = maxError_byBining(mean_tot, col_tot, workingNdim-6)

                output  = '{blocks:^12d}{beads:^10d}{var:^10.6f}{pot:^16.6f}{tot:^16.6f}{er1:^12.6f}{er2:^12.6f}'.format(blocks=ncol_block,beads=numbbeads, var=variable, pot=mean_pot, tot=mean_tot, er1=error_pot, er2=error_tot)
                output  += "\n"

        return output
        
def fmtAverageOrderParam(status,variable):
        '''
        This function gives us the output 
        '''
        if variable == "Rpt":
                unit = "(Angstrom)"
        else:
                unit = "(1/K)"

        if status == "analysis":
                output    ="# "
                output    += '{blocks:^10}{beads:^10}{var:^10}{eiz:^12}{eiejz:^12}{er1:^12}{er2:^12}'.format(blocks='nBlocks',beads='nBeads', var=variable+' invK', eiz='<eiz>', eiejz='<eiejz>', er1='Err-eiz', er2='Err-eiejz')
                output    +="\n"
                output    += '{0:=<80}'.format('#')
                output    +="\n"
                return output

def GetAverageOrderParam(TypeCal,numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip,numbblocks):
        '''
        See PRL 118, 027402 (2017) 
        '''
        if (TypeCal != 'PIMC'):
                axis_index = {"cost":0, "phi":1, "chi":2}
                axis_read = "cost"
                ndofs = 3
                beads_pos = int((numbbeads-1)/2)
                collist=[0]
                for i in range(numbmolecules):
                        ncol1 = beads_pos+i*numbbeads
                        ncol = axis_index[axis_read]+ncol1*ndofs
                        ncol = ncol+1
                        collist.append(ncol)
                        #print(str(ncol)+'th column')
                col=tuple(collist)

                if (os.path.isdir(final_dir_in_work)):
                        file_old = final_dir_in_work+"/results/output.xyz_old"
                        if (os.path.isfile(file_old) == True):
                                file_new = final_dir_in_work+"/results/output.xyz"
                                if (os.path.isfile(file_new) == True):
                                        print(final_dir_in_work + " -- Restarted data")
                                        if "H2O1" in open(file_new).read():
                                                file_temp = final_dir_in_work+"/results/output_temp.xyz"
                                                rmstr = int(numbmolecules*numbbeads+3)
                                                cmd1="tail -n +"+str(rmstr)+" "+file_new+">"+file_temp
                                                os.system(cmd1)
                                                col_data_new = np.genfromtxt(file_temp, usecols = col)
                                                call(["rm", file_temp])
                                        else:
                                                col_data_new = np.genfromtxt(file_new, usecols = col)
                                        index = int(col_data_new[0,0])
                                        col_data_old = np.genfromtxt(file_old,usecols = col)
                                        marged_data  = np.concatenate((col_data_old[:index-1], col_data_new), axis=0)
                                        aa = col_data_new[:,0]
                                        final_data_set = marged_data[preskip:(int(aa[-1])-postskip),:]
                                else:
                                        final_data_set = genfromtxt(file_old, usecols = col, skip_header=preskip, skip_footer=postskip)
                        else:
                                final_data_set = genfromtxt(final_dir_in_work+"/results/output.xyz", usecols = col, skip_header=preskip, skip_footer=postskip)

        ncol_block = len(final_data_set[:,0])
        if (int(ncol_block) != numbblocks-(preskip+postskip)):
                print(ncol_block)
                print(final_dir_in_work)


        workingNdim   = int(math.log(ncol_block)/math.log(2))
        trunc         = int(ncol_block-2**workingNdim)
        raw_data=np.delete(final_data_set, 0, 1)

        if (numbmolecules == 2):
                raw_data1 = np.absolute(raw_data)
                eiz=np.sum(raw_data1[trunc:,:],axis=1)/numbmolecules
        if (numbmolecules == 11):
                raw_data1 = np.absolute(raw_data)
                eiz=np.sum(raw_data1[trunc:,2:numbmolecules-3],axis=1)/len([i for i in range(2,numbmolecules-2)])
        mean_eiz = np.mean(eiz)
        error_eiz = maxError_byBining(mean_eiz, eiz, workingNdim-6)

        if (numbmolecules == 2):
                paireiej = [i for i in range(numbmolecules-1)]
        if (numbmolecules == 11):
                paireiej = [i for i in range(2,numbmolecules-3)]
        norm_eiejz = len(paireiej)
        eiejz=np.zeros(ncol_block-trunc,dtype=float)
        for i in paireiej:
                eiejz += np.multiply(raw_data[trunc:,i],raw_data[trunc:,i+1])/norm_eiejz
        mean_eiejz = np.mean(eiejz)
        error_eiejz = maxError_byBining(mean_eiejz, eiejz, workingNdim-6)

        output  = '{blocks:^12d}{beads:^10d}{var:^10.6f}{eiz:^12.6f}{eiejz:^12.6f}{er1:^12.6f}{er2:^12.6f}'.format(blocks=ncol_block,beads=numbbeads, var=variable, eiz=mean_eiz,eiejz=mean_eiejz,er1=error_eiz,er2=error_eiejz)
        output  += "\n"
        return output

'''
def GetAverageOrderParam(TypeCal,numbmolecules,numbbeads,variable,final_dir_in_work,preskip,postskip,numbblocks):
        #See PRL 118, 027402 (2017) 
        if (os.path.isdir(final_dir_in_work)):
                condition = True

                file_old = final_dir_in_work+"/results/outputOrderPara.corr_old"
                if (os.path.isfile(file_old) == True):
                        file_old_1 = final_dir_in_work+"/results/outputOrderPara.corr_old_1"
                        file_new = final_dir_in_work+"/results/outputOrderPara.corr"
                        if (os.path.isfile(file_old_1) == True):
                                col_data_new = genfromtxt(final_dir_in_work+"/results/outputOrderPara.corr_old_1")
                                index = int(col_data_new[0,0])
                                col_data_old = genfromtxt(final_dir_in_work+"/results/outputOrderPara.corr_old")
                                marged_data  = np.concatenate((col_data_old[:index-1], col_data_new), axis=0)

                                col_data_new = genfromtxt(final_dir_in_work+"/results/outputOrderPara.corr")
                                index = int(col_data_new[0,0])
                                marged_data  = np.concatenate((marged_data[:index-1], col_data_new), axis=0)
                                aa = col_data_new[:,0]
                                final_data_set = marged_data[preskip:(int(aa[-1])-postskip),:]
                        elif ((os.path.isfile(file_new) == True) and (os.path.isfile(file_old_1) == False)):
                                print(final_dir_in_work + " -- Restarted data")
                                col_data_new = genfromtxt(final_dir_in_work+"/results/outputOrderPara.corr")
                                index = int(col_data_new[0,0])
                                col_data_old = genfromtxt(final_dir_in_work+"/results/outputOrderPara.corr_old")
                                marged_data  = np.concatenate((col_data_old[:index-1], col_data_new), axis=0)
                                aa = col_data_new[:,0]
                                final_data_set = marged_data[preskip:(int(aa[-1])-postskip),:]
                        elif ((os.path.isfile(file_new) == False) and (os.path.isfile(file_old_1) == False)):
                                final_data_set = genfromtxt(final_dir_in_work+"/results/outputOrderPara.corr_old", skip_header=preskip, skip_footer=postskip)
                else:
                        final_data_set = genfromtxt(final_dir_in_work+"/results/outputOrderPara.corr", skip_header=preskip, skip_footer=postskip)
        if (TypeCal != 'PIMC'):
                if (os.path.isdir(final_dir_in_work)):
                        condition = True

                        file_old = final_dir_in_work+"/results/output.xyz_old"
                        if (os.path.isfile(file_old) == True):
                                file_new = final_dir_in_work+"/results/output.xyz"
                                if (os.path.isfile(file_new) == True):
                                        print(final_dir_in_work + " -- Restarted data")
                                        col_data_new = genfromtxt(final_dir_in_work+"/results/output.xyz")
                                        index = int(col_data_new[0,0])
                                        col_data_old = genfromtxt(final_dir_in_work+"/results/output.xyz_old")
                                        marged_data  = np.concatenate((col_data_old[:index-1], col_data_new), axis=0)
                                        aa = col_data_new[:,0]
                                        final_data_set = marged_data[preskip:(int(aa[-1])-postskip),:]
                                elif ((os.path.isfile(file_new) == False) and (os.path.isfile(file_old_1) == False)):
                                        final_data_set = genfromtxt(final_dir_in_work+"/results/output.xyz_old", skip_header=preskip, skip_footer=postskip)
                        else:
                                final_data_set = genfromtxt(final_dir_in_work+"/results/output.xyz", skip_header=preskip, skip_footer=postskip)


        axis_index = {"cost":0, "phi":1, "chi":2}
        axis_read = "cost"
        ndofs = 3
        beads_pos = int((numbbeads-1)/2)
        for i in range(numbmolecules):
                ncol1 = beads_pos+i*numbbeads
                ncol = axis_index[axis_read]+ncol1*ndofs
                ncol = ncol+1
                #print(str(ncol)+'th column')
                save_data[i,:] = final_data_set[(preskip+trunc):(nlen-postskip),ncol]

        col_block = final_data_set[:,0] 
        col_eiejx = final_data_set[:,1] 
        col_eiejy = final_data_set[:,2] 
        col_eiejz = final_data_set[:,3] 
        col_eiej  = final_data_set[:,4] 
        col_eix   = final_data_set[:,5] 
        col_eiy   = final_data_set[:,6] 
        col_eiz   = final_data_set[:,7] 

        ncol_block = len(col_block)
        if (int(len(col_block)) != numbblocks-(preskip+postskip)):
                print(len(col_block))
                print(final_dir_in_work)

        workingNdim   = int(math.log(len(col_eiz))/math.log(2))
        trunc         = int(len(col_eiz)-2**workingNdim)
        
        col_eiejx = col_eiejx[trunc:]
        col_eiejy = col_eiejy[trunc:]
        col_eiejz = col_eiejz[trunc:]
        col_eiej = col_eiej[trunc:]

        mean_eiejx = np.mean(col_eiejx)
        mean_eiejy = np.mean(col_eiejy)
        mean_eiejz = np.mean(col_eiejz)
        mean_eiej = np.mean(col_eiej)

        error_eiejx = maxError_byBining(mean_eiejx, col_eiejx, workingNdim-6)
        error_eiejy = maxError_byBining(mean_eiejy, col_eiejy, workingNdim-6)
        error_eiejz = maxError_byBining(mean_eiejz, col_eiejz, workingNdim-6)
        error_eiej = maxError_byBining(mean_eiej, col_eiej, workingNdim-6)

        col_eix = col_eix[trunc:]
        col_eiy = col_eiy[trunc:]
        col_eiz = col_eiz[trunc:]

        mean_eix = np.mean(col_eix)
        mean_eiy = np.mean(col_eiy)
        mean_eiz = np.mean(col_eiz)

        error_eix = maxError_byBining(mean_eix, col_eix, workingNdim-6)
        error_eiy = maxError_byBining(mean_eiy, col_eiy, workingNdim-6)
        error_eiz = maxError_byBining(mean_eiz, col_eiz, workingNdim-6)

        output  = '{blocks:^12d}{beads:^10d}{var:^10.6f}{eiejx:^12.6f}{eiejy:^12.6f}{eiejz:^12.6f}{eiej:^12.6f}{eix:^12.6f}{eiy:^12.6f}{eiz:^12.6f}{er1:^12.6f}{er2:^12.6f}{er3:^12.6f}{er4:^12.6f}{er5:^12.6f}{er6:^12.6f}{er7:^12.6f}'.format(blocks=ncol_block,beads=numbbeads, var=variable, eiejx=mean_eiejx, eiejy=mean_eiejy, eiejz=mean_eiejz, eiej=mean_eiej, eix=mean_eix, eiy=mean_eiy, eiz=mean_eiz, er1=error_eiejx, er2=error_eiejy, er3=error_eiejz, er4=error_eiej, er5=error_eix, er6=error_eiy, er7=error_eiz)
        output  += "\n"
        return output
'''

def GetAverageEntropy(numbbeads,variable,final_dir_in_work,preskip,postskip,numbblocks,ENT_TYPE):
        '''
        This function gives us the output 
        '''
        if (ENT_TYPE == "EXTENDED_ENSMBL"):
                col_block, col_nm, col_dm = genfromtxt(final_dir_in_work+"/results/output.rden",unpack=True, usecols=[0,1,2], skip_header=preskip, skip_footer=postskip)
                ncol_block = len(col_block)
                if (int(len(col_block)) != numbblocks-(preskip+postskip)):
                        print(len(col_block))
                        print(final_dir_in_work)

                workingNdim  = int(math.log(len(col_nm))/math.log(2))
                trunc        = int(len(col_nm)-2**workingNdim)
        
                col_nm       = col_nm[trunc:]
                col_dm       = col_dm[trunc:]
                mean_nm      = np.mean(col_nm)
                mean_dm      = np.mean(col_dm)
                purity       = mean_nm/mean_dm
                mean_EN      = -log(purity)

                error_nm     = maxError_byBining(mean_nm, col_nm, workingNdim-6) 
                error_dm     = maxError_byBining(mean_dm, col_dm, workingNdim-6)
                error_purity = abs(purity)*sqrt((error_dm/mean_dm)*(error_dm/mean_dm) + (error_nm/mean_nm)*(error_nm/mean_nm))
                error_EN     = sqrt((error_dm/mean_dm)*(error_dm/mean_dm) + (error_nm/mean_nm)*(error_nm/mean_nm))

                output  = '{blocks:^8d}{beads:^10d}{var:^10.6f}{nm:^12.6f}{dm:^12.6f}{pur:^12.6f}{ent:^12.6f}{er1:^12.6f}{er2:^12.6f}{er3:^12.6f}{er4:^12.6f}'.format(blocks=ncol_block,beads=numbbeads, var=variable, nm=mean_nm, dm=mean_dm, pur=purity, ent=mean_EN, er1=error_nm, er2=error_dm, er3=error_purity, er4=error_EN)
                output  += "\n"

        if (ENT_TYPE == 'BROKENPATH'):
                col_block, col_nm, col_dm = genfromtxt(final_dir_in_work+"/results/output.rden",unpack=True, usecols=[0,1,2], skip_header=preskip, skip_footer=postskip)
                workingNdim  = int(math.log(len(col_nm))/math.log(2))
                trunc        = int(len(col_nm)-2**workingNdim)
        
                col_nm       = col_nm[trunc:]
                col_dm       = col_dm[trunc:]
                mean_nm      = np.mean(col_nm)
                mean_dm      = np.mean(col_dm)
                mean_EN      = -log(mean_nm/mean_dm)

                error_nm     = maxError_byBining(mean_nm, col_nm, workingNdim-6) 
                error_dm     = maxError_byBining(mean_dm, col_dm, workingNdim-6)
                error_EN     = sqrt((error_dm/mean_dm)*(error_dm/mean_dm) + (error_nm/mean_nm)*(error_nm/mean_nm))

                output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}'.format(numbbeads, variable, mean_nm, mean_dm, mean_EN, error_nm, error_dm, error_EN)
                output  += "\n"

        if (ENT_TYPE == "SWAP"):
                col_block, col_nm, col_dm, col_TrInv = genfromtxt(final_dir_in_work+"/results/output.rden",unpack=True, usecols=[0,1,2,3], skip_header=preskip, skip_footer=postskip)
                workingNdim  = int(math.log(len(col_nm))/math.log(2))
                trunc        = int(len(col_nm)-2**workingNdim)
        
                col_nm       = col_nm[trunc:]
                col_dm       = col_dm[trunc:]
                mean_nm      = np.mean(col_nm)
                mean_dm      = np.mean(col_dm)
                mean_TrInv   = np.mean(col_TrInv)
                purity       = 1.0/mean_TrInv
                mean_EN      = -log(purity)

                error_nm     = maxError_byBining(mean_nm, col_nm, workingNdim-6) 
                error_dm     = maxError_byBining(mean_dm, col_dm, workingNdim-6)
                error_TrInv  = np.std(col_TrInv,ddof=1)/sqrt(len(col_block))
                error_purity = abs(1.0/(mean_TrInv*mean_TrInv))*error_TrInv
                error_EN     = abs(1.0/mean_TrInv)*error_TrInv #Write the proper equation

                output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(numbbeads, variable, mean_nm, mean_dm, purity, mean_EN, error_nm, error_dm, error_purity, error_EN)
                output  += "\n"

        if (ENT_TYPE == "UNSWAP"):
                col_block, col_nm, col_dm, col_Tr = genfromtxt(final_dir_in_work+"/results/output.rden",unpack=True, usecols=[0,1,2,3], skip_header=preskip, skip_footer=postskip)
                workingNdim  = int(math.log(len(col_nm))/math.log(2))
                trunc        = int(len(col_nm)-2**workingNdim)
        
                col_nm       = col_nm[trunc:]
                col_dm       = col_dm[trunc:]
                mean_nm      = np.mean(col_nm)
                mean_dm      = np.mean(col_dm)
                mean_purity  = np.mean(col_Tr)
                mean_EN      = -log(mean_Tr)

                error_nm     = maxError_byBining(mean_nm, col_nm, workingNdim-6) 
                error_dm     = maxError_byBining(mean_dm, col_dm, workingNdim-6)
                error_purity = np.std(col_Tr,ddof=1)/sqrt(len(col_block))
                error_EN     = abs(1.0/mean_TrInv)*error_TrInv #Write the proper equation

                output  = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(numbbeads, variable, mean_nm, mean_dm, mean_purity, mean_EN, error_nm, error_dm, error_purity, error_EN)
                output  += "\n"

        return output

def fmtAverageEntropy(status,variable,ENT_TYPE):
        '''
        This function gives us the output 
        '''
        if variable == "Rpt":
                unit = "(Angstrom)"
        else:
                unit = "(1/K)"

        if status == "analysis":
                output     ="#"
                if (ENT_TYPE == "EXTENDED_ENSMBL"):
                        output    += '{blocks:^8}{beads:^10}{var:^10}{nm:^12}{dm:^12}{pur:^12}{ent:^12}{er1:^12}{er2:^12}{er3:^12}{er4:^12}'.format(blocks='nBlocks',beads='nBeads', var=variable+' invK', nm='<Nm>', dm='<Dm>', pur='<Purity>', ent='<S2>', er1='Err-Nm', er2='Err-Dm', er3='Err-Purity', er4='Err-S2')
                if (ENT_TYPE == "BROKENPATH"):
                        output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Entropy', 'Error of Nm', 'Error of Dm', 'Error of Entropy')
                if (ENT_TYPE == "SWAP"):
                        output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Avg. Purity', 'Entropy', 'Error of Nm', 'Error of Dm', 'Error of Purity', 'Error of Entropy')
                if (ENT_TYPE == "REGULARPATH"):
                        output    += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format('Beads', variable+'  (1/K)', '<Nm>', '<Dm>', 'Avg. Purity', 'Entropy', 'Error of Nm', 'Error of Dm', 'Error of Purity', 'Error of Entropy')
                output    +="\n"
                output    += '{0:=<123}'.format('#')
                output    +="\n"
                return output

def GetAverageEntropyRT(particleAList, TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, preskip, postskip, extra_file_name, dir_output, variable, crystal, final_results_path,impurity, ext_ent):
        '''
        This function gives us Renyi entropy of a many-rotors system simulated by Ratio trick algorithm.
        '''
        list_nb           = inputFile.Getbeads(TypeCal, variableName)
        ndim_beads        = int(len(list_nb))

        purity_combo      = np.zeros(ndim_beads,dtype = 'f')
        err_purity_combo  = np.zeros(ndim_beads,dtype = 'f')
        entropy_combo     = np.zeros(ndim_beads,dtype = 'f')
        err_entropy_combo = np.zeros(ndim_beads,dtype = 'f')
        col_beads         = np.zeros(ndim_beads,dtype = int)
        col_var           = np.zeros(ndim_beads,dtype = 'f')
        ii = 0

        for iBead in list_nb:
                if ((iBead % 2) != 0):
                        value     = iBead
                else:
                        value     = iBead+1

                if (variableName == "beta"):
                        beta      = parameter*(value - 1)
                        variable  = beta
                if (variableName == "tau"):
                        tau       = parameter/(value-1)
                        variable  = tau

                numbbeads     = value

                col_purity= np.zeros(len(particleAList),dtype = 'f')
                col_err_purity= np.zeros(len(particleAList),dtype = 'f')

                contition = False
                iPartition = 0
                for partition in particleAList:
                        file1_name = GetFileNameSubmission(TypeCal, molecule_rot, TransMove, RotMove, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, partition, extra_file_name, crystal, impurity, ext_ent)
                        folder_run = file1_name+str(numbbeads)
                        final_dir_in_work = dir_output+folder_run
                        if os.path.isdir(final_dir_in_work):
                                condition = True

                                file_old = final_dir_in_work+"/results/output.rden_old"
                                if os.path.isfile(file_old) == True:
                                        file_old_1 = final_dir_in_work+"/results/output.rden_old_1"
                                        file_new = final_dir_in_work+"/results/output.rden"
                                        if os.path.isfile(file_old_1) == True:
                                                col_data_new = genfromtxt(final_dir_in_work+"/results/output.rden_old_1")
                                                index = int(col_data_new[0,0])
                                                col_data_old = genfromtxt(final_dir_in_work+"/results/output.rden_old")
                                                marged_data  = np.concatenate((col_data_old[:index-1], col_data_new), axis=0)

                                                col_data_new = genfromtxt(final_dir_in_work+"/results/output.rden")
                                                index = int(col_data_new[0,0])
                                                marged_data  = np.concatenate((marged_data[:index-1], col_data_new), axis=0)
                                                aa = col_data_new[:,0]
                                                col_block    = marged_data[preskip:(int(aa[-1])-postskip),0]
                                                col_nm       = marged_data[preskip:(int(aa[-1])-postskip),1]
                                                col_dm       = marged_data[preskip:(int(aa[-1])-postskip),2]
                                        elif ((os.path.isfile(file_new) == True) and (os.path.isfile(file_old_1) == False)):
                                                col_data_new = genfromtxt(final_dir_in_work+"/results/output.rden")
                                                index = int(col_data_new[0,0])
                                                col_data_old = genfromtxt(final_dir_in_work+"/results/output.rden_old")
                                                marged_data  = np.concatenate((col_data_old[:index-1], col_data_new), axis=0)
                                                aa = col_data_new[:,0]
                                                col_block    = marged_data[preskip:(int(aa[-1])-postskip),0]
                                                col_nm       = marged_data[preskip:(int(aa[-1])-postskip),1]
                                                col_dm       = marged_data[preskip:(int(aa[-1])-postskip),2]
                                        elif ((os.path.isfile(file_new) == False) and (os.path.isfile(file_old_1) == False)):
                                                col_block, col_nm, col_dm = genfromtxt(final_dir_in_work+"/results/output.rden_old",unpack=True, usecols=[0,1,2], skip_header=preskip, skip_footer=postskip)
                                else:
                                        col_block, col_nm, col_dm = genfromtxt(final_dir_in_work+"/results/output.rden",unpack=True, usecols=[0,1,2], skip_header=preskip, skip_footer=postskip)
                                if (int(len(col_block)) != numbblocks-(preskip+postskip)):
                                        print(len(col_block))
                                        print(folder_run)

                                workingNdim  = int(math.log(len(col_nm))/math.log(2))
                                trunc        = int(len(col_nm)-2**workingNdim)
                                mean_nm = np.mean(col_nm[trunc:])
                                mean_dm = np.mean(col_dm[trunc:])
                                err_nm  = maxError_byBining(mean_nm, col_nm[trunc:], workingNdim-6) 
                                err_dm  = maxError_byBining(mean_dm, col_dm[trunc:], workingNdim-6)

                                col_purity[iPartition]      = mean_nm/mean_dm
                                col_err_purity[iPartition]  = abs(mean_nm/mean_dm)*sqrt((err_dm/mean_dm)*(err_dm/mean_dm) + (err_nm/mean_nm)*(err_nm/mean_nm))
                        else:
                                condition = False
                        iPartition += 1
                        
                if (condition == True):
                        if (len(particleAList) > 1):
                                purity       = np.prod(col_purity,axis=0)
                                entropy      = -log(purity)
                                err_purity   = abs(purity)*sqrt(np.sum(np.square(np.divide(col_err_purity,col_purity))))
                                err_entropy  = abs(err_purity)/purity
                        else:
                                purity       = col_purity[0]
                                entropy      = -log(purity)
                                err_purity   = col_err_purity[0]
                                err_entropy  = abs(err_purity)/purity

                        purity_combo[ii]      = purity
                        entropy_combo[ii]     = entropy
                        col_beads[ii]         = numbbeads
                        col_var[ii]           = variable
                        err_purity_combo[ii]  = err_purity
                        err_entropy_combo[ii] = err_entropy
                        ii = ii+1

        #extra_file_name = 'Ratio-Trick-'       
        FileAnalysis = GetFileNameAnalysis(TypeCal, True, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, partition, ext_ent)
        headerString = fmtAverageEntropyRT('tau')
        np.savetxt(FileAnalysis.SaveEntropyRT, np.transpose([col_beads[:ii], col_var[:ii], purity_combo[:ii], entropy_combo[:ii], err_purity_combo[:ii], err_entropy_combo[:ii]]), fmt=['%4d','%10.6f', '%10.6f', '%10.6f', '%10.6f', '%10.6f'],header=headerString)
        SavedFile = FileAnalysis.SaveEntropyRT
        FileCheck(TypeCal,list_nb,variableName,SavedFile)
        call(["cat", FileAnalysis.SaveEntropyRT])
        print("Successful execution!")

def fmtAverageEntropyRT(variable):
        '''
        This function gives us the output 
        '''
        output    = '{0:5}{1:10}{2:10}{3:10}{4:10}{5:30}'.format('P', variable+'(1/K)', 'Purity', 'Entropy', 'ErrorPurity', '  ErrorEntropy')
        return output

def GetInput(TypeCal, ENT_TYPE, ENT_ALGR, temperature,numbbeads,numbblocks,numbpass,molecule_rot,numbmolecules,distance,level,step,step_trans,gfact,dipolemoment,particleA, Restart1, numbblocks_Restart1, crystal, RotorType, TransMove, RotMove, path_Density,impurity,step_trans1, level1):
        '''
        This function modifies parameters in qmc_run.input
        '''
        src_dir = os.getcwd()+"/"

        replace("temperature_input", str(temperature), "qmc_run.input", "qmc2.input")
        replace("numbbeads_input", str(numbbeads), "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])
        if (molecule_rot == "CH3F"):
                replace("potread_input", "pesch3fph2-180", "qmc2.input", "qmc3.input")
        elif (molecule_rot == "H2O"):
                replace("potread", "nopot", "qmc2.input", "qmc3.input")
        elif (molecule_rot == "HF"):
                replace("potread_input", "nopot", "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        replace("PathToPot", src_dir, "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        replace("sim_input", TypeCal+"_SIM", "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        replace("sim_ensmbl_input", ENT_TYPE, "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        if (ENT_ALGR == "WR"):
                replace("sim_algr_input", "RATIOTRICK", "qmc2.input", "qmc3.input")
        elif (ENT_ALGR == "WOR"):
                replace("sim_algr_input", "NORATIOTRICK", "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        if TransMove:
                replace("#TRANSLATION", "TRANSLATION", "qmc2.input", "qmc3.input")
        else:
                replace("#TRANSLATION", "#TRANSLATION", "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        if RotMove:
                replace("cal_type_input", "ROTATION", "qmc2.input", "qmc3.input")
        else:
                replace("cal_type_input", "#ROTATION", "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        if RotMove:
                replace("den_path_input", path_Density, "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        if Restart1:
                replace("numbblocks_input", str(numbblocks_Restart1), "qmc2.input", "qmc3.input")
        else:
                replace("numbblocks_input", str(numbblocks),          "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        replace("numbmolecules_input", str(numbmolecules),        "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        if (distance >= 0.0):
                replace("distanceArg_input", "DISTANCE",              "qmc2.input", "qmc3.input")
                call(["mv", "qmc3.input", "qmc2.input"])
                replace("distance_input", str(distance),              "qmc2.input", "qmc3.input")
        else:
                replace("distanceArg_input", "",                      "qmc2.input", "qmc3.input")
                call(["mv", "qmc3.input", "qmc2.input"])
                replace("distance_input", "",                         "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        replace("type_input", RotorType,                          "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        replace("molecule_input", str(molecule_rot),              "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        replace("level_input", str(level),                        "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        replace("dstep_input", str(step),                         "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        replace("dstep_tr_input", str(step_trans),                "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        if ((dipolemoment < 0.0) and (gfact < 0.0)):
                replace("dipolemomentArg_input", "",                  "qmc2.input", "qmc3.input")
                call(["mv", "qmc3.input", "qmc2.input"])
                replace("dipolemoment_input", "",                     "qmc2.input", "qmc3.input")
        else:
                if (dipolemoment >= 0.0):
                        replace("dipolemomentArg_input", "DIPOLEMOMENT",  "qmc2.input", "qmc3.input")
                        call(["mv", "qmc3.input", "qmc2.input"])
                        replace("dipolemoment_input", str(dipolemoment),  "qmc2.input", "qmc3.input")

                if (gfact >= 0.0):
                        replace("dipolemomentArg_input", "DIPOLEMOMENT",  "qmc2.input", "qmc3.input")
                        call(["mv", "qmc3.input", "qmc2.input"])
                        dipolemoment = GetDipoleMomentFromGFactor(molecule_rot, distance, gfact)
                        replace("dipolemoment_input", str(dipolemoment),  "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        replace("numbpass_input", str(numbpass),                  "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        mcskip = numbbeads*numbpass
        replace("mskip_input", str(mcskip),                       "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        mcskip_avg = numbbeads*numbpass
        replace("mskip_avg_input", str(mcskip_avg),               "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        replace("numbparticle_input", str(particleA),             "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        if (crystal == True):
                replace("read_input", "READMCCOORDS",                 "qmc2.input", "qmc3.input")
        else:
                replace("read_input", "",                             "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])


        if Restart1:
                replace("job_input", "RESTART",                       "qmc2.input", "qmc3.input")
        else:
                replace("job_input", "START",                         "qmc2.input", "qmc3.input")
        call(["mv", "qmc3.input", "qmc2.input"])

        if (len(impurity) == 4):
                replace("#type_impurity", impurity[0], "qmc2.input", "qmc3.input")
                call(["mv", "qmc3.input", "qmc2.input"])

                replace("molecule_impurity", impurity[1], "qmc2.input", "qmc3.input")
                call(["mv", "qmc3.input", "qmc2.input"])

                replace("numbmolecules_impurity", str(impurity[2]), "qmc2.input", "qmc3.input")
                call(["mv", "qmc3.input", "qmc2.input"])

                replace("stat_impurity", impurity[3], "qmc2.input", "qmc3.input")
                call(["mv", "qmc3.input", "qmc2.input"])

                replace("dstep_tr_impurity", str(step_trans1), "qmc2.input", "qmc3.input")
                call(["mv", "qmc3.input", "qmc2.input"])

                replace("level_impurity", str(level1), "qmc2.input", "qmc3.input")
                call(["mv", "qmc3.input", "qmc2.input"])

                replace("potread_impurity", "nopot", "qmc2.input", "qmc3.input")
                call(["mv", "qmc3.input", "qmc2.input"])
                
        call(["mv", "qmc2.input", "qmc.input"])

def rotmat(TypeCal,molecule,temperature,numbbeads,source_dir_exe):
        '''
        This function generates rotational density matrix - linden.dat
        '''
        #temperature1    = dropzeros(temperature)
        temperature1    = "%8.6f" % temperature
        if (TypeCal == 'PIMC'):
                numbbeads1              = numbbeads
        else:
                numbbeads1              = numbbeads - 1
        command_linden_run = source_dir_exe+"linear_prop/linden.x "+str(temperature)+" "+str(numbbeads1)+" "+str(GetBconst(molecule))+" 150000 -1"
        system(command_linden_run)
        file_rotdens    = molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
        call(["mv", "linden.out", file_rotdens])

def GetTwoBodyDensity(Rpt, DipoleMoment, numbbeads, lmax, ltotalmax, tau, molecule):
        '''
        This function generates pair density of two body Hamiltonian
        command_line = "python script_PairDensityGenerator.py -d 1.0 -R 10.05 -P 20 -l-max 2 tau HF 0.2"
        system(command_line)
        '''
        srcCodePath         = os.path.expanduser("~")+"/DipoleChain.jl-master/examples/"
        Units               = GetUnitConverter()
        BConstant           = GetBconst(molecule)  # in wavenumber
        BConstantK          = BConstant*Units.CMRECIP2KL
        ################################################################################
        RFactorList         = GetrAndgFactor(molecule, Rpt, DipoleMoment)
        RFactor             = RFactorList[0]
        tau                     = tau*BConstantK
        if ltotalmax == 0:
                commandRun              = "julia "+srcCodePath+"pair_density.jl -R "+str(RFactor)+" --l-max "+str(lmax)+" --tau "+str(tau)
        else:
                commandRun              = "julia "+srcCodePath+"pair_density.jl -R "+str(RFactor)+" --l-max "+str(lmax)+" --l-total-max "+str(ltotalmax)+" --tau "+str(tau)
        system(commandRun)

def cagepot(source_dir_exe):
        '''
        This function generates tabulated potential - cagepot.dat
        '''
        command_cagepot_run = source_dir_exe+"tabulated_potential/hfc60.x 100 360"
        system(command_cagepot_run)
        file_cagepot    = "hfc60.pot"
        call(["mv", "cagepot.out", file_cagepot])

def jobstring_scratch(file_name, value, thread, run_dir, molecule, temperature, numbbeads, final_dir, dir_run_input_pimc,status_cagepot):
        '''
        This function creats jobstring for #PBS script
        '''
        if (thread > 4):
                thread = 4
        job_name       = "job_"+str(file_name)+str(value)
        walltime       = "200:00:00"
        processors     = "nodes=1:ppn="+str(thread)
        omp_thread     = str(thread)
        output_dir     = run_dir+"/results"
        temperature1   = "%5.3f" % temperature
        logpath        = final_dir+"/"

        input_file     = dir_run_input_pimc+"/qmc"+file_name+str(value)+".input"
        exe_file       = dir_run_input_pimc+"/pimc"
        qmcinp         = "qmc"+file_name+str(value)+".input"

        job_string     = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
##PBS -q medium
#PBS -l %s
#PBS -o %s%s.out
#PBS -e %s%s.err
export OMP_NUM_THREADS=%s
rm -rf %s
mkdir -p %s
cd $PBS_O_WORKDIR
mv %s %s
mv %s %s 
cd %s
cp %s qmc.input
cp %s %s
./pimc
mv %s /work/tapas/linear_rotors
""" % (job_name, walltime, processors, logpath, job_name, logpath, job_name, omp_thread, run_dir, output_dir, input_file, run_dir, file_rotdens, run_dir, run_dir, qmcinp, exe_file, run_dir, run_dir)
        return job_string

def Submission(NameOfServer,status, TransMove, RotMove, RUNDIR, dir_run_job, folder_run, src_dir, execution_file, Rpt, numbbeads, i, step, step_trans, level, temperature, numbblocks, numbpass, molecule_rot, numbmolecules, gfact, dipolemoment, TypeCal, ENT_TYPE, ENT_ALGR, dir_output, dir_run_input_pimc, RUNIN, particleA, NameOfPartition, status_cagepot, iStep, PPA1, user_name, out_dir, source_dir_exe, Restart1,numbblocks_Restart1,crystal,RotorType, spin_isomer,impurity, step_trans_impurity, level_impurity):
        argument1     = Rpt
        step1_trans   = step_trans[iStep]
        level1        = level[iStep]
        step1         = step[iStep]
        step_trans_impurity1   = step_trans_impurity[iStep]
        level_impurity1        = level_impurity[iStep]
        final_dir_in_work = dir_output + folder_run

        if (TypeCal == 'PIGS'):
                argument2     = "pgR"+str(Rpt)+'n'+str(numbmolecules)+"b"
        if (TypeCal == 'PIMC'):
                argument2     = "pm"+str(numbmolecules)+"b"
        if (TypeCal == 'ENT'):
                argument2     = "et"+str(numbmolecules)+"a"+str(particleA)+"b"
        job_name       = argument2+str(i)+".out"

        if not Restart1:
                os.chdir(dir_output)
                if (os.path.isdir(folder_run) == True):
                        print("")
                        print("")
                        print("Error message")
                        print("")
                        print("")
                        printingMessage = "Remove "+str(dir_output)+str(folder_run)
                        print(printingMessage)
                        os.chdir(src_dir)
                        return

                os.chdir(src_dir)
                if (RotorType == "LINEAR"):
                        rotmat(TypeCal,molecule_rot,temperature,numbbeads,source_dir_exe)

                        #call(["rm", "-rf", folder_run])
                
                temperature1    = "%8.6f" % temperature
                if (RotorType == "LINEAR"):
                        file_rotdens    = molecule_rot+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
                        call(["mv", file_rotdens, dir_run_input_pimc])
                        path_Density = dir_run_input_pimc+"/"
                else:
                        iodevn   = spin_isomer
                        jmax = 70
                        if (TypeCal == 'PIMC'):
                                numbbeads2 = int(numbbeads)
                        else:
                                numbbeads2 = int(numbbeads-1)
                        if (molecule_rot == "H2O"):     
                                if (NameOfServer == "graham"):
                                        dir_dens =  "/scratch/"+user_name+"/rot-dens-asymmetric-top/"+asym.GetDirNameSubmission(molecule_rot,temperature1, numbbeads2, iodevn, jmax)
                                if (NameOfServer == "nlogn"):
                                        dir_dens =  "/work/"+user_name+"/rot-dens-asymmetric-top/"+asym.GetDirNameSubmission(molecule_rot,temperature1, numbbeads2, iodevn, jmax)
                        if (molecule_rot == "CH3F"):    
                                if (NameOfServer == "nlogn"):
                                        dir_dens =  "/work/"+user_name+"/rot-dens-symmetric-top/"+sym.GetDirNameSubmission(molecule_rot,temperature1, numbbeads2, 3, jmax)
                                if (NameOfServer == "graham"):
                                        dir_dens =  "/scratch/"+user_name+"/rot-dens-symmetric-top/"+sym.GetDirNameSubmission(molecule_rot,temperature1, numbbeads2, 3, jmax)
                        file_rotdens = "/"+molecule_rot+"_T"+str(dropzeros(temperature1))+"t"+str(numbbeads2)
                        file_rotdens_mod = "/"+molecule_rot+"_T"+str(temperature1)+"t"+str(numbbeads)
                        if (os.path.isfile(dir_dens+file_rotdens_mod+".rho") == False):
                                call(["cp", dir_dens+file_rotdens+".rho", dir_dens+file_rotdens_mod+".rho"])
                                call(["cp", dir_dens+file_rotdens+".eng", dir_dens+file_rotdens_mod+".eng"])
                                call(["cp", dir_dens+file_rotdens+".esq", dir_dens+file_rotdens_mod+".esq"])
                        path_Density = dir_dens+"/"
                        call(["cp","initial_euler_angles_and_com.txt", dir_run_input_pimc])

                if (crystal == True):
                        call(["cp", "LatticeConfig.xyz", dir_run_input_pimc])
        else:
                os.chdir(dir_output)


                #
                if (os.path.isdir(folder_run) == False):
                        print("")
                        print("")
                        print("Error message")
                        print("")
                        print("")
                        printingMessage = dir_output+folder_run+"  --- This directory is absent."
                        print(printingMessage)
                        os.chdir(src_dir)
                        return

                
                #
                logout_file = dir_run_input_pimc+"/"+job_name

                if not "slurmstepd" in open(logout_file).read():
                        if not "real" in open(logout_file).read():
                                printingMessage = dir_output+folder_run+"  --- This job is running."
                                print(printingMessage)
                                os.chdir(src_dir)
                                return


                #
                if ((TypeCal == 'PIGS') or (TypeCal == "PIMC")):
                        file_count = "output.eng"
                if (TypeCal == 'ENT'):
                        file_count = "output.rden"


                #
                os.chdir(dir_output+folder_run+"/results")
                if (os.path.isfile(file_count) == True):
                        col_data_new = genfromtxt(file_count)
                        lastIndex = int(col_data_new[-1,0])
                        if ((numbblocks-lastIndex) <= 0):
                                print(dir_output+folder_run)
                                print(" Done. Do not need to resubmit it.")
                                os.chdir(src_dir)
                                return


                #
                os.chdir(src_dir)

                flag=file_operations(TypeCal,final_dir_in_work,numbmolecules,numbbeads)
                if (flag == False):
                        print(dir_output+folder_run)
                        print("Already resubmitted.")
                        return


                print(dir_output+folder_run)
                print(" Just resubmitted.")


                # Rotational density matrices
                temperature1    = "%8.6f" % temperature
                if (RotorType != "LINEAR"):
                        iodevn   = spin_isomer
                        jmax = 66
                        if (TypeCal == 'PIMC'):
                                numbbeads2 = int(numbbeads)
                        else:
                                numbbeads2 = int(numbbeads-1)
                if (molecule_rot == "H2O"):     
                        if (NameOfServer == "graham"):
                                dir_dens =  "/scratch/"+user_name+"/rot-dens-asymmetric-top/"+asym.GetDirNameSubmission(molecule_rot,temperature1, numbbeads2, iodevn, jmax)
                        if (NameOfServer == "nlogn"):
                                dir_dens =  "/work/"+user_name+"/rot-dens-asymmetric-top/"+asym.GetDirNameSubmission(molecule_rot,temperature1, numbbeads2, iodevn, jmax)
                if (molecule_rot == "CH3F"):    
                        if (NameOfServer == "nlogn"):
                                dir_dens =  "/work/"+user_name+"/rot-dens-symmetric-top/"+sym.GetDirNameSubmission(molecule_rot,temperature1, numbbeads2, 3, jmax)
                        if (NameOfServer == "graham"):
                                dir_dens =  "/scratch/"+user_name+"/rot-dens-symmetric-top/"+sym.GetDirNameSubmission(molecule_rot,temperature1, numbbeads2, 3, jmax)
                file_rotdens = "/"+molecule_rot+"_T"+str(dropzeros(temperature1))+"t"+str(numbbeads2)
                file_rotdens_mod = "/"+molecule_rot+"_T"+str(temperature1)+"t"+str(numbbeads)
                if (os.path.isfile(dir_dens+file_rotdens_mod+".rho") == False):
                        call(["cp", dir_dens+file_rotdens+".rho", dir_dens+file_rotdens_mod+".rho"])
                        call(["cp", dir_dens+file_rotdens+".eng", dir_dens+file_rotdens_mod+".eng"])
                        call(["cp", dir_dens+file_rotdens+".esq", dir_dens+file_rotdens_mod+".esq"])
                path_Density = dir_dens+"/"

                
        # For qmc.imput
        GetInput(TypeCal,ENT_TYPE,ENT_ALGR, temperature,numbbeads,numbblocks,numbpass,molecule_rot,numbmolecules,argument1,level1,step1,step1_trans,gfact,dipolemoment,particleA, Restart1, numbblocks_Restart1, crystal, RotorType, TransMove, RotMove, path_Density,impurity, step_trans_impurity1, level_impurity1)

        input_file    = "qmcbeads"+str(i)+".input"
        call(["mv", "qmc.input", dir_run_input_pimc+"/"+input_file])
        folder_run_path = dir_run_job + folder_run 


        #job submission
        fname         = 'job-for-P'+str(numbbeads)
        fwrite        = open(fname, 'w')

        if RUNDIR == "scratch":
                if RUNIN == "CPU":
                        fwrite.write(jobstring_scratch_cpu(argument2,i,numbmolecules, folder_run_path, molecule_rot, temperature, numbbeads, final_dir_in_work, dir_run_input_pimc, src_dir, Restart1, dir_run_job,status_cagepot, dir_output))
                else:
                        fwrite.write(jobstring_sbatch(NameOfServer, RUNDIR, argument2,i,numbmolecules, folder_run_path, molecule_rot, temperature, numbbeads, final_dir_in_work, dir_run_input_pimc, PPA1, user_name, out_dir, Restart1, dir_run_job,status_cagepot, dir_output, numbblocks))
        else: 
                fwrite.write(jobstring_sbatch(NameOfServer, RUNDIR, argument2, i, numbmolecules, folder_run_path, molecule_rot, temperature, numbbeads, final_dir_in_work, dir_run_input_pimc, PPA1, user_name, out_dir, Restart1, dir_run_job,status_cagepot, dir_output, numbblocks))

        fwrite.close()
        call(["mv", fname, dir_run_input_pimc])
        os.chdir(dir_run_input_pimc)

        if (RUNIN == "CPU"):
                call(["chmod", "755", fname])
                #command_pimc_run = "./"+fname + ">"+ final_dir_in_work+"/outpimc"+str(i)+" & "
                command_pimc_run = "./"+fname + ">outpimc"+str(i)+" & "
                print(command_pimc_run)
                system(command_pimc_run)
        else:
                #call(["qsub", fname])
                if (NameOfPartition == user_name):
                        call(["sbatch", "-p", user_name, fname])
                else:
                        call(["sbatch", fname])

        os.chdir(src_dir)

def jobstring_scratch_cpu(file_name, value, thread, run_dir, molecule, temperature, numbbeads, final_dir, dir_run_input_pimc, src_dir, Restart1, dir_run_job,status_cagepot, dir_output):
        '''
        This function creats jobstring for #PBS script
        '''
        omp_thread     = str(thread)
        output_dir     = run_dir+"/results"
        temperature1   = "%5.3f" % temperature
        file_rotdens   = dir_run_input_pimc+"/"+molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"

        input_file     = dir_run_input_pimc+"/qmc"+file_name+str(value)+".input"
        exe_file       = dir_run_input_pimc+"/pimc"
        qmcinp         = "qmc"+file_name+str(value)+".input"

        job_string     = """#!/bin/bash
export OMP_NUM_THREADS=%s
rm -rf %s
mkdir -p %s
cd %s
mv %s %s
mv %s %s 
cd %s
cp %s qmc.input
cp %s %s
./pimc 
mv %s /work/tapas/linear_rotors
""" % (omp_thread, run_dir, output_dir, src_dir, input_file, run_dir, file_rotdens, run_dir, run_dir, qmcinp, exe_file, run_dir, run_dir)
        return job_string

def jobstring_sbatch(NameOfServer, RUNDIR, file_name, value, numbmolecules, folder_run_path, molecule, temperature, numbbeads, final_dir_in_work, dir_run_input_pimc, PPA1, user_name, out_dir, Restart1, dir_run_job,status_cagepot, dir_output, numbblocks):
        '''
        This function creats jobstring for #SBATCH script
        '''
        if (numbblocks <= 1000):
                walltime   = "00-00:30"
                thread     = 1
        else:
                if (numbbeads >= 160):
                        thread     = 8
                        walltime   = "3-00:00"
                elif ((numbbeads >= 50) and (numbbeads < 160)):
                        thread     = 8
                        walltime   = "03-00:00"
                elif ((numbbeads >= 30) and (numbbeads < 50)):
                        thread     = 4
                        walltime   = "03-00:00"
                else:
                        thread     = 1
                        walltime   = "03-00:00"
        
        job_name       = file_name+str(value)
        omp_thread     = str(thread)
        output_dir     = folder_run_path+"/results"
        temperature1   = "%8.6f" % temperature
        file_rotdens   = dir_run_input_pimc+"/"+molecule+"_T"+str(temperature1)+"t"+str(numbbeads)+".*"
        logpath        = dir_run_input_pimc+"/"+job_name

        input_file     = dir_run_input_pimc+"/qmcbeads"+str(value)+".input"
        exe_file       = dir_run_input_pimc+"/pimc"
        input_file1    = dir_run_input_pimc+"/initial_euler_angles_and_com.txt"
        qmcinp         = "qmcbeads"+str(value)+".input"

        if (status_cagepot == True):
                cagepot_file   = dir_run_input_pimc+"/hfc60.pot"
                cagepot_cp     = "cp "+cagepot_file+" "+folder_run_path
        else:
                cagepot_cp     = ""

        if NameOfServer == "graham":
                CommandForMove = " " 
        else:
                if (RUNDIR == "scratch"):
                        CommandForMove = "mv "+folder_run_path+" "+dir_output
                if (RUNDIR == "work"):
                        CommandForMove = " "

        if NameOfServer == "graham":
                account = "#SBATCH --account=rrg-pnroy"
        else:
                account = ""

        if not PPA1:
                CommandForPPA = "#"
        else:
                CommandForPPA = ""
        file_PPA       = dir_run_input_pimc+"/PairDensity.txt"

        print("")
        print("")
        print("")
        print("*****************Important Notice***********************")
        print("")
        print("Full path of the directory where the submitted job is running - ")
        print(folder_run_path)
        print("")
        print("Full path of the directory where all the outputs of MoRiBs are stored - ")
        print(output_dir)
        print("")
        print("Informations about all types of acceptance ratios of Monte Carlo simulations are saved at - ")
        print(logpath+".out")
        print("")
        print("Number of OpenMP thread used = "+str(thread))
        print("")
        print("After the job is completed successfully, the directory where the job was running is moved to - ")
        print(final_dir_in_work)
        print("")
        print("********************************************************")
        print("")
        print("")
        print("")


        job_string     = """#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --output=%s.out
#SBATCH --time=%s
%s
#SBATCH --constraint=broadwell
#SBATCH --mem-per-cpu=600mb
#SBATCH --cpus-per-task=%s
export OMP_NUM_THREADS=%s
export OMP_STACKSIZE=600M
export GOMP_STACKSIZE=600M
rm -rf %s
mkdir -p %s
mv %s %s
cd %s
cp %s qmc.input
cp %s %s
cp %s %s
%s cp %s %s
#valgrind --tool=memcheck --leak-check=yes --show-reachable=yes ./pimc
time ./pimc 
%s
""" % (job_name, logpath, walltime, account, omp_thread, omp_thread, folder_run_path, output_dir, input_file, folder_run_path, folder_run_path, qmcinp, exe_file, folder_run_path, input_file1, folder_run_path, CommandForPPA, file_PPA, folder_run_path, CommandForMove)

        job_string_restart = """#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --output=%s.out
#SBATCH --time=%s
%s
#SBATCH --mem-per-cpu=8192mb
#SBATCH --cpus-per-task=%s
export OMP_NUM_THREADS=%s
mv %s %s
mv %s %s
cd %s
cp %s qmc.input
cp %s %s
####valgrind --leak-check=full -v --show-leak-kinds=all ./pimc 
time ./pimc 
%s
""" % (job_name, logpath, walltime, account, omp_thread, omp_thread,final_dir_in_work, dir_run_job, input_file, folder_run_path, folder_run_path, qmcinp, exe_file, folder_run_path, CommandForMove)

        if Restart1:
                return job_string_restart
        else:
                return job_string

def GetRotEnergy(molecule,jrot):
        Energy = GetBconst(molecule)*jrot*(jrot+1.0)
        return Energy

def GetAvgRotEnergy(molecule,beta):
        CMRECIP2KL = 1.4387672
        Zsum = 0.0
        Nsum = 0.0
        for jrot in range(0,10000,1):
                BoltzmannProb = exp(-beta*GetRotEnergy(molecule,jrot)*CMRECIP2KL)
                if (BoltzmannProb > 10e-16):
                        Zsum += (2*jrot+1.0)*BoltzmannProb
                        Nsum += (2*jrot+1.0)*GetRotEnergy(molecule,jrot)*BoltzmannProb
                else:
                        break
        AvgEnergy = Nsum*CMRECIP2KL/Zsum
        return AvgEnergy

def GetFileNameSubmission(TypeCal, molecule_rot, TransMove, RotMove, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules, molecule, ENT_TYPE, particleA, extra, crystal,impurity,ext_ent):
        #add                     = "-NumTimes"
        add                     = ""
        if (TypeCal == "ENT"):
                add1                = "-ParticleA"+str(particleA)
                add2                = "-"
        else:
                if (len(impurity) == 4):
                        add1            = "-"+str(impurity[2])+"-"+impurity[1]
                else:
                        add1            = ""
                add2                = ""

        if (crystal == True):
                mainFileName        = parameterName+str(parameter)+"Kinv-Blocks"+str(numbblocks)+"-Passes"+str(numbpass)+add+"-System"+str(molecule)+add1+"-e0vsbeads"+add2 
        else:
                mainFileName        = parameterName+str(parameter)+"Kinv-Blocks"+str(numbblocks)+"-Passes"+str(numbpass)+add+"-System"+str(numbmolecules)+str(molecule)+add1+"-e0vsbeads"+add2 

        if (TypeCal == "PIGS"):
                frontName           = "PIGS-"
        if (TypeCal == "PIMC"):
                frontName           = "PIMC-"
        if (TypeCal == "ENT"):
                frontName           = "ENT-"

        frontName              += extra

        if ((TransMove == True) and (RotMove == True)):
                frontName      += "TransAndRotDOFs-"
        if ((TransMove == False) and (RotMove == True)):
                frontName      += "RotDOFs-"
        if ((TransMove == True) and (RotMove == False)):
                frontName      += "TransDOFs-"


        if (Rpt >= 0.0):
                FragmentRpt = "Rpt"+str(Rpt)+"Angstrom-"
        else:
                FragmentRpt = ""

        if (dipolemoment >= 0.0):
                if (numbmolecules > 1):
                        FragmentDipoleMoment = "DipoleMoment"+str(dipolemoment)+"Debye-"
                else:
                        FragmentDipoleMoment = "Field"+str(dipolemoment)+"Kinv-"
        else:
                FragmentDipoleMoment = ""

        if (gfact >= 0.0):
                FragmentGFactor = "gFactor"+str(gfact)+"-"
        else:
                FragmentGFactor = ""

        file1_name      = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+mainFileName
        if (TypeCal == "ENT"):
                file1_name += ENT_TYPE+"-"+ext_ent+"-"
        
        return file1_name

class GetFileNameAnalysis:
        def __init__(self, TypeCal1,TypeCal2, molecule_rot1, TransMove1, RotMove1, variableName1, Rpt1, gfact1, dipolemoment1, parameterName1, parameter1, numbblocks1, numbpass1, numbmolecules1, molecule1, ENT_TYPE1, preskip1, postskip1, extra1, src_dir1, particleA1, ENT_ALGR):
                self.TypeCal      = TypeCal1
                self.molecule_rot = molecule_rot1
                self.TransMove    = TransMove1
                self.RotMove      = RotMove1
                self.variableName = variableName1
                self.Rpt          = Rpt1
                self.dipolemoment = dipolemoment1
                self.gfact        = gfact1
                self.parameter    = parameter1
                self.parameterName= parameterName1
                self.numbblocks   = numbblocks1
                self.numbpass     = numbpass1
                self.numbmolecules= numbmolecules1
                self.molecule     = molecule1
                self.ENT_TYPE     = ENT_TYPE1
                self.preskip      = preskip1
                self.postskip     = postskip1
                self.extra        = extra1
                self.src_dir      = src_dir1
                self.particleA    = particleA1

                if (self.TypeCal == "ENT"):
                        frontName             = "ENT-"+self.extra
                        add1                  = "-ParticleA"+str(self.particleA)
                        add2                  = "-"+self.ENT_TYPE+"-"+ENT_ALGR
                else:
                        frontName             = self.TypeCal+"-"+self.extra
                        add1                  = ""
                        add2                  = ""

                if ((self.TransMove == True) and (self.RotMove == True)):
                        frontName            += "TransAndRotDOFs-"
                if ((self.TransMove == False) and (self.RotMove == True)):
                        frontName            += "RotDOFs-"
                if ((self.TransMove == True) and (self.RotMove == False)):
                        frontName            += "TransDOFs-"

                if (self.Rpt >= 0.0):
                        FragmentRpt           = "Rpt"+str(self.Rpt)+"Angstrom-"
                else:
                        FragmentRpt           = ""

                if (self.dipolemoment >= 0.0):
                        if (self.numbmolecules > 1):
                                FragmentDipoleMoment = "DipoleMoment"+str(self.dipolemoment)+"Debye-"
                        else:
                                FragmentDipoleMoment = "Field"+str(self.dipolemoment)+"Kinv-"
                else:
                        FragmentDipoleMoment  = ""

                if (self.gfact >= 0.0):
                        FragmentGFactor       = "gFactor"+str(self.gfact)+"-"
                else:
                        FragmentGFactor       = ""

                mainFileName  = "vs-"+str(self.variableName)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-Blocks"+str(self.numbblocks)
                mainFileName += "-Passes"+str(self.numbpass)+"-System"+str(self.numbmolecules)+str(self.molecule)+add1+add2+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)

                file_output1  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"Energy-"    
                file_output2  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"correlation-"
                file_output3  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"total-correlation-function-"
                file_output4  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"X-component-correlation-function-"
                file_output5  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"Y-component-correlation-function-"
                file_output6  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"Z-component-correlation-function-"
                file_output7  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"XandY-component-correlation-function-"
                file_output8  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"Entropy-"

                self.SaveEnergy       = self.src_dir+file_output1+mainFileName+".txt"
                self.SaveCorr         = self.src_dir+file_output2+mainFileName+".txt"
                self.SaveEntropy      = self.src_dir+file_output8+mainFileName+".txt"

                if (TypeCal2 == False):
                        if os.path.exists(self.SaveEntropy):   os.remove(self.SaveEntropy)
                        if os.path.exists(self.SaveEnergy):    os.remove(self.SaveEnergy)
                        if os.path.exists(self.SaveCorr):      os.remove(self.SaveCorr)

                if (self.TypeCal != "ENT"):
                        print("#-------------------------------------#")
                        print("Final analyzed results are stored in - ")
                        print(self.src_dir)
                        print("")
                        print("Final results - Energy vs "+str(self.variableName))
                        print(file_output1+mainFileName+".txt")
                        print(file_output2+mainFileName+".txt")
                        print("#------------------------------------------------------------------------#")

                if (self.TypeCal == "ENT"):
                        mainFileNameRT    = "vs-"+str(self.variableName)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-Blocks"+str(self.numbblocks)
                        mainFileNameRT   += "-Passes"+str(self.numbpass)+"-System"+str(self.numbmolecules)+str(self.molecule)+add1+add2+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)

                        self.SaveEntropyRT= self.src_dir+file_output8+mainFileNameRT+".txt"
                        if (ENT_ALGR != "WR"):
                                #print(self.SaveEntropyRT)
                                print(self.src_dir)
                                print("")
                                print("Final results - Entropy vs "+str(self.variableName))
                                print(self.SaveEntropy)
                                
                                print("#------------------------------------------------------------------------#")

class GetFileNamePlot:
        def __init__(self, TypeCal1, molecule_rot1, TransMove1, RotMove1, variableName1, Rpt1, gfact1, dipolemoment1, parameterName1, parameter1, numbblocks1, numbpass1, numbmolecules1, molecule1, ENT_TYPE1, preskip1, postskip1, extra1, src_dir1, particleA1, var1):
                self.TypeCal      = TypeCal1
                self.molecule_rot = molecule_rot1
                self.TransMove    = TransMove1
                self.RotMove      = RotMove1
                self.variableName = variableName1
                self.var          = var1
                self.Rpt          = Rpt1
                self.dipolemoment = dipolemoment1
                self.gfact        = gfact1
                self.parameter    = parameter1
                self.parameterName= parameterName1
                self.numbblocks   = numbblocks1
                self.numbpass     = numbpass1
                self.numbmolecules= numbmolecules1
                self.molecule     = molecule1
                self.ENT_TYPE     = ENT_TYPE1
                self.preskip      = preskip1
                self.postskip     = postskip1
                self.extra        = extra1
                self.src_dir      = src_dir1
                self.particleA    = particleA1

                if (self.TypeCal == "ENT"):
                        frontName             = "ENT-"+self.extra
                        add1                  = "-ParticleA"+str(self.particleA)
                        if (self.ENT_TYPE):     
                                add2                  = "-"+self.ENT_TYPE
                        else:
                                add2                  = ""
                else:
                        frontName             = self.TypeCal+"-"+self.extra
                        add1                  = ""
                        add2                  = ""

                if ((self.TransMove == True) and (self.RotMove == True)):
                        frontName            += "TransAndRotDOFs-"
                if ((self.TransMove == False) and (self.RotMove == True)):
                        frontName            += "RotDOFs-"
                if ((self.TransMove == True) and (self.RotMove == False)):
                        frontName            += "TransDOFs-"

                if (self.Rpt >= 0.0):
                        FragmentRpt           = "Rpt"+str(self.Rpt)+"Angstrom-"
                else:
                        FragmentRpt           = ""

                if (self.dipolemoment >= 0.0):
                        FragmentDipoleMoment  = "DipoleMoment"+str(self.dipolemoment)+"Debye-"
                else:
                        FragmentDipoleMoment  = ""

                if (self.gfact >= 0.0):
                        FragmentGFactor       = "gFactor"+str(self.gfact)+"-"
                else:
                        FragmentGFactor       = ""

                mainFileName  = "vs-"+str(self.variableName)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-Blocks"+str(self.numbblocks)
                mainFileName += "-Passes"+str(self.numbpass)+"-System"+str(self.numbmolecules)+str(self.molecule)+add1+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)+add2

                file_output1  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"Energy-"
                file_output2  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"correlation-"
                file_output3  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"total-correlation-function-"
                file_output4  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"X-component-correlation-function-"
                file_output5  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"Y-component-correlation-function-"
                file_output6  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"Z-component-correlation-function-"
                file_output7  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"XandY-component-correlation-function-"
                file_output8  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"Chemical-Potential-"
                file_output9  = frontName+FragmentRpt+FragmentDipoleMoment+FragmentGFactor+"Entropy-"
                file_output10 = frontName+FragmentDipoleMoment+FragmentGFactor+"Energy-vs-R-"
                file_output11 = frontName+FragmentDipoleMoment+FragmentGFactor+"Order-parameters-vs-R-"
                mainFileNameFitvsR  = "fixed-"+self.parameterName+str(self.parameter)+"Kinv-Blocks"+str(self.numbblocks)
                mainFileNameFitvsR += "-Passes"+str(self.numbpass)+"-System"+str(self.numbmolecules)+str(self.molecule)+add1+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)+add2

                self.SaveEnergy       = self.src_dir+file_output1+mainFileName
                self.SaveCorr         = self.src_dir+file_output2+mainFileName
                self.SaveChemPot      = self.src_dir+file_output8+mainFileName
                self.SaveEntropy      = self.src_dir+file_output9+mainFileName
                self.SaveEnergyFitvsR = self.src_dir+file_output10+mainFileNameFitvsR
                self.SaveCorrFitvsR   = self.src_dir+file_output11+mainFileNameFitvsR

#---------------------------------------------------------------------------#
#       special cases                                                           #
#---------------------------------------------------------------------------#
                '''
                mainFileNameCP        = "vs-number-of-"+str(self.molecule)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-Blocks"+str(self.numbblocks)
                mainFileNameCP       += "-Passes"+str(self.numbpass)+add1+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)
                mainFileNameCONV      = "vs-beta-and-tau-Blocks"+str(self.numbblocks)
                mainFileNameCONV     += "-Passes"+str(self.numbpass)+"-System"+str(self.numbmolecules)+str(self.molecule)+add1+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)
                self.SaveChemPot      = self.src_dir+file_output8+mainFileNameCP
                self.SaveEntropyCONV  = self.src_dir+"/ResultsOfPIGSENT/"+file_output1+mainFileNameCONV+"-"+self.ENT_TYPE

                '''
#
                mainFileNameGFAC      = "vs-gFactor-of-"+str(self.molecule)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-numbbeads"+str(self.var)+"-Blocks"+str(self.numbblocks)
                mainFileNameGFAC     += "-Passes"+str(self.numbpass)+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)+add2
                mainFileNameGFACFit   = "vs-gFactor-of-"+str(self.molecule)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-Blocks"+str(self.numbblocks)
                mainFileNameGFACFit  += "-Passes"+str(self.numbpass)+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)+add2+"-Fit"
#
                mainFileNameRFAC      = "vs-RFactor-of-"+str(self.molecule)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-numbbeads"+str(self.var)+"-Blocks"+str(self.numbblocks)
                mainFileNameRFAC     += "-Passes"+str(self.numbpass)+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)+add2
                mainFileNameMM            = "vs-"+str(self.variableName)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv"
                mainFileNameMM       += "-System"+str(self.numbmolecules)+str(self.molecule)+add1
                mainFileNameCOMBO     = "vs-gFactor-fixed-"+self.parameterName+str(self.parameter)+"Kinv-Blocks"+str(self.numbblocks)
                mainFileNameCOMBO    += "-Passes"+str(self.numbpass)+"-System"+str(self.numbmolecules)+str(self.molecule)+add1+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)
                mainFileNameED            = "of-System"+str(self.numbmolecules)+str(self.molecule)+add1
#
                self.SaveEntropyGFAC  = self.src_dir+frontName+FragmentRpt+"Entropy-"+mainFileNameGFAC
                self.SaveEntropyRFAC  = self.src_dir+frontName+FragmentRpt+"Entropy-"+mainFileNameRFAC
                self.SaveEnergyGFAC   = self.src_dir+frontName+FragmentRpt+"Energy-"+mainFileNameGFAC
                self.SaveEnergyRFAC   = self.src_dir+frontName+FragmentRpt+"Energy-"+mainFileNameRFAC
                self.SaveEnergyED     = self.src_dir+file_output1+mainFileNameED+add2+"-ED"
                self.SaveEntropyED    = self.src_dir+file_output9+mainFileNameED+add2+"-ED"
                self.SaveCorrGFAC     = self.src_dir+frontName+FragmentRpt+"correlation-"+mainFileNameGFAC
                self.SaveCorrRFAC     = self.src_dir+frontName+FragmentRpt+"correlation-"+mainFileNameRFAC
                self.SaveEnergyMM     = self.src_dir+file_output1+mainFileNameMM+add2+"-MM"
                self.SaveEntropyMM    = self.src_dir+file_output9+mainFileNameMM+add2+"-MM"
                self.SaveEntropyCOMBO = self.src_dir+frontName+FragmentRpt+"Entropy-"+mainFileNameGFACFit+"-COMBINE"
                self.SaveEntropyGFACFit = self.src_dir+frontName+FragmentRpt+"Entropy-"+mainFileNameGFACFit
                self.SaveEnergyGFACFit= self.src_dir+frontName+FragmentRpt+"Energy-"+mainFileNameGFACFit

                if (self.TypeCal == "ENT"):
                        mainFileNameRT    = "vs-"+str(self.variableName)+"-fixed-"+self.parameterName+str(self.parameter)+"Kinv-Blocks"+str(self.numbblocks)
                        mainFileNameRT   += "-Passes"+str(self.numbpass)+"-System"+str(self.numbmolecules)+str(self.molecule)+"-preskip"+str(self.preskip)+"-postskip"+str(self.postskip)+add2

                        self.SaveEntropyRT     = self.src_dir+file_output9+mainFileNameRT

def FileCheck(TypeCal,list_nb,variableName,SavedFile):
        for i in list_nb:
                if (TypeCal == "PIMC"):
                        if ((i%2) == 0):
                                bead = i
                        else:
                                bead = i+1
                else:
                        if ((i%2) != 0):
                                bead = i
                        else:
                                bead = i+1

                string = str(bead)
                if string in open(SavedFile).read():
                        return

        if not string in open(SavedFile).read():
                call(["rm", SavedFile])

class GetUnitConverter:
        def __init__(self):
                self.BOHRRADIUS = 0.5291772108;      # angstrom
                self.HARTREE2JL = 4.359748e-18;         # hartree to joule  conversion factor
                self.HARTREE2KL = 3.157732e+05;         # hartree to Kelvin conversion factor
                self.CMRECIP2KL = 1.4387672;            # cm^-1 to Kelvin conversion factor
                self.MHZ2RCM    = 3.335640952e-5;       # MHz to cm^-1 conversion factor

                self.AuToDebye     = 1.0/0.39343;
                self.AuToCmInverse = 219474.63137;
                self.AuToKelvin    = 315777.0;
                self.KCalperMolToCmInverse = 349.75509;

                self.HBAR  = 1.05457266;                        #  (10^-34 Js)     Planck constant
                self.AMU   = 1.6605402;                         #  (10^-27 kg)     atomic mass unit
                self.K_B   = 1.380658;                          #  (10^-23 JK^-1)  Boltzmann constant
                self.WNO2K = 0.6950356;                                 # conversion from CM-1 to K
                self.kcalmoleinvToKelvin = 503.228

def GetrAndgFactor(molecule, RCOM, DipoleMoment):
        '''
        It calculates g and R value 
        '''
        Units          = GetUnitConverter()
        BConstant      = GetBconst(molecule)  # in wavenumber
        DipoleMomentAU = DipoleMoment/Units.AuToDebye
        RCOMAU         = RCOM/Units.BOHRRADIUS
        BConstantAU    = BConstant/Units.AuToCmInverse
        rFactor        = RCOMAU/((DipoleMomentAU*DipoleMomentAU/BConstantAU)**(1.0/3.0))
        gFactor        = (DipoleMomentAU*DipoleMomentAU)/(RCOMAU*RCOMAU*RCOMAU*BConstantAU)
        printingmessage   = " DipoleMoment = "+str(DipoleMoment)+" gFactor = " + str(gFactor)+ " rFactor = "+str(rFactor)
        print(printingmessage)
        returnList = [rFactor, gFactor]
        return returnList

def GetDipoleMomentFromGFactor(molecule, RCOM, gFactor):
        '''
        It extracts dipole moment from a g value - g 
        '''
        Units          = GetUnitConverter()
        BConstant      = GetBconst(molecule)  # in wavenumber
        RCOMAU         = RCOM/Units.BOHRRADIUS
        BConstantAU    = BConstant/Units.AuToCmInverse
        DipoleMomentAU = sqrt(gFactor*RCOMAU*RCOMAU*RCOMAU*BConstantAU)
        DipoleMoment   = DipoleMomentAU*Units.AuToDebye
        return DipoleMoment

def GetEDResults(TypeCal, FilePlotName, srcCodePath, RFactor, numbmolecules, particleA, lmax, ltotalmax):
        '''
        It will give ground state energy, von Neuman and Renyi entropies computed by diagonalizing full Hamiltonian matrix. It is developed by Dmitri https://github.com/0/DipoleChain.jl
        '''
        if (TypeCal == "ENT"):
                FileToBeSavedED = FilePlotName.SaveEntropyED+".txt"
                FileToBeSavedMM = FilePlotName.SaveEntropyMM+".txt"
        if (TypeCal == "PIGS"):
                FileToBeSavedED = FilePlotName.SaveEnergyED+".txt"
                FileToBeSavedMM = FilePlotName.SaveEnergyMM+".txt"
        print(FileToBeSavedED)

        commandRunED    = "julia "+srcCodePath+"diagonalization.jl -R "+str(RFactor)+" -N "+str(numbmolecules)+" --l-max "+str(lmax)+" --l-total-max "+str(ltotalmax)+" --A-start 1 --A-size "+str(particleA)
        call(["rm", "outputED.txt"])
        system(commandRunED)
        call(["mv", "outputED.txt", FileToBeSavedED])
        '''
        if (numbmolecules >6):
                print("It is not computing Matrix multiplication stuffs")
                return

        if (numbmolecules <= 4):
                call(["rm", "outputMM.txt"])

                for numbbeads in loop:
                        print(numbbeads)
                        RFactor      = GetrAndgFactor(molecule_rot, Rpt, dipolemoment)
                        if (variableName == "beta"):
                                parameterR    = parameter*BConstantK
                                commandRun   = "julia "+srcCodePath+"path_integral.jl -R "+str(RFactor)+" -N "+str(numbmolecules)+" --l-max 6 --tau "+str(parameterR)+" -P "+str(numbbeads)+" --pigs --A-start 1"+" --A-size "+str(particleA)
                        if (variableName == "tau"):
                                parameterR    = parameter*BConstantK
                                commandRun   = "julia "+srcCodePath+"path_integral.jl -R "+str(RFactor)+" -N "+str(numbmolecules)+" --l-max 6 --beta "+str(parameterR)+" -P "+str(numbbeads)+" --pigs --A-start 1"+" --A-size "+str(particleA)
                        system(commandRun)
                call(["mv", "outputMM.txt", FileToBeSavedMM])
        '''

def GetPairDensity(FilePlotName, srcCodePath, RFactor, numbmolecules, loop, particleA, molecule_rot, Rpt, dipolemoment, parameter, BConstantK, variableName, TypeCal):
        FileToBeSavedDensity = FilePlotName+".txt"
        for numbbeads in loop:
                print(numbbeads)
                parameterR    = parameter*BConstantK
                commandRun    = "julia "+srcCodePath+"pair_density.jl -R "+str(RFactor)+" --l-max 2 --l-total-max 2 --tau "+str(parameterR)
                print(commandRun)
                call(["rm", "outputDensity.txt"])
                system(commandRun)
                call(["mv", "outputDensity.txt", FileToBeSavedDensity])

#def GetEntropyRT(status, maxloop, TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA, variable):
        #FileAnalysis = GetFileNameAnalysis(TypeCal, True, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA)
#       GetAverageEntropyRT(maxloop, TypeCal, molecule_rot, TransMove, RotMove, variableName, Rpt, gfact, dipolemoment, parameterName, parameter, numbblocks, numbpass, numbmolecules1, molecule, ENT_TYPE, preskip, postskip, extra_file_name, final_results_path, particleA, variable)

def GetPreFactDDPot(molecule, RCOM, DipoleMoment):
        '''
        It calculates the pre factor in inverse temperature of dipole - dipole interaction potential
        '''
        Units          = GetUnitConverter()
        DipoleMomentAU = DipoleMoment/Units.AuToDebye
        RCOMAU         = RCOM/Units.BOHRRADIUS
        preFact        = (DipoleMomentAU*DipoleMomentAU)/(RCOMAU*RCOMAU*RCOMAU)
        preFact        = preFact*Units.HARTREE2KL
        printingmessage= " DipoleMoment = "+str(DipoleMoment)+" Debye and the prefactor of the dipole-dipole interaction potential = " + str(preFact)+ " K^-1"
        print(printingmessage)

def GetRenamingFunc(dir_run_input_pimc, dir_input_pimc_renamed, dir_output, folder_run, folder_renamed, src_dir):
        #final_dir_in_work = dir_output + folder_run

        call(["mkdir", "-p", dir_input_pimc_renamed])
        cmd_run= "cp -r "+dir_run_input_pimc+"/*  "+dir_input_pimc_renamed+"/"
        os.system(cmd_run)
        os.chdir(dir_output)
        if (os.path.isdir(folder_run) == True):
                call(["mv", folder_run, folder_renamed])
                print(dir_input_pimc_renamed)
                printingMessage = "move "+str(dir_output)+str(folder_run)
                print(printingMessage)
        os.chdir(src_dir)

def RemoveFiles(TypeCal, numbbeads, temperature, molecule_rot, RotorType, preskip, postskip, numbblocks, final_dir_in_work):

        col_block = genfromtxt(final_dir_in_work+"/results/output.eng",unpack=True, usecols=[0], skip_header=preskip, skip_footer=postskip)
        if (int(len(col_block)) == numbblocks-(preskip+postskip)):
        
                temperature1 = "%8.6f" % temperature
                if (RotorType == "LINEAR"):
                        file_rotdens = molecule_rot+"_T"+str(temperature1)+"t"+str(numbbeads)+".rot"
                        call(["rm", final_dir_in_work+"/"+file_rotdens])
                else:
                        if (TypeCal == 'PIMC'):
                                numbbeads2 = int(numbbeads)
                        else:
                                numbbeads2 = int(numbbeads-1)
                        file_rotdens_mod = molecule_rot+"_T"+str(temperature1)+"t"+str(numbbeads)
                        if (os.path.exists(final_dir_in_work+"/"+file_rotdens_mod+".rho") == True):
                                call(["rm", final_dir_in_work+"/"+file_rotdens_mod+".rho"])
                                call(["rm", final_dir_in_work+"/"+file_rotdens_mod+".eng"])
                                call(["rm", final_dir_in_work+"/"+file_rotdens_mod+".esq"])
