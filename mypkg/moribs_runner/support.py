import time
from subprocess import call
from os import system
import os
import decimal
import numpy as np
import math
from pathlib import Path
import get_beads_and_mc_steps as mc
import mypkg.asymrho_runner.support as asym
import mypkg.symrho_runner.support as sym


def error_message(number):
	if (number <= 1):
		print("Warning!!!")
		print("Increase the number of blocks in your computation!!!")
		exit(0)


def error_propagation(mean, data):
	ndim = len(data)
	error = np.std(data, ddof=0) / math.sqrt(ndim)
	return error


def get_error(mean, data, binary_exponent):
	error_message(binary_exponent)
	error = np.zeros(binary_exponent)
	error[0] = error_propagation(mean, data)
	for i in range(1, binary_exponent):
		ndim = int(len(data) / 2)
		data1 = np.zeros(ndim)

		for j in range(ndim):
			data1[j] = 0.5 * (data[2 * j] + data[2 * j + 1])

		data = data1
		error[i] = error_propagation(mean, data)
	return np.max(error)


def get_execution_file(input_dir, method, ent_method, source_dir_exe):
	execution_file_dir = source_dir_exe
	os.chdir(execution_file_dir)
	call(["make", "-f", "Makefile-GNU", "clean"])
	call(["make", "-f", "Makefile-GNU"])
	path = Path(execution_file_dir + "pimc")

	print("")
	print("")
	print("****************** Important message *******************")
	print("")
	if path.is_file():
		print(" ---- The compilation is successful. ----")
		print("")
		print("********************************************************")
	else:
		print(f' ---- The file {execution_file_dir}pimc does not exist ---- ')
		os.chdir(input_dir)
		print("********************************************************")
		exit()
	os.chdir(input_dir)


def compile_rotmat(source_dir_exe, input_dir):
	path_enter_linden = source_dir_exe + "linear_prop/"
	os.chdir(path_enter_linden)
	call(["make", "clean"])
	call(["make"])
	path_exit_linden = input_dir
	os.chdir(path_exit_linden)


def compile_cagepot(source_dir_exe, input_dir):
	path_enter_cagepot = source_dir_exe + "tabulated_potential/"
	os.chdir(path_enter_cagepot)
	call(["make", "clean"])
	call(["make"])
	path_exit_cagepot = source_dir_exe + input_dir
	os.chdir(path_exit_cagepot)


def jack_knife(mean, data):
	ai = [((np.sum(data) - data[j]) / (len(data) - 1.0))
		  for j in range(len(data))]
	deviation = ai - mean
	devsq = np.multiply(deviation, deviation)
	error = sqrt(np.sum(devsq) * (len(data) - 1.0) / len(data))
	return error


def get_rotational_bconstant(rotor):
	'''
	The unit is cm inverse.
	'''
	if (rotor == "HF"):
		bconst = 20.559 # https://opg.optica.org/josa/abstract.cfm?URI=josa-54-1-20
	if (rotor == "H2"):
		bconst = 60.853
	return bconst


def replace(string_old, string_new, file1, file2):
	'''
	This function replaces old string by new string
	'''
	f1 = open(file1, 'r')
	f2 = open(file2, 'w')
	for line in f1:
		f2.write(line.replace(string_old, string_new))
	f1.close()
	f2.close()


def beads(tau, beta):
	'''
	This function determins number of beads
	'''
	numb_bead1 = beta / tau + 1
	numb_bead2 = int(round(numb_bead1, 0))
	if (numb_bead2 % 2 == 0):
		numb_bead2 = numb_bead2 + 1
	return numb_bead2


def file_operations(method, final_dir_in_work, numb_molecule, numb_bead):

	if (method == "ENT"):
		fileList = ["output.rden", "output.xyz"]
		file_old = final_dir_in_work + "/results/output.rden_old"
	else:
		fileList = ["output.eng", "output.xyz"]
		file_old = final_dir_in_work + "/results/output.eng_old"

	flag = False
	if (os.path.isfile(final_dir_in_work + "/results/" + fileList[0])):
		flag = True
		if (os.path.isfile(file_old)):
			for filecat in fileList:
				if (filecat != "output.xyz"):
					col_data_new = np.genfromtxt(
						final_dir_in_work + "/results/" + filecat)
					index = int(col_data_new[0, 0])
					col_data_old = np.genfromtxt(
						final_dir_in_work + "/results/" + filecat + "_old")
					merged_data = np.concatenate(
						(col_data_old[:index - 1], col_data_new), axis=0)

					if (filecat == "output.eng"):
						np.savetxt(
							final_dir_in_work +
							"/results/" +
							filecat +
							"_old",
							merged_data,
							fmt='%d	%.6e	%.6e	%.6e	%.6e')
					else:
						np.savetxt(
							final_dir_in_work +
							"/results/" +
							filecat +
							"_old",
							merged_data,
							fmt='%.6e',
							delimiter='	')

				if (filecat == "output.xyz"):
					if "H2O1" in open(
							final_dir_in_work + "/results/" + filecat).read():
						rmstr = int(numb_molecule * numb_bead + 3)
						file_temp = final_dir_in_work + "/results/" + filecat + "_temp"
						cmd1 = "tail -n +" + \
							str(rmstr) + " " + final_dir_in_work + \
							"/results/" + filecat + ">" + file_temp
						os.system(cmd1)
						col_data_new = np.genfromtxt(file_temp)
						call(["rm", file_temp])
					else:
						col_data_new = np.genfromtxt(
							final_dir_in_work + "/results/" + filecat)
					index = int(col_data_new[0, 0])
					col_data_old = np.genfromtxt(
						final_dir_in_work + "/results/" + filecat + "_old")
					merged_data = np.concatenate(
						(col_data_old[:index - 1], col_data_new), axis=0)
					np.savetxt(final_dir_in_work + "/results/" + filecat +
							   "_old", merged_data, fmt='%.6e', delimiter='	')

				call(["rm", final_dir_in_work + "/results/" + filecat])
		else:
			for filemv in fileList:
				call(["mv", final_dir_in_work + "/results/" + filemv,
					  final_dir_in_work + "/results/" + filemv + "_old"])
	return flag


def fmt_energy_data(method, parameter_name):
	'''
	This function produces formatted output.
	'''
	if (parameter_name == "beta"):
		variable_name = "tau"
	if (parameter_name == "tau"):
		variable_name = "beta"
	unit = "(1/K)"

	print("")
	output = "# The units of all kinds of energies are in Kelvin."
	output += "\n"
	output += "# "
	if (method == "PIMC"):
		output += '{blocks:^10}{beads:^10}{var:^10}{rot:^16}{pot:^16}{tot:^16}{rotsq:^16}{cv:^16}{er1:^12}{er2:^12}{er3:^12}{er4:^12}{er5:^12}'.format(
			blocks='nBlocks',
			beads='nBeads',
			var=variable_name +
			' invK',
			rot='<K>',
			pot='<V>',
			tot='<E>',
			rotsq='<Ksq>',
			cv='<Cv>',
			er1='Err-K',
			er2='Err-V',
			er3='Err-E',
			er4='Err-Ksq',
			er5='Err-Cv',
		)
		output += "\n"
		output += '{0:=<170}'.format('#')
		output += "\n"

	if (method == "PIGS"):
		output += '{blocks:^10}{beads:^10}{var:^10}{pot:^16}{tot:^16}{er1:^12}{er2:^12}'.format(
			blocks='nBlocks',
			beads='nBeads',
			var=variable_name +
			' invK',
			pot='<V>',
			tot='<E>',
			er1='Err-V',
			er2='Err-E')
		output += "\n"
		output += '{0:=<90}'.format('#')
		output += "\n"

	if (method == "ENT"):
		output += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format(
			'Beads',
			variable,
			'Avg. rotational',
			'Avg. (E - V)',
			'Avg. Potential',
			'Avg. Total',
			'Error of Rotational',
			'Error of (E - V)',
			'Error of Potential',
			'Error of Total')
		output += "\n"
		output += '{0:=<90}'.format('#')
		output += "\n"

	return output


def get_average_energy(
		method,
		numb_bead,
		parameter_name,
		parameter_value,
		final_dir_in_work,
		preskip,
		postskip,
		numb_block):
	'''
	This function gives us the output
	'''
	if (parameter_name == "beta"):
		variable_value = parameter_value/(numb_bead-1)
	if (parameter_name == "tau"):
		variable_value = parameter_value*(numb_bead-1)

	if (os.path.isdir(final_dir_in_work)):
		condition = True

		file_old = final_dir_in_work + "/results/output.eng_old"
		if (os.path.isfile(file_old)):
			file_old_1 = final_dir_in_work + "/results/output.eng_old_1"
			file_new = final_dir_in_work + "/results/output.eng"
			if (os.path.isfile(file_old_1)):
				col_data_new = np.genfromtxt(
					final_dir_in_work + "/results/output.eng_old_1")
				index = int(col_data_new[0, 0])
				col_data_old = np.genfromtxt(
					final_dir_in_work + "/results/output.eng_old")
				marged_data = np.concatenate(
					(col_data_old[:index - 1], col_data_new), axis=0)

				col_data_new = np.genfromtxt(
					final_dir_in_work + "/results/output.eng")
				index = int(col_data_new[0, 0])
				marged_data = np.concatenate(
					(marged_data[:index - 1], col_data_new), axis=0)
				aa = col_data_new[:, 0]
				final_data_set = marged_data[preskip:(
					int(aa[-1]) - postskip), :]
			elif ((os.path.isfile(file_new)) and (os.path.isfile(file_old_1) == False)):
				print(final_dir_in_work + " -- Restarted data")
				col_data_new = np.genfromtxt(
					final_dir_in_work + "/results/output.eng")
				index = int(col_data_new[0, 0])
				col_data_old = np.genfromtxt(
					final_dir_in_work + "/results/output.eng_old")
				marged_data = np.concatenate(
					(col_data_old[:index - 1], col_data_new), axis=0)
				aa = col_data_new[:, 0]
				final_data_set = marged_data[preskip:(
					int(aa[-1]) - postskip), :]
			elif ((os.path.isfile(file_new) == False) and (os.path.isfile(file_old_1) == False)):
				final_data_set = np.genfromtxt(
					final_dir_in_work +
					"/results/output.eng_old",
					skip_header=preskip,
					skip_footer=postskip)
		else:
			final_data_set = np.genfromtxt(
				final_dir_in_work +
				"/results/output.eng",
				skip_header=preskip,
				skip_footer=postskip)

	if (method == "PIMC"):
		col_block = final_data_set[:, 0]
		col_rot = final_data_set[:, 1]
		col_pot = final_data_set[:, 2]
		col_tot = final_data_set[:, 3]
		col_rotsq = final_data_set[:, 4]
		col_Cv1 = final_data_set[:, 5]
		col_Cv2 = final_data_set[:, 6]

		ncol_block = len(col_block)
		if (int(len(col_block)) != numb_block - (preskip + postskip)):
			print(len(col_block))
			print(final_dir_in_work)

		binary_exponent = int(math.log(len(col_tot)) / math.log(2))
		trunc = int(len(col_tot) - pow(2,binary_exponent))

		col_rot = col_rot[trunc:]
		col_pot = col_pot[trunc:]
		col_tot = col_tot[trunc:]
		col_rotsq = col_rotsq[trunc:]
		col_Cv1 = col_Cv1[trunc:]
		col_Cv2 = col_Cv2[trunc:]

		mean_rot = np.mean(col_rot)
		mean_pot = np.mean(col_pot)
		mean_tot = np.mean(col_tot)
		mean_rotsq = np.mean(col_rotsq)
		mean_Cv1 = np.mean(col_Cv1)
		mean_Cv2 = np.mean(col_Cv2)
		mcbeta = variable_value * numb_bead
		mean_Cv = mcbeta * mcbeta * (mean_Cv1 - mean_Cv2 - mean_tot * mean_tot)

		error_rot = get_error(mean_rot, col_rot, binary_exponent - 6)
		error_pot = get_error(mean_pot, col_pot, binary_exponent - 6)
		error_tot = get_error(mean_tot, col_tot, binary_exponent - 6)
		error_rotsq = get_error(mean_rotsq, col_rotsq, binary_exponent - 6)
		error_Cv1 = get_error(mean_Cv1, col_Cv1, binary_exponent - 6)
		error_Cv2 = get_error(mean_Cv2, col_Cv2, binary_exponent - 6)
		error_Cv = mcbeta * mcbeta * \
			(math.sqrt(error_Cv1 * error_Cv1 + error_Cv2 * error_Cv2 +
					   4.0 * mean_tot * mean_tot * error_tot * error_tot))

		output = '{blocks:^12d}{beads:^10d}{var:^10.6f}{rot:^16.6f}{pot:^16.6f}{tot:^16.6f}{rotsq:^16.6f}{Cv:^16.6f}{er1:^12.6f}{er2:^12.6f}{er3:^12.6f}{er4:^12.6f}{er5:^12.6f}'.format(
			blocks=ncol_block, beads=numb_bead, var=variable_value, rot=mean_rot, pot=mean_pot, tot=mean_tot, rotsq=mean_rotsq, Cv=mean_Cv, er1=error_rot, er2=error_pot, er3=error_tot, er4=error_rotsq, er5=error_Cv)
		output += "\n"

	if (method == "PIGS"):
		col_block = final_data_set[:, 0]
		col_pot = final_data_set[:, 1]
		col_tot = final_data_set[:, 2]

		ncol_block = len(col_block)
		if (int(len(col_block)) != numb_block - (preskip + postskip)):
			print(len(col_block))
			print(final_dir_in_work)

		binary_exponent = int(math.log(len(col_tot)) / math.log(2))
		trunc = int(len(col_tot) - pow(2, binary_exponent))

		col_pot = col_pot[trunc:]
		col_tot = col_tot[trunc:]

		mean_pot = np.mean(col_pot)
		mean_tot = np.mean(col_tot)

		error_pot = get_error(mean_pot, col_pot, binary_exponent - 6)
		error_tot = get_error(mean_tot, col_tot, binary_exponent - 6)

		output = '{blocks:^12d}{beads:^10d}{var:^10.6f}{pot:^16.6f}{tot:^16.6f}{er1:^12.6f}{er2:^12.6f}'.format(
			blocks=ncol_block, beads=numb_bead, var=variable_value, pot=mean_pot, tot=mean_tot, er1=error_pot, er2=error_tot)
		output += "\n"

	return output


def fmtAverageOrderParam(status, variable):
	'''
	This function gives us the output
	'''
	if variable == "rpt_val":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output = "# "
		output += '{blocks:^10}{beads:^10}{var:^10}{eiz:^12}{eiejz:^12}{er1:^12}{er2:^12}'.format(
			blocks='nBlocks',
			beads='nBeads',
			var=variable +
			' invK',
			eiz='<eiz>',
			eiejz='<eiejz>',
			er1='Err-eiz',
			er2='Err-eiejz')
		output += "\n"
		output += '{0:=<80}'.format('#')
		output += "\n"
		return output


def GetAverageOrderParam(
		method,
		numb_molecule,
		numb_bead,
		variable,
		final_dir_in_work,
		preskip,
		postskip,
		numb_block):
	'''
	See PRL 118, 027402 (2017)
	'''
	if (method != 'PIMC'):
		axis_index = {"cost": 0, "phi": 1, "chi": 2}
		axis_read = "cost"
		ndofs = 3
		beads_pos = int((numb_bead - 1) / 2)
		collist = [0]
		for i in range(numb_molecule):
			ncol1 = beads_pos + i * numb_bead
			ncol = axis_index[axis_read] + ncol1 * ndofs
			ncol = ncol + 1
			collist.append(ncol)
			#print(str(ncol)+'th column')
		col = tuple(collist)

		if (os.path.isdir(final_dir_in_work)):
			file_old = final_dir_in_work + "/results/output.xyz_old"
			if (os.path.isfile(file_old)):
				file_new = final_dir_in_work + "/results/output.xyz"
				if (os.path.isfile(file_new)):
					print(final_dir_in_work + " -- Restarted data")
					if "H2O1" in open(file_new).read():
						file_temp = final_dir_in_work + "/results/output_temp.xyz"
						rmstr = int(numb_molecule * numb_bead + 3)
						cmd1 = "tail -n +" + str(rmstr) + \
							" " + file_new + ">" + file_temp
						os.system(cmd1)
						col_data_new = np.genfromtxt(file_temp, usecols=col)
						call(["rm", file_temp])
					else:
						col_data_new = np.genfromtxt(file_new, usecols=col)
					index = int(col_data_new[0, 0])
					col_data_old = np.genfromtxt(file_old, usecols=col)
					marged_data = np.concatenate(
						(col_data_old[:index - 1], col_data_new), axis=0)
					aa = col_data_new[:, 0]
					final_data_set = marged_data[preskip:(
						int(aa[-1]) - postskip), :]
				else:
					final_data_set = np.genfromtxt(
						file_old, usecols=col, skip_header=preskip, skip_footer=postskip)
			else:
				final_data_set = np.genfromtxt(
					final_dir_in_work +
					"/results/output.xyz",
					usecols=col,
					skip_header=preskip,
					skip_footer=postskip)

	ncol_block = len(final_data_set[:, 0])
	if (int(ncol_block) != numb_block - (preskip + postskip)):
		print(ncol_block)
		print(final_dir_in_work)

	binary_exponent = int(math.log(ncol_block) / math.log(2))
	trunc = int(ncol_block - 2**binary_exponent)
	raw_data = np.delete(final_data_set, 0, 1)

	if (numb_molecule == 2):
		raw_data1 = np.absolute(raw_data)
		eiz = np.sum(raw_data1[trunc:, :], axis=1) / numb_molecule
	if (numb_molecule == 11):
		raw_data1 = np.absolute(raw_data)
		eiz = np.sum(raw_data1[trunc:, 2:numb_molecule - 3],
					 axis=1) / len([i for i in range(2, numb_molecule - 2)])
	mean_eiz = np.mean(eiz)
	error_eiz = get_error(mean_eiz, eiz, binary_exponent - 6)

	if (numb_molecule == 2):
		paireiej = [i for i in range(numb_molecule - 1)]
	if (numb_molecule == 11):
		paireiej = [i for i in range(2, numb_molecule - 3)]
	norm_eiejz = len(paireiej)
	eiejz = np.zeros(ncol_block - trunc, dtype=float)
	for i in paireiej:
		eiejz += np.multiply(raw_data[trunc:, i],
							 raw_data[trunc:, i + 1]) / norm_eiejz
	mean_eiejz = np.mean(eiejz)
	error_eiejz = get_error(mean_eiejz, eiejz, binary_exponent - 6)

	output = '{blocks:^12d}{beads:^10d}{var:^10.6f}{eiz:^12.6f}{eiejz:^12.6f}{er1:^12.6f}{er2:^12.6f}'.format(
		blocks=ncol_block, beads=numb_bead, var=variable, eiz=mean_eiz, eiejz=mean_eiejz, er1=error_eiz, er2=error_eiejz)
	output += "\n"
	return output


'''
def GetAverageOrderParam(method,numb_molecule,numb_bead,variable,final_dir_in_work,preskip,postskip,numb_block):
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
		if (method != 'PIMC'):
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
		beads_pos = int((numb_bead-1)/2)
		for i in range(numb_molecule):
				ncol1 = beads_pos+i*numb_bead
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
		if (int(len(col_block)) != numb_block-(preskip+postskip)):
				print(len(col_block))
				print(final_dir_in_work)

		binary_exponent   = int(math.log(len(col_eiz))/math.log(2))
		trunc		 = int(len(col_eiz)-2**binary_exponent)

		col_eiejx = col_eiejx[trunc:]
		col_eiejy = col_eiejy[trunc:]
		col_eiejz = col_eiejz[trunc:]
		col_eiej = col_eiej[trunc:]

		mean_eiejx = np.mean(col_eiejx)
		mean_eiejy = np.mean(col_eiejy)
		mean_eiejz = np.mean(col_eiejz)
		mean_eiej = np.mean(col_eiej)

		error_eiejx = get_error(mean_eiejx, col_eiejx, binary_exponent-6)
		error_eiejy = get_error(mean_eiejy, col_eiejy, binary_exponent-6)
		error_eiejz = get_error(mean_eiejz, col_eiejz, binary_exponent-6)
		error_eiej = get_error(mean_eiej, col_eiej, binary_exponent-6)

		col_eix = col_eix[trunc:]
		col_eiy = col_eiy[trunc:]
		col_eiz = col_eiz[trunc:]

		mean_eix = np.mean(col_eix)
		mean_eiy = np.mean(col_eiy)
		mean_eiz = np.mean(col_eiz)

		error_eix = get_error(mean_eix, col_eix, binary_exponent-6)
		error_eiy = get_error(mean_eiy, col_eiy, binary_exponent-6)
		error_eiz = get_error(mean_eiz, col_eiz, binary_exponent-6)

		output  = '{blocks:^12d}{beads:^10d}{var:^10.6f}{eiejx:^12.6f}{eiejy:^12.6f}{eiejz:^12.6f}{eiej:^12.6f}{eix:^12.6f}{eiy:^12.6f}{eiz:^12.6f}{er1:^12.6f}{er2:^12.6f}{er3:^12.6f}{er4:^12.6f}{er5:^12.6f}{er6:^12.6f}{er7:^12.6f}'.format(blocks=ncol_block,beads=numb_bead, var=variable, eiejx=mean_eiejx, eiejy=mean_eiejy, eiejz=mean_eiejz, eiej=mean_eiej, eix=mean_eix, eiy=mean_eiy, eiz=mean_eiz, er1=error_eiejx, er2=error_eiejy, er3=error_eiejz, er4=error_eiej, er5=error_eix, er6=error_eiy, er7=error_eiz)
		output  += "\n"
		return output
'''


def GetAverageEntropy(
		numb_bead,
		variable,
		final_dir_in_work,
		preskip,
		postskip,
		numb_block,
		ent_method):
	'''
	This function gives us the output
	'''
	if (ent_method == "EXTENDED_ENSMBL"):
		col_block, col_nm, col_dm = genfromtxt(
			final_dir_in_work + "/results/output.rden", unpack=True, usecols=[
				0, 1, 2], skip_header=preskip, skip_footer=postskip)
		ncol_block = len(col_block)
		if (int(len(col_block)) != numb_block - (preskip + postskip)):
			print(len(col_block))
			print(final_dir_in_work)

		binary_exponent = int(math.log(len(col_nm)) / math.log(2))
		trunc = int(len(col_nm) - 2**binary_exponent)

		col_nm = col_nm[trunc:]
		col_dm = col_dm[trunc:]
		mean_nm = np.mean(col_nm)
		mean_dm = np.mean(col_dm)
		purity = mean_nm / mean_dm
		mean_EN = -log(purity)

		error_nm = get_error(mean_nm, col_nm, binary_exponent - 6)
		error_dm = get_error(mean_dm, col_dm, binary_exponent - 6)
		error_purity = abs(purity) * sqrt((error_dm / mean_dm) *
										  (error_dm / mean_dm) + (error_nm / mean_nm) * (error_nm / mean_nm))
		error_EN = sqrt((error_dm / mean_dm) * (error_dm / mean_dm) +
						(error_nm / mean_nm) * (error_nm / mean_nm))

		output = '{blocks:^8d}{beads:^10d}{var:^10.6f}{nm:^12.6f}{dm:^12.6f}{pur:^12.6f}{ent:^12.6f}{er1:^12.6f}{er2:^12.6f}{er3:^12.6f}{er4:^12.6f}'.format(
			blocks=ncol_block, beads=numb_bead, var=variable, nm=mean_nm, dm=mean_dm, pur=purity, ent=mean_EN, er1=error_nm, er2=error_dm, er3=error_purity, er4=error_EN)
		output += "\n"

	if (ent_method == 'BROKENPATH'):
		col_block, col_nm, col_dm = genfromtxt(
			final_dir_in_work + "/results/output.rden", unpack=True, usecols=[
				0, 1, 2], skip_header=preskip, skip_footer=postskip)
		binary_exponent = int(math.log(len(col_nm)) / math.log(2))
		trunc = int(len(col_nm) - 2**binary_exponent)

		col_nm = col_nm[trunc:]
		col_dm = col_dm[trunc:]
		mean_nm = np.mean(col_nm)
		mean_dm = np.mean(col_dm)
		mean_EN = -log(mean_nm / mean_dm)

		error_nm = get_error(mean_nm, col_nm, binary_exponent - 6)
		error_dm = get_error(mean_dm, col_dm, binary_exponent - 6)
		error_EN = sqrt((error_dm / mean_dm) * (error_dm / mean_dm) +
						(error_nm / mean_nm) * (error_nm / mean_nm))

		output = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}'.format(
			numb_bead, variable, mean_nm, mean_dm, mean_EN, error_nm, error_dm, error_EN)
		output += "\n"

	if (ent_method == "SWAP"):
		col_block, col_nm, col_dm, col_TrInv = genfromtxt(
			final_dir_in_work + "/results/output.rden", unpack=True, usecols=[
				0, 1, 2, 3], skip_header=preskip, skip_footer=postskip)
		binary_exponent = int(math.log(len(col_nm)) / math.log(2))
		trunc = int(len(col_nm) - 2**binary_exponent)

		col_nm = col_nm[trunc:]
		col_dm = col_dm[trunc:]
		mean_nm = np.mean(col_nm)
		mean_dm = np.mean(col_dm)
		mean_TrInv = np.mean(col_TrInv)
		purity = 1.0 / mean_TrInv
		mean_EN = -log(purity)

		error_nm = get_error(mean_nm, col_nm, binary_exponent - 6)
		error_dm = get_error(mean_dm, col_dm, binary_exponent - 6)
		error_TrInv = np.std(col_TrInv, ddof=1) / sqrt(len(col_block))
		error_purity = abs(1.0 / (mean_TrInv * mean_TrInv)) * error_TrInv
		# Write the proper equation
		error_EN = abs(1.0 / mean_TrInv) * error_TrInv

		output = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(
			numb_bead, variable, mean_nm, mean_dm, purity, mean_EN, error_nm, error_dm, error_purity, error_EN)
		output += "\n"

	if (ent_method == "UNSWAP"):
		col_block, col_nm, col_dm, col_Tr = genfromtxt(
			final_dir_in_work + "/results/output.rden", unpack=True, usecols=[
				0, 1, 2, 3], skip_header=preskip, skip_footer=postskip)
		binary_exponent = int(math.log(len(col_nm)) / math.log(2))
		trunc = int(len(col_nm) - 2**binary_exponent)

		col_nm = col_nm[trunc:]
		col_dm = col_dm[trunc:]
		mean_nm = np.mean(col_nm)
		mean_dm = np.mean(col_dm)
		mean_purity = np.mean(col_Tr)
		mean_EN = -log(mean_Tr)

		error_nm = get_error(mean_nm, col_nm, binary_exponent - 6)
		error_dm = get_error(mean_dm, col_dm, binary_exponent - 6)
		error_purity = np.std(col_Tr, ddof=1) / sqrt(len(col_block))
		# Write the proper equation
		error_EN = abs(1.0 / mean_TrInv) * error_TrInv

		output = '{0:10d}{1:20.5f}{2:20.5f}{3:20.5f}{4:20.5f}{5:20.5f}{6:20.5f}{7:20.5f}{8:20.5f}{9:20.5f}'.format(
			numb_bead, variable, mean_nm, mean_dm, mean_purity, mean_EN, error_nm, error_dm, error_purity, error_EN)
		output += "\n"

	return output


def fmtAverageEntropy(status, variable, ent_method):
	'''
	This function gives us the output
	'''
	if variable == "rpt_val":
		unit = "(Angstrom)"
	else:
		unit = "(1/K)"

	if status == "analysis":
		output = "#"
		if (ent_method == "EXTENDED_ENSMBL"):
			output += '{blocks:^8}{beads:^10}{var:^10}{nm:^12}{dm:^12}{pur:^12}{ent:^12}{er1:^12}{er2:^12}{er3:^12}{er4:^12}'.format(
				blocks='nBlocks',
				beads='nBeads',
				var=variable +
				' invK',
				nm='<Nm>',
				dm='<Dm>',
				pur='<Purity>',
				ent='<S2>',
				er1='Err-Nm',
				er2='Err-Dm',
				er3='Err-Purity',
				er4='Err-S2')
		if (ent_method == "BROKENPATH"):
			output += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}'.format(
				'Beads',
				variable +
				'  (1/K)',
				'<Nm>',
				'<Dm>',
				'Entropy',
				'Error of Nm',
				'Error of Dm',
				'Error of Entropy')
		if (ent_method == "SWAP"):
			output += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format(
				'Beads',
				variable +
				'  (1/K)',
				'<Nm>',
				'<Dm>',
				'Avg. Purity',
				'Entropy',
				'Error of Nm',
				'Error of Dm',
				'Error of Purity',
				'Error of Entropy')
		if (ent_method == "REGULARPATH"):
			output += '{0:^15}{1:^20}{2:^20}{3:^20}{4:^20}{5:^20}{6:^20}{7:^20}{8:^20}{9:^20}'.format(
				'Beads',
				variable +
				'  (1/K)',
				'<Nm>',
				'<Dm>',
				'Avg. Purity',
				'Entropy',
				'Error of Nm',
				'Error of Dm',
				'Error of Purity',
				'Error of Entropy')
		output += "\n"
		output += '{0:=<123}'.format('#')
		output += "\n"
		return output


def GetAverageEntropyRT(
		particleAList,
		method,
		rotor,
		translational_move,
		rotational_move,
		variable_name,
		rpt_val,
		gfactor,
		dipole_moment,
		parameter_name,
		parameter,
		numb_block,
		numb_pass,
		numbmolecules1,
		molecule,
		ent_method,
		preskip,
		postskip,
		extra_file_name,
		output_dir_path,
		variable,
		crystal,
		final_results_path,
		impurity,
		ext_ent):
	'''
	This function gives us Renyi entropy of a many-rotors system simulated by Ratio trick algorithm.
	'''
	list_nb = mc.Getbeads(method, variable_name)
	ndim_beads = int(len(list_nb))

	purity_combo = np.zeros(ndim_beads, dtype='f')
	err_purity_combo = np.zeros(ndim_beads, dtype='f')
	entropy_combo = np.zeros(ndim_beads, dtype='f')
	err_entropy_combo = np.zeros(ndim_beads, dtype='f')
	col_beads = np.zeros(ndim_beads, dtype=int)
	col_var = np.zeros(ndim_beads, dtype='f')
	ii = 0

	for iBead in list_nb:
		if ((iBead % 2) != 0):
			value = iBead
		else:
			value = iBead + 1

		if (variable_name == "beta"):
			beta = parameter * (value - 1)
			variable = beta
		if (variable_name == "tau"):
			tau = parameter / (value - 1)
			variable = tau

		numb_bead = value

		col_purity = np.zeros(len(particleAList), dtype='f')
		col_err_purity = np.zeros(len(particleAList), dtype='f')

		contition = False
		iPartition = 0
		for partition in particleAList:
			final_file_name = get_working_file(
				method,
				rotor,
				translational_move,
				rotational_move,
				rpt_val,
				gfactor,
				dipole_moment,
				parameter_name,
				parameter,
				numb_block,
				numb_pass,
				numbmolecules1,
				molecule,
				ent_method,
				partition,
				extra_file_name,
				crystal,
				impurity,
				ext_ent)
			execution_bead_dir_name = final_file_name + str(numb_bead)
			final_dir_in_work = output_dir_path + execution_bead_dir_name
			if os.path.isdir(final_dir_in_work):
				condition = True

				file_old = final_dir_in_work + "/results/output.rden_old"
				if os.path.isfile(file_old):
					file_old_1 = final_dir_in_work + "/results/output.rden_old_1"
					file_new = final_dir_in_work + "/results/output.rden"
					if os.path.isfile(file_old_1):
						col_data_new = genfromtxt(
							final_dir_in_work + "/results/output.rden_old_1")
						index = int(col_data_new[0, 0])
						col_data_old = genfromtxt(
							final_dir_in_work + "/results/output.rden_old")
						marged_data = np.concatenate(
							(col_data_old[:index - 1], col_data_new), axis=0)

						col_data_new = genfromtxt(
							final_dir_in_work + "/results/output.rden")
						index = int(col_data_new[0, 0])
						marged_data = np.concatenate(
							(marged_data[:index - 1], col_data_new), axis=0)
						aa = col_data_new[:, 0]
						col_block = marged_data[preskip:(
							int(aa[-1]) - postskip), 0]
						col_nm = marged_data[preskip:(
							int(aa[-1]) - postskip), 1]
						col_dm = marged_data[preskip:(
							int(aa[-1]) - postskip), 2]
					elif ((os.path.isfile(file_new)) and (os.path.isfile(file_old_1) == False)):
						col_data_new = genfromtxt(
							final_dir_in_work + "/results/output.rden")
						index = int(col_data_new[0, 0])
						col_data_old = genfromtxt(
							final_dir_in_work + "/results/output.rden_old")
						marged_data = np.concatenate(
							(col_data_old[:index - 1], col_data_new), axis=0)
						aa = col_data_new[:, 0]
						col_block = marged_data[preskip:(
							int(aa[-1]) - postskip), 0]
						col_nm = marged_data[preskip:(
							int(aa[-1]) - postskip), 1]
						col_dm = marged_data[preskip:(
							int(aa[-1]) - postskip), 2]
					elif ((os.path.isfile(file_new) == False) and (os.path.isfile(file_old_1) == False)):
						col_block, col_nm, col_dm = genfromtxt(
							final_dir_in_work + "/results/output.rden_old", unpack=True, usecols=[
								0, 1, 2], skip_header=preskip, skip_footer=postskip)
				else:
					col_block, col_nm, col_dm = genfromtxt(
						final_dir_in_work + "/results/output.rden", unpack=True, usecols=[
							0, 1, 2], skip_header=preskip, skip_footer=postskip)
				if (int(len(col_block)) != numb_block - (preskip + postskip)):
					print(len(col_block))
					print(execution_bead_dir_name)

				binary_exponent = int(math.log(len(col_nm)) / math.log(2))
				trunc = int(len(col_nm) - 2**binary_exponent)
				mean_nm = np.mean(col_nm[trunc:])
				mean_dm = np.mean(col_dm[trunc:])
				err_nm = get_error(
					mean_nm, col_nm[trunc:], binary_exponent - 6)
				err_dm = get_error(
					mean_dm, col_dm[trunc:], binary_exponent - 6)

				col_purity[iPartition] = mean_nm / mean_dm
				col_err_purity[iPartition] = abs(mean_nm / mean_dm) * sqrt(
					(err_dm / mean_dm) * (err_dm / mean_dm) + (err_nm / mean_nm) * (err_nm / mean_nm))
			else:
				condition = False
			iPartition += 1

		if (condition):
			if (len(particleAList) > 1):
				purity = np.prod(col_purity, axis=0)
				entropy = -log(purity)
				err_purity = abs(
					purity) * sqrt(np.sum(np.square(np.divide(col_err_purity, col_purity))))
				err_entropy = abs(err_purity) / purity
			else:
				purity = col_purity[0]
				entropy = -log(purity)
				err_purity = col_err_purity[0]
				err_entropy = abs(err_purity) / purity

			purity_combo[ii] = purity
			entropy_combo[ii] = entropy
			col_beads[ii] = numb_bead
			col_var[ii] = variable
			err_purity_combo[ii] = err_purity
			err_entropy_combo[ii] = err_entropy
			ii = ii + 1

	#extra_file_name = 'Ratio-Trick-'
	FileAnalysis = GetAnalysisFileName(
		method,
		True,
		rotor,
		translational_move,
		rotational_move,
		rpt_val,
		gfactor,
		dipole_moment,
		parameter_name,
		parameter,
		numb_block,
		numb_pass,
		numbmolecules1,
		molecule,
		ent_method,
		preskip,
		postskip,
		extra_file_name,
		final_results_path,
		partition,
		ext_ent)
	headerString = fmtAverageEntropyRT('tau')
	np.savetxt(FileAnalysis.SaveEntropyRT,
			   np.transpose([col_beads[:ii],
							 col_var[:ii],
							 purity_combo[:ii],
							 entropy_combo[:ii],
							 err_purity_combo[:ii],
							 err_entropy_combo[:ii]]),
			   fmt=['%4d',
					'%10.6f',
					'%10.6f',
					'%10.6f',
					'%10.6f',
					'%10.6f'],
			   header=headerString)
	SavedFile = FileAnalysis.SaveEntropyRT
	FileCheck(method, list_nb, variable_name, SavedFile)
	call(["cat", FileAnalysis.SaveEntropyRT])
	print("Successful execution!")


def fmtAverageEntropyRT(variable):
	'''
	This function gives us the output
	'''
	output = '{0:5}{1:10}{2:10}{3:10}{4:10}{5:30}'.format(
		'P', variable + '(1/K)', 'Purity', 'Entropy', 'ErrorPurity', '  ErrorEntropy')
	return output


def get_input_file(
		method,
		ent_method,
		ent_algorithm,
		temperature,
		numb_bead,
		numb_block,
		numb_pass,
		rotor,
		numb_molecule,
		distance,
		level_bisection,
		step_rot_move,
		step_com_move,
		gfactor,
		dipole_moment,
		particle_a,
		restart_bool,
		numb_block_restart,
		crystal,
		rotor_type,
		translational_move,
		rotational_move,
		propagator_path,
		impurity,
		step_trans1,
		level1):
	'''
	This function modifies parameters in qmc_run.input.
	'''
	input_dir = os.getcwd() + "/"

	if restart_bool:
		replace("job_input", "RESTART",
				"qmc_run.input", "qmc_temp.input")
	else:
		replace("job_input", "START",
				"qmc_run.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("temperature_input", str(temperature),
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("numb_bead_input", str(numb_bead), "qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	if (rotor == "CH3F"):
		replace("pot_read_input", "pesch3fph2-180", "qmc_temp1.input", "qmc_temp.input")
	elif (rotor == "H2O"):
		replace("pot_read", "nopot", "qmc_temp1.input", "qmc_temp.input")
	elif (rotor == "HF"):
		replace("pot_read_input", "nopot", "qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("pot_path_input", input_dir, "qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("method_input", method + "_SIM", "qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("pigs_ent_input", ent_method, "qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	if (ent_algorithm == "WR"):
		replace("pigs_ent_algorithm_input", "RATIOTRICK", "qmc_temp1.input", "qmc_temp.input")
	elif (ent_algorithm == "WOR"):
		replace("pigs_ent_algorithm_input", "NORATIOTRICK", "qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	if translational_move:
		replace("#TRANSLATION", "TRANSLATION", "qmc_temp1.input", "qmc_temp.input")
	else:
		replace("translation_flag_input", "#TRANSLATION", "qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	if rotational_move:
		replace("rotation_flag_input", "ROTATION", "qmc_temp1.input", "qmc_temp.input")
	else:
		replace("rotation_flag_input", "#ROTATION", "qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	if rotational_move:
		replace("molecule_type_input", rotor_type, "qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	if rotational_move:
		replace("propagator_path_input", propagator_path, "qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("numb_molecule_input", str(numb_molecule),
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("numb_rotor_input", str(numb_molecule),
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	if (distance >= 0.0):
		replace("distance_flag_input", "DISTANCE",
				"qmc_temp1.input", "qmc_temp.input")
		call(["mv", "qmc_temp.input", "qmc_temp1.input"])
		replace("distance_input", str(distance),
				"qmc_temp1.input", "qmc_temp.input")
	else:
		replace("distance_flag_input", "",
				"qmc_temp1.input", "qmc_temp.input")
		call(["mv", "qmc_temp.input", "qmc_temp1.input"])
		replace("distance_input", "",
				"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("rotor_type_input", rotor_type,
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("molecule_input", str(rotor),
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("level_bisection_com_input", str(level_bisection),
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("step_rot_input", str(step_rot_move),
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("step_com_input", str(step_com_move),
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	if ((dipole_moment < 0.0) and (gfactor < 0.0)):
		replace("dipole_moment_flag_input", "",
				"qmc_temp1.input", "qmc_temp.input")
		call(["mv", "qmc_temp.input", "qmc_temp1.input"])
		replace("dipole_moment_input", "",
				"qmc_temp1.input", "qmc_temp.input")
	else:
		if (dipole_moment >= 0.0):
			replace("dipole_moment_flag_input", "DIPOLE_MOMENT",
					"qmc_temp1.input", "qmc_temp.input")
			call(["mv", "qmc_temp.input", "qmc_temp1.input"])
			replace("dipole_moment_input", str(
				dipole_moment), "qmc_temp1.input", "qmc_temp.input")

		if (gfactor >= 0.0):
			replace("dipole_moment_flag_input", "DIPOLE_MOMENT",
					"qmc_temp1.input", "qmc_temp.input")
			call(["mv", "qmc_temp.input", "qmc_temp1.input"])
			dipole_moment = GetDipoleMomentFromGFactor(
				rotor, distance, gfactor)
			replace("dipole_moment_input", str(
				dipole_moment), "qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("numb_block_input", str(numb_block),
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("numb_pass_input", str(numb_pass),
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	mcskip = numb_bead * numb_pass
	replace("mskip_input", str(mcskip),
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	mcskip_avg = numb_bead * numb_pass
	replace("mskip_avg_input", str(mcskip_avg),
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	replace("numb_molecule_input", str(particle_a),
			"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	if (crystal):
		replace("read_input", "READMCCOORDS",
				"qmc_temp1.input", "qmc_temp.input")
	else:
		replace("read_input", "",
				"qmc_temp1.input", "qmc_temp.input")
	call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	if (len(impurity) == 4):
		replace("#type_impurity", impurity[0], "qmc_temp1.input", "qmc_temp.input")
		call(["mv", "qmc_temp.input", "qmc_temp1.input"])

		replace("molecule_impurity", impurity[1], "qmc_temp1.input", "qmc_temp.input")
		call(["mv", "qmc_temp.input", "qmc_temp1.input"])

		replace("numb_molecule_impurity", str(
			impurity[2]), "qmc_temp1.input", "qmc_temp.input")
		call(["mv", "qmc_temp.input", "qmc_temp1.input"])

		replace("stat_impurity", impurity[3], "qmc_temp1.input", "qmc_temp.input")
		call(["mv", "qmc_temp.input", "qmc_temp1.input"])

		replace("step_com_impurity", str(step_trans1),
				"qmc_temp1.input", "qmc_temp.input")
		call(["mv", "qmc_temp.input", "qmc_temp1.input"])

		replace("level_bisection_com_impurity", str(level1), "qmc_temp1.input", "qmc_temp.input")
		call(["mv", "qmc_temp.input", "qmc_temp1.input"])

		replace("pot_read_impurity", "nopot", "qmc_temp1.input", "qmc_temp.input")
		call(["mv", "qmc_temp.input", "qmc_temp1.input"])

	call(["mv", "qmc_temp1.input", "qmc.input"])


def get_rotmat(method, molecule, temperature, numb_bead, source_dir_exe):
	'''
	This function generates rotational density matrix - linden.dat
	'''
	print("********************************************************")
	print("          Generation of rotational propagator           ")
	print("********************************************************")
	temperature1 = "%8.6f" % temperature
	if (method == 'PIMC'):
		numb_bead1 = numb_bead
	else:
		numb_bead1 = numb_bead - 1
	command_linden_run = source_dir_exe + "linear_prop/linden.x " + \
		str(temperature) + " " + str(numb_bead1) + " " + \
		str(get_rotational_bconstant(molecule)) + " 15000 -1"
	system(command_linden_run)
	file_rotdens = molecule + "_T" + \
		str(temperature1) + "t" + str(numb_bead1) + ".rot"
	call(["mv", "linden.out", file_rotdens])


def GetTwoBodyDensity(
		rpt_val,
		dipole_moment,
		numb_bead,
		lmax,
		ltotalmax,
		tau,
		molecule):
	'''
	This function generates pair density of two body Hamiltonian
	command_line = "python script_PairDensityGenerator.py -d 1.0 -R 10.05 -P 20 -l-max 2 tau HF 0.2"
	system(command_line)
	'''
	srcCodePath = os.path.expanduser("~") + "/DipoleChain.jl-master/examples/"
	Units = GetUnitConverter()
	BConstant = get_rotational_bconstant(molecule)  # in wavenumber
	BConstantK = BConstant * Units.CMRECIP2KL
	##########################################################################
	RFactorList = GetrAndgFactor(molecule, rpt_val, dipole_moment)
	RFactor = RFactorList[0]
	tau = tau * BConstantK
	if ltotalmax == 0:
		commandRun = "julia " + srcCodePath + "pair_density.jl -R " + \
			str(RFactor) + " --l-max " + str(lmax) + " --tau " + str(tau)
	else:
		commandRun = "julia " + srcCodePath + "pair_density.jl -R " + \
			str(RFactor) + " --l-max " + str(lmax) + \
			" --l-total-max " + str(ltotalmax) + " --tau " + str(tau)
	system(commandRun)


def cagepot(source_dir_exe):
	'''
	This function generates tabulated potential - cagepot.dat
	'''
	command_cagepot_run = source_dir_exe + "tabulated_potential/hfc60.x 100 360"
	system(command_cagepot_run)
	file_cagepot = "hfc60.pot"
	call(["mv", "cagepot.out", file_cagepot])


def jobstring_scratch(
		file_name,
		value,
		thread,
		run_dir,
		molecule,
		temperature,
		numb_bead,
		final_dir,
		dir_run_input_pimc,
		status_cagepot):
	'''
	This function creats jobstring for #PBS script
	'''
	if (thread > 4):
		thread = 4
	job_name = "job_" + str(file_name) + str(value)
	wall_time = "200:00:00"
	processors = "nodes=1:ppn=" + str(thread)
	omp_thread = str(thread)
	output_dir = run_dir + "/results"
	temperature1 = "%5.3f" % temperature
	log_file_path = final_dir + "/"

	input_file = dir_run_input_pimc + "/qmc" + \
		file_name + str(value) + ".input"
	execution_file = dir_run_input_pimc + "/pimc"
	qmc_input = "qmc" + file_name + str(value) + ".input"

	job_string = """#!/bin/bash
#PBS -N %s
#PBS -l wall_time=%s
##PBS -q medium
#PBS -l %s
#PBS -o %s%s.log
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
""" % (job_name, wall_time, processors, log_file_path, job_name, log_file_path, job_name, omp_thread, run_dir, output_dir, input_file, run_dir, file_rotdens, run_dir, run_dir, qmc_input, execution_file, run_dir, run_dir)
	return job_string


def job_submission(
		server_name,
		status,
		translational_move,
		rotational_move,
		dir_run,
		dir_run_job,
		execution_bead_dir_name,
		input_dir,
		execution_file,
		rpt_val,
		ibead,
		numb_bead,
		step_rot_move,
		step_com_move,
		level_bisection,
		temperature,
		numb_block,
		numb_pass,
		rotor,
		numb_molecule,
		gfactor,
		dipole_moment,
		method,
		ent_method,
		ent_algorithm,
		output_dir_path,
		dir_run_input_pimc,
		cpu_run,
		particle_a,
		partition_name,
		status_cagepot,
		user_name,
		output_file_dir,
		source_dir_exe,
		restart_bool,
		numb_block_restart,
		crystal,
		rotor_type,
		spin_isomer,
		impurity,
		step_com_impurity,
		level_impurity):

	arg1 = rpt_val
	step1_trans = step_com_move[ibead]
	level1 = level_bisection[ibead]
	step1_rot = step_rot_move[ibead]
	step_com_impurity1 = step_com_impurity[ibead]
	level_impurity1 = level_impurity[ibead]
	final_dir_in_work = output_dir_path + execution_bead_dir_name

	if (method == 'PIGS'):
		job_name_temp = "pgR" + str(rpt_val) + 'n' + str(numb_molecule) + "b"
	if (method == 'PIMC'):
		job_name_temp = "pm" + str(numb_molecule) + "b"
	if (method == 'ENT'):
		job_name_temp = "et" + str(numb_molecule) + "a" + str(particle_a) + "b"
	job_name = job_name_temp + str(numb_bead) + ".log"

	if not restart_bool:
		os.chdir(output_dir_path)
		if (os.path.isdir(execution_bead_dir_name)):
			print("====")
			print("Error message")
			print("")
			print_message = "Remove " + str(output_dir_path) + str(execution_bead_dir_name)
			print(print_message)
			os.chdir(input_dir)
			return

		os.chdir(input_dir)

		if (rotor_type == "LINEAR"):
			get_rotmat(method, rotor, temperature,
				   numb_bead, source_dir_exe)

			#call(["rm", "-rf", execution_bead_dir_name])

		temperature1 = "%8.6f" % temperature
		if (rotor_type == "LINEAR"):
			if (method == 'PIMC'):
				numb_bead1 = numb_bead
			else:
				numb_bead1 = numb_bead - 1
			file_rotdens = rotor + "_T" + \
				str(temperature1) + "t" + str(numb_bead1) + ".rot"
			call(["mv", file_rotdens, dir_run_input_pimc])
			propagator_path = dir_run_input_pimc + "/"
		else:
			iodevn = spin_isomer
			jmax = 70
			if (method == 'PIMC'):
				numb_bead2 = int(numb_bead)
			else:
				numb_bead2 = int(numb_bead - 1)
			if (rotor == "H2O"):
				if (server_name == "graham"):
					dir_dens = "/scratch/" + user_name + "/rot-dens-asymmetric-top/" + \
						asym.GetDirNameSubmission(
							rotor, temperature1, numb_bead2, iodevn, jmax)
				if (server_name == "nlogn"):
					dir_dens = "/work/" + user_name + "/rot-dens-asymmetric-top/" + \
						asym.GetDirNameSubmission(
							rotor, temperature1, numb_bead2, iodevn, jmax)
			if (rotor == "CH3F"):
				if (server_name == "nlogn"):
					dir_dens = "/work/" + user_name + "/rot-dens-symmetric-top/" + \
						sym.GetDirNameSubmission(
							rotor, temperature1, numb_bead2, 3, jmax)
				if (server_name == "graham"):
					dir_dens = "/scratch/" + user_name + "/rot-dens-symmetric-top/" + \
						sym.GetDirNameSubmission(
							rotor, temperature1, numb_bead2, 3, jmax)
			file_rotdens = "/" + rotor + "_T" + \
				str(dropzeros(temperature1)) + "t" + str(numb_bead2)
			file_rotdens_mod = "/" + rotor + "_T" + \
				str(temperature1) + "t" + str(numb_bead)
			if (os.path.isfile(dir_dens + file_rotdens_mod + ".rho") == False):
				call(["cp", dir_dens + file_rotdens + ".rho",
					  dir_dens + file_rotdens_mod + ".rho"])
				call(["cp", dir_dens + file_rotdens + ".eng",
					  dir_dens + file_rotdens_mod + ".eng"])
				call(["cp", dir_dens + file_rotdens + ".esq",
					  dir_dens + file_rotdens_mod + ".esq"])
			propagator_path = dir_dens + "/"
			call(["cp", "initial_euler_angles_and_com.txt", dir_run_input_pimc])

		if (crystal):
			call(["cp", "LatticeConfig.xyz", dir_run_input_pimc])
	else:
		os.chdir(output_dir_path)

		#
		if not os.path.isdir(execution_bead_dir_name):
			print("")
			print("")
			print("Error message")
			print("")
			print("")
			print_message = output_dir_path + execution_bead_dir_name + "  --- This directory is absent."
			print(print_message)
			os.chdir(input_dir)
			return

		#
		logout_file = dir_run_input_pimc + "/" + job_name

		if "slurmstepd" not in open(logout_file).read():
			if "real" not in open(logout_file).read():
				print_message = output_dir_path + execution_bead_dir_name + "  --- This job is running."
				print(print_message)
				os.chdir(input_dir)
				return

		#
		if ((method == 'PIGS') or (method == "PIMC")):
			file_count = "output.eng"
		if (method == 'ENT'):
			file_count = "output.rden"

		#
		os.chdir(output_dir_path + execution_bead_dir_name + "/results")
		if (os.path.isfile(file_count)):
			col_data_new = genfromtxt(file_count)
			lastIndex = int(col_data_new[-1, 0])
			if ((numb_block - lastIndex) <= 0):
				print(output_dir_path + execution_bead_dir_name)
				print(" Done. Do not need to resubmit it.")
				os.chdir(input_dir)
				return

		#
		os.chdir(input_dir)

		flag = file_operations(method, final_dir_in_work,
							   numb_molecule, numb_bead)
		if (flag == False):
			print(output_dir_path + execution_bead_dir_name)
			print("Already resubmitted.")
			return

		print(output_dir_path + execution_bead_dir_name)
		print(" Just resubmitted.")

		# Rotational density matrices
		temperature1 = "%8.6f" % temperature
		if (rotor_type != "LINEAR"):
			iodevn = spin_isomer
			jmax = 66
			if (method == 'PIMC'):
				numb_bead2 = int(numb_bead)
			else:
				numb_bead2 = int(numb_bead - 1)
		if (rotor == "H2O"):
			if (server_name == "graham"):
				dir_dens = "/scratch/" + user_name + "/rot-dens-asymmetric-top/" + \
					asym.GetDirNameSubmission(
						rotor, temperature1, numb_bead2, iodevn, jmax)
			if (server_name == "nlogn"):
				dir_dens = "/work/" + user_name + "/rot-dens-asymmetric-top/" + \
					asym.GetDirNameSubmission(
						rotor, temperature1, numb_bead2, iodevn, jmax)
		if (rotor == "CH3F"):
			if (server_name == "nlogn"):
				dir_dens = "/work/" + user_name + "/rot-dens-symmetric-top/" + \
					sym.GetDirNameSubmission(
						rotor, temperature1, numb_bead2, 3, jmax)
			if (server_name == "graham"):
				dir_dens = "/scratch/" + user_name + "/rot-dens-symmetric-top/" + \
					sym.GetDirNameSubmission(
						rotor, temperature1, numb_bead2, 3, jmax)
		file_rotdens = "/" + rotor + "_T" + \
			str(dropzeros(temperature1)) + "t" + str(numb_bead2)
		file_rotdens_mod = "/" + rotor + "_T" + \
			str(temperature1) + "t" + str(numb_bead)
		if (os.path.isfile(dir_dens + file_rotdens_mod + ".rho") == False):
			call(["cp", dir_dens + file_rotdens + ".rho",
				  dir_dens + file_rotdens_mod + ".rho"])
			call(["cp", dir_dens + file_rotdens + ".eng",
				  dir_dens + file_rotdens_mod + ".eng"])
			call(["cp", dir_dens + file_rotdens + ".esq",
				  dir_dens + file_rotdens_mod + ".esq"])
		propagator_path = dir_dens + "/"

	# For qmc.imput
	get_input_file(
		method,
		ent_method,
		ent_algorithm,
		temperature,
		numb_bead,
		numb_block,
		numb_pass,
		rotor,
		numb_molecule,
		arg1,
		level1,
		step1_rot,
		step1_trans,
		gfactor,
		dipole_moment,
		particle_a,
		restart_bool,
		numb_block_restart,
		crystal,
		rotor_type,
		translational_move,
		rotational_move,
		propagator_path,
		impurity,
		step_com_impurity,
		level_impurity)

	input_file = "qmc_trotter_number" + str(numb_bead) + ".input"
	call(["mv", "qmc.input", dir_run_input_pimc + "/" + input_file])
	execution_bead_dir_name_path = dir_run_job + execution_bead_dir_name

	# job submission
	fname = 'job-for-P' + str(numb_bead)
	fwrite = open(fname, 'w')

	if ((server_name == "graham") or (server_name == "nlogn")):
		if (cpu_run == "CPU"):
			fwrite.write(
				jobstring_scratch_cpu(
					job_name_temp,
					numb_molecule,
					execution_bead_dir_name_path,
					rotor,
					temperature,
					numb_bead,
					final_dir_in_work,
					dir_run_input_pimc,
					input_dir,
					restart_bool,
					dir_run_job,
					status_cagepot,
					output_dir_path))
		else:
			fwrite.write(
				jobstring_sbatch(
					server_name,
					dir_run,
					job_name_temp,
					numb_molecule,
					execution_bead_dir_name_path,
					rotor,
					temperature,
					numb_bead,
					final_dir_in_work,
					dir_run_input_pimc,
					user_name,
					output_file_dir,
					restart_bool,
					dir_run_job,
					status_cagepot,
					output_dir_path,
					numb_block))
	else:
		fwrite.write(
			jobstring_sbatch(
				server_name,
				dir_run,
				job_name_temp,
				numb_molecule,
				execution_bead_dir_name_path,
				rotor,
				temperature,
				numb_bead,
				final_dir_in_work,
				dir_run_input_pimc,
				user_name,
				output_file_dir,
				restart_bool,
				dir_run_job,
				status_cagepot,
				output_dir_path,
				numb_block))

	fwrite.close()
	call(["mv", fname, dir_run_input_pimc])
	os.chdir(dir_run_input_pimc)

	if (cpu_run == "CPU"):
		call(["chmod", "755", fname])
		#command_pimc_run = "./"+fname + ">"+ final_dir_in_work+"/outpimc"+str(numb_bead)+" & "
		command_pimc_run = "./" + fname + ">outpimc" + str(numb_bead) + " & "
		print(command_pimc_run)
		system(command_pimc_run)
	elif ((server_name == "graham") or (server_name == "nlogn")):
		if (partition_name == user_name):
			call(["sbatch", "-p", user_name, fname])
		else:
			call(["sbatch", fname])
	else:
		call(["chmod", "+x", fname])
		os.system("./" + fname + " &" )
	print("")
	print("***************** Successfully submitted ***************")
	print("")

	os.chdir(input_dir)


def jobstring_scratch_cpu(
		file_name,
		value,
		thread,
		run_dir,
		molecule,
		temperature,
		numb_bead,
		final_dir,
		dir_run_input_pimc,
		input_dir,
		restart_bool,
		dir_run_job,
		status_cagepot,
		output_dir_path):
	'''
	This function creats jobstring for #PBS script
	'''
	omp_thread = str(thread)
	output_dir = run_dir + "/results"
	temperature1 = "%5.3f" % temperature
	file_rotdens = dir_run_input_pimc + "/" + molecule + \
		"_T" + str(temperature1) + "t" + str(numb_bead) + ".rot"

	input_file = dir_run_input_pimc + "/qmc" + \
		file_name + str(value) + ".input"
	execution_file = dir_run_input_pimc + "/pimc"
	qmc_input = "qmc" + file_name + str(value) + ".input"

	job_string = """#!/bin/bash
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
""" % (omp_thread, run_dir, output_dir, input_dir, input_file, run_dir, file_rotdens, run_dir, run_dir, qmc_input, execution_file, run_dir, run_dir)
	return job_string


def jobstring_sbatch(
		server_name,
		dir_run,
		file_name,
		numb_molecule,
		execution_bead_dir_name_path,
		molecule,
		temperature,
		numb_bead,
		final_dir_in_work,
		dir_run_input_pimc,
		user_name,
		output_file_dir,
		restart_bool,
		dir_run_job,
		status_cagepot,
		output_dir_path,
		numb_block):
	'''
	This function creats jobstring for #SBATCH script
	'''
	if (numb_block <= 1000):
		wall_time = "00-00:30"
		thread = 1
	else:
		if (numb_bead >= 160):
			thread = 8
			wall_time = "7-00:00"
		elif ((numb_bead >= 50) and (numb_bead < 160)):
			thread = 1
			wall_time = "03-00:00"
		elif ((numb_bead >= 30) and (numb_bead < 50)):
			thread = 1
			wall_time = "03-00:00"
		else:
			thread = 1
			wall_time = "03-00:00"

	job_name = file_name + str(numb_bead)
	omp_thread = str(thread)
	output_dir = execution_bead_dir_name_path + "/results"
	temperature1 = "%8.6f" % temperature
	file_rotdens = dir_run_input_pimc + "/" + molecule + \
		"_T" + str(temperature1) + "t" + str(numb_bead) + ".*"
	log_file_path = dir_run_input_pimc + "/" + job_name

	qmc_input = "qmc_trotter_number" + str(numb_bead) + ".input"
	input_file = dir_run_input_pimc + "/" + qmc_input
	execution_file = dir_run_input_pimc + "/pimc"
	input_file1 = dir_run_input_pimc + "/initial_euler_angles_and_com.txt"

	if (status_cagepot):
		cagepot_file = dir_run_input_pimc + "/hfc60.pot"
		cagepot_cp = "cp " + cagepot_file + " " + execution_bead_dir_name_path
	else:
		cagepot_cp = ""

	if server_name == "graham":
		mv_cmd = " "
	else:
		if (dir_run == "scratch"):
			mv_cmd = "mv " + execution_bead_dir_name_path + " " + output_dir_path
		if (dir_run == "work"):
			mv_cmd = " "

	if (server_name == "graham"):
		account = "#SBATCH --account=rrg-pnroy"
	else:
		account = ""

	print("")
	print("")
	print("")
	print("**************** Important Note ************************")
	print("")
	print("The full path of the directory where the submitted job is running is given below:")
	print(execution_bead_dir_name_path)
	print("")
	print("The full path of the directory where all the outputs are stored is geven below:")
	print(output_dir)
	print("")
	print("Detailed information about the system and all the Monte Carlo acceptance ratios are saved in the below-mentioned file.")
	print(log_file_path + ".log")
	print("")
	print("The number of threads used = " + str(thread))
	print("")
	print("The outputs exist in the below-mentioned directory after the job is executed successfully.")
	print(final_dir_in_work)
	print("")
	print("********************************************************")
	print("")

	job_string = """#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --output=%s.log
#SBATCH --time=%s
%s
#SBATCH --constraint=broadwell
#SBATCH --mem-per-cpu=1024mb
#SBATCH --cpus-per-task=%s
export OMP_NUM_THREADS=%s
#export OMP_STACKSIZE=1024M
export GOMP_STACKSIZE=1024M
rm -rf %s
mkdir -p %s
mv %s %s
cd %s
cp %s qmc.input
cp %s %s
#cp %s %s
#valgrind --tool=memcheck --leak-check=yes --show-reachable=yes ./pimc
#valgrind --leak-check=full --show-leak-kinds=all --leak-check=yes ./pimc
time ./pimc
%s
""" % (job_name, log_file_path, wall_time, account, omp_thread, omp_thread, execution_bead_dir_name_path, output_dir, input_file, execution_bead_dir_name_path, execution_bead_dir_name_path, qmc_input, execution_file, execution_bead_dir_name_path, input_file1, execution_bead_dir_name_path, mv_cmd)

	job_string_desktop = """#!/bin/bash
export OMP_NUM_THREADS=%s
export GOMP_STACKSIZE=1024M
rm -rf %s
mkdir -p %s
mv %s %s
cd %s
cp %s qmc.input
cp %s %s
./pimc 1>> %s.log 2>&1
""" % (omp_thread, execution_bead_dir_name_path, output_dir, input_file, execution_bead_dir_name_path, execution_bead_dir_name_path, qmc_input, execution_file, execution_bead_dir_name_path, log_file_path)


	job_string_restart = """#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --output=%s.log
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
""" % (job_name, log_file_path, wall_time, account, omp_thread, omp_thread, final_dir_in_work, dir_run_job, input_file, execution_bead_dir_name_path, execution_bead_dir_name_path, qmc_input, execution_file, execution_bead_dir_name_path, mv_cmd)

	if restart_bool:
		return job_string_restart
	else:
		if (server_name == "graham"):
			return job_string
		else:
			return job_string_desktop


def GetRotEnergy(molecule, jrot):
	Energy = get_rotational_bconstant(molecule) * jrot * (jrot + 1.0)
	return Energy


def GetAvgRotEnergy(molecule, beta):
	CMRECIP2KL = 1.4387672
	Zsum = 0.0
	Nsum = 0.0
	for jrot in range(0, 10000, 1):
		BoltzmannProb = exp(-beta * GetRotEnergy(molecule, jrot) * CMRECIP2KL)
		if (BoltzmannProb > 10e-16):
			Zsum += (2 * jrot + 1.0) * BoltzmannProb
			Nsum += (2 * jrot + 1.0) * \
				GetRotEnergy(molecule, jrot) * BoltzmannProb
		else:
			break
	AvgEnergy = Nsum * CMRECIP2KL / Zsum
	return AvgEnergy

def get_working_file(
		method,
		molecular_system,
		numb_molecule,
		translational_move,
		rotational_move,
		parameter_name,
		parameter_value,
		rpt_val,
		dipole_moment,
		gfactor,
		numb_block,
		numb_pass,
		rotor_name,
		extra_file_name,
		crystal,
		impurity,
		ent_method,
		particle_a,
		ent_algorithm):

	#add					 = "-NumTimes"
	add = ""
	if (method == "ENT"):
		add1 = "-ParticleA" + str(particle_a)
		add2 = "-"
	else:
		if (len(impurity) == 4):
			add1 = "-" + str(impurity[2]) + "-" + impurity[1]
		else:
			add1 = ""
		add2 = ""

	name_layer1 = parameter_name + str(parameter_value) + "Kinv-Blocks" + str(numb_block) + "-Passes" + str(
		numb_pass) + add + "-System-" + str(numb_molecule) + str(rotor_name) + add1 + add2

	front_layer = method+"-"+extra_file_name

	if ((translational_move) and (rotational_move)):
		front_layer += "TransAndRotDOFs-"
	if ((translational_move == False) and (rotational_move)):
		front_layer += "RotDOFs-"
	if ((translational_move) and (rotational_move == False)):
		front_layer += "TransDOFs-"

	if (rpt_val >= 0.0):
		name_rpt = "Rpt" + str(rpt_val) + "Angstrom-"
	else:
		name_rpt = ""

	if (dipole_moment >= 0.0):
		if (numb_molecule > 1):
			name_dipole_moment = "Dipole-Moment" + \
				str(dipole_moment) + "Debye-"
		else:
			name_dipole_moment = "Field" + str(dipole_moment) + "Kinv-"
	else:
		name_dipole_moment = ""

	if (gfactor >= 0.0):
		name_gfactor = "gFactor" + str(gfactor) + "-"
	else:
		name_gfactor = ""

	final_file_name = front_layer + name_rpt + \
		name_dipole_moment + name_gfactor + name_layer1
	if (method == "ENT"):
		final_file_name += ent_method + "-" + ent_algorithm + "-"

	return final_file_name


def get_dmrg_working_file(
		method,
		rotor_name,
		numb_rotor,
		rpt_value,
		dipole_moment,
		l_max,
		l_total_max):

	system_fragment = "-of-" + str(numb_rotor) + str(rotor_name) 

	parameter_fragment = ""
	if (rpt_value >= 0.0):
		parameter_fragment += "-Rpt" + str(rpt_value) + "Angstrom"

	if (dipole_moment >= 0.0):
		parameter_fragment += "-Dipole-Moment" + str(dipole_moment) + "Debye"

	parameter_fragment += "-lmax" + str(l_max) + "-ltotalmax" + str(l_total_max)
	first_fragment = method + "-output" 
	final_file_name = first_fragment + system_fragment + parameter_fragment

	return final_file_name


class GetAnalysisFileName:
	def __init__(
			self,
			input_dir1,
			method1,
			method2,
			molecular_system1,
			molecule_rot1,
			numb_molecule1,
			translational_move1,
			rotational_move1,
			parameter_name1,
			parameter_value1,
			rpt_val1,
			dipole_moment1,
			gfacor1,
			numb_block1,
			numb_pass1,
			preskip1,
			postskip1,
			extra_file_name1,
			particle_a1,
			ent_method1,
			ent_algorithm):
		self.method = method1
		self.rotor = molecule_rot1
		self.translational_move = translational_move1
		self.rotational_move = rotational_move1
		self.rpt_val = rpt_val1
		self.dipole_moment = dipole_moment1
		self.gfactor = gfacor1
		self.parameter_value = parameter_value1
		self.parameter_name = parameter_name1
		self.numb_block = numb_block1
		self.numb_pass = numb_pass1
		self.numb_molecule = numb_molecule1
		self.molecular_system = molecular_system1
		self.ent_method = ent_method1
		self.preskip = preskip1
		self.postskip = postskip1
		self.extra_file_name = extra_file_name1
		self.input_dir = input_dir1
		self.particle_a = particle_a1

		if (self.method == "ENT"):
			front_layer = "ENT-" + self.extra_file_name
			add1 = "-ParticleA" + str(self.particle_a)
			add2 = "-" + self.ent_method + "-" + ent_algorithm
		else:
			front_layer = self.method + "-" + self.extra_file_name
			add1 = ""
			add2 = ""

		if ((self.translational_move) and (self.rotational_move)):
			front_layer += "TransAndRotDOFs-"
		if ((self.translational_move == False) and (self.rotational_move)):
			front_layer += "RotDOFs-"
		if ((self.translational_move) and (self.rotational_move == False)):
			front_layer += "TransDOFs-"

		if (self.rpt_val >= 0.0):
			name_rpt = "Rpt" + str(self.rpt_val) + "Angstrom-"
		else:
			name_rpt = ""

		if (self.dipole_moment >= 0.0):
			if (self.numb_molecule > 1):
				name_dipole_moment = "Dipole-Moment" + \
					str(self.dipole_moment) + "Debye-"
			else:
				name_dipole_moment = "Field" + \
					str(self.dipole_moment) + "Kinv-"
		else:
			name_dipole_moment = ""

		if (self.gfactor >= 0.0):
			name_gfactor = "gFactor" + str(self.gfactor) + "-"
		else:
			name_gfactor = ""

		if (self.parameter_name == "beta"):
			variable_name = "tau"
		if (self.parameter_name == "tau"):
			variable_name = "beta"

		name_layer1 = "vs-" + str(variable_name) + "-fixed-" + self.parameter_name + \
			str(self.parameter_value) + "Kinv-Blocks" + str(self.numb_block)
		name_layer1 += "-Passes" + str(self.numb_pass) + "-System" + str(self.numb_molecule) + str(
			self.molecular_system) + add1 + add2 + "-preskip" + str(self.preskip) + "-postskip" + str(self.postskip)

		file_output1 = front_layer + name_rpt + \
			name_dipole_moment + name_gfactor + "Energy-"
		file_output2 = front_layer + name_rpt + \
			name_dipole_moment + name_gfactor + "correlation-"
		file_output3 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "total-correlation-function-"
		file_output4 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "X-component-correlation-function-"
		file_output5 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "Y-component-correlation-function-"
		file_output6 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "Z-component-correlation-function-"
		file_output7 = front_layer + name_rpt + name_dipole_moment + \
			name_gfactor + "XandY-component-correlation-function-"
		file_output8 = front_layer + name_rpt + \
			name_dipole_moment + name_gfactor + "Entropy-"

		self.save_file_energy = self.input_dir + file_output1 + name_layer1 + ".txt"
		self.save_file_correlation = self.input_dir + file_output2 + name_layer1 + ".txt"
		self.SaveEntropy = self.input_dir + file_output8 + name_layer1 + ".txt"

		if (method2 == False):
			if os.path.exists(self.SaveEntropy):
				os.remove(self.SaveEntropy)
			if os.path.exists(self.save_file_energy):
				os.remove(self.save_file_energy)
			if os.path.exists(self.save_file_correlation):
				os.remove(self.save_file_correlation)

		if (self.method != "ENT"):
			print("")
			print("")
			print("****************** Important message *******************")
			print("")
			print(f'Name of the file where the potential and the total energies computed for several {variable_name} values are stored are given below:')
			print(file_output1 + name_layer1 + ".txt")
			print("")
			print("Final analyzed results are stored in - ")
			print(self.input_dir)
			print("")
			print("********************************************************")

		if (self.method == "ENT"):
			name_layer1RT = "vs-" + str(self.variable_name) + "-fixed-" + self.parameter_name + str(
				self.parameter) + "Kinv-Blocks" + str(self.numb_block)
			name_layer1RT += "-Passes" + str(self.numb_pass) + "-System" + str(self.numb_molecule) + str(
				self.molecule) + add1 + add2 + "-preskip" + str(self.preskip) + "-postskip" + str(self.postskip)

			self.SaveEntropyRT = self.input_dir + file_output8 + name_layer1RT + ".txt"
			if (ent_algorithm != "WR"):
				# print(self.SaveEntropyRT)
				print(self.input_dir)
				print("")
				print("Final results - Entropy vs " + str(self.variable_name))
				print(self.SaveEntropy)
				print("********************************************************")


def FileCheck(method, list_nb, variable_name, SavedFile):
	for i in list_nb:
		if (method == "PIMC"):
			if ((i % 2) == 0):
				bead = i
			else:
				bead = i + 1
		else:
			if ((i % 2) != 0):
				bead = i
			else:
				bead = i + 1

		string = str(bead)
		if string in open(SavedFile).read():
			return

	if string not in open(SavedFile).read():
		call(["rm", SavedFile])


class UnitConverter:
	def __init__(self):
		self.bohr_radius = 0.5291772108			# angstrom
		self.hartree_to_jule = 4.359748e-18		
		self.hartree_to_kelvin = 3.157732e+05	
		self.cm_recip_to_kelvin = 1.4387672	

		self.au_to_debye = 1.0 / 0.39343
		self.au_to_cm_inverse = 219474.63137
		self.au_to_kelvin = 315777.0
		self.kcal_per_mol_to_cm_inverse = 349.75509
		self.kcal_per_mole_to_kelvin = 503.228


def get_gfactor_rfactor(rotor, rcom, dipole_moment):
	'''
	It calculates g and R value
	'''
	unit_conversion_factor = UnitConverter()
	rotational_b_constant = get_rotational_bconstant(rotor)  # in wavenumber
	dipole_moment_in_au = dipole_moment / unit_conversion_factor.au_to_debye
	rcom_in_au = rcom / unit_conversion_factor.bohr_radius
	rotational_b_constant_in_au = rotational_b_constant / unit_conversion_factor.au_to_cm_inverse
	rfactor = rcom_in_au / \
		((dipole_moment_in_au * dipole_moment_in_au / rotational_b_constant_in_au)**(1.0 / 3.0))
	gfactor = (dipole_moment_in_au * dipole_moment_in_au) / \
		(rcom_in_au * rcom_in_au * rcom_in_au * rotational_b_constant_in_au)
	print_message = " Dipole Moment = " + \
		str(dipole_moment) + " gFactor = " + \
		str(gfactor) + " Dmitri's rFactor = " + str(rfactor)
	print(print_message)
	return_list = {"g":gfactor,"R":rfactor}
	return return_list


def GetDipoleMomentFromGFactor(molecule, RCOM, gFactor):
	'''
	It extracts dipole moment from a g value - g
	'''
	Units = GetUnitConverter()
	BConstant = get_rotational_bconstant(molecule)  # in wavenumber
	RCOMAU = RCOM / Units.BOHRRADIUS
	BConstantAU = BConstant / Units.AuToCmInverse
	DipoleMomentAU = sqrt(gFactor * RCOMAU * RCOMAU * RCOMAU * BConstantAU)
	dipole_moment = DipoleMomentAU * Units.AuToDebye
	return dipole_moment


def get_dmrg_result(
	server_name,
	root_dir_execution,
	method,
	rotor_name,
	dir_run_job,
	job_execution_dir,
	input_dir,
	source_code_dir,
	numb_rotor,
	rfactor,
	l_max,
	l_total_max):
	'''
	It will give ground state energy, von Neuman and Renyi entropies computed by diagonalizing full Hamiltonian matrix. It is developed by Dmitri https://github.com/0/DipoleChain.jl
	'''
	log_file = method + ".log"

	os.chdir(dir_run_job)
	if (os.path.isdir(job_execution_dir)):
		print("====")
		print("Error message")
		print("")
		print_message = "Remove " + str(dir_run_job) + str(job_execution_dir)
		print(print_message)
		os.chdir(input_dir)
		return

	os.chdir(input_dir)

	job_execution_dir_path = dir_run_job + job_execution_dir
	print(job_execution_dir_path)
	temp_dir = os.path.dirname(job_execution_dir_path)
	if not os.path.exists(temp_dir):
		os.makedirs(temp_dir)

	job_submission_file = job_execution_dir_path + "/job-submission-script.sh"

	fwrite = open(submission_file_name, 'w')


	fwrite.close()
	exit()
	os.chdir(job_execution_dir_path)
	if ((server_name == "graham") or (server_name == "nlogn")):
		call(["sbatch", submission_file_name])
	else:
		call(["chmod", "+x", submission_file_name])
		os.system("./" + submission_file_name + " &" )
	print("")
	print("***************** Successfully submitted ***************")
	print("")

	os.chdir(input_dir)


def get_job_submission_script(
		server_name,
		root_dir_execution,
		file_name,
		numb_rotor,
		job_execution_dir_path,
		molecule,
		final_dir_in_work,
		dir_run_input_pimc,
		user_name,
		output_file_dir,
		dir_run_job,
		final_output_dir):

	if (numb_rotor <= 5):
		wall_time = "00-03:00"
	else:
		wall_time = "01-00:00"

	job_name = file_name + str(numb_bead)
	omp_thread = str(thread)
	log_file_path = dir_run_input_pimc + "/" + job_name

	if (server_name == "graham"):
		mv_cmd = " "
	else:
		if (dir_run == "scratch"):
			mv_cmd = "mv " + job_execution_dir_path + " " + final_output_dir
		if (dir_run == "work"):
			mv_cmd = " "

	if (server_name == "graham"):
		account = "#SBATCH --account=rrg-pnroy"
	else:
		account = ""

	cmd_run = "julia " + source_code_dir + "diagonalization.jl -R " + str(rfactor) + " -N " + str(
		numb_molecule) + " --l-max " + str(l_max) + " --l-total-max " + str(l_total_max) + " --A-start 1 --A-size 1"

	print("")
	print("")
	print("")
	print("**************** Important Note ************************")
	print("")
	print("The full path of the directory where the submitted job is running is given below:")
	print(job_execution_dir_path)
	print("")
	print("Detailed information about the system and all the Monte Carlo acceptance ratios are saved in the below-mentioned file.")
	print(log_file_path + ".log")
	print("")
	print("The outputs exist in the below-mentioned directory after the job is executed successfully.")
	print(final_output_dir)
	print("")
	print("********************************************************")
	print("")

	job_string = """#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --output=%s.log
#SBATCH --time=%s
%s
#SBATCH --constraint=broadwell
#SBATCH --mem-per-cpu=1024mb
#SBATCH --cpus-per-task=%s
export OMP_NUM_THREADS=%s
#export OMP_STACKSIZE=1024M
export GOMP_STACKSIZE=1024M
rm -rf %s
mkdir -p %s
cd %s
time ./pimc
%s
""" % (job_name, log_file_path, wall_time, account, omp_thread, omp_thread, execution_bead_dir_name_path, output_dir, input_file, execution_bead_dir_name_path, execution_bead_dir_name_path, qmc_input, execution_file, execution_bead_dir_name_path, input_file1, execution_bead_dir_name_path, mv_cmd)

	job_string_desktop = """#!/bin/bash
export OMP_NUM_THREADS=%s
export GOMP_STACKSIZE=1024M
rm -rf %s
mkdir -p %s
cd %s
./pimc 1>> %s.log 2>&1
""" % (omp_thread, execution_bead_dir_name_path, output_dir, input_file, execution_bead_dir_name_path, execution_bead_dir_name_path, qmc_input, execution_file, execution_bead_dir_name_path, log_file_path)

def GetPairDensity(
		FilePlotName,
		srcCodePath,
		RFactor,
		numb_molecule,
		loop,
		particle_a,
		rotor,
		rpt_val,
		dipole_moment,
		parameter,
		BConstantK,
		variable_name,
		method):
	FileToBeSavedDensity = FilePlotName + ".txt"
	for numb_bead in loop:
		print(numb_bead)
		parameterR = parameter * BConstantK
		commandRun = "julia " + srcCodePath + "pair_density.jl -R " + \
			str(RFactor) + " --l-max 2 --l-total-max 2 --tau " + str(parameterR)
		print(commandRun)
		call(["rm", "outputDensity.txt"])
		system(commandRun)
		call(["mv", "outputDensity.txt", FileToBeSavedDensity])

# def GetEntropyRT(status, maxloop, method, rotor, translational_move, rotational_move, variable_name, rpt_val, gfactor, dipole_moment, parameter_name, parameter, numb_block, numb_pass, numbmolecules1, molecule, ent_method, preskip, postskip, extra_file_name, final_results_path, particle_a, variable):
	#FileAnalysis = GetAnalysisFileName(method, True, rotor, translational_move, rotational_move, variable_name, rpt_val, gfactor, dipole_moment, parameter_name, parameter, numb_block, numb_pass, numbmolecules1, molecule, ent_method, preskip, postskip, extra_file_name, final_results_path, particle_a)
#	   GetAverageEntropyRT(maxloop, method, rotor, translational_move, rotational_move, variable_name, rpt_val, gfactor, dipole_moment, parameter_name, parameter, numb_block, numb_pass, numbmolecules1, molecule, ent_method, preskip, postskip, extra_file_name, final_results_path, particle_a, variable)


def GetPreFactDDPot(molecule, RCOM, dipole_moment):
	'''
	It calculates the pre factor in inverse temperature of dipole - dipole interaction potential
	'''
	Units = GetUnitConverter()
	DipoleMomentAU = dipole_moment / Units.AuToDebye
	RCOMAU = RCOM / Units.BOHRRADIUS
	preFact = (DipoleMomentAU * DipoleMomentAU) / (RCOMAU * RCOMAU * RCOMAU)
	preFact = preFact * Units.HARTREE2KL
	printingmessage = " dipole_moment = " + \
		str(dipole_moment) + " Debye and the prefactor of the dipole-dipole interaction potential = " + \
		str(preFact) + " K^-1"
	print(printingmessage)


def GetRenamingFunc(
		dir_run_input_pimc,
		dir_input_pimc_renamed,
		output_dir_path,
		execution_bead_dir_name,
		folder_renamed,
		input_dir):
	#final_dir_in_work = output_dir_path + execution_bead_dir_name

	call(["mkdir", "-p", dir_input_pimc_renamed])
	cmd_run = "cp -r " + dir_run_input_pimc + "/*  " + dir_input_pimc_renamed + "/"
	os.system(cmd_run)
	os.chdir(output_dir_path)
	if (os.path.isdir(execution_bead_dir_name)):
		call(["mv", execution_bead_dir_name, folder_renamed])
		print(dir_input_pimc_renamed)
		print_message = "move " + str(output_dir_path) + str(execution_bead_dir_name)
		print(print_message)
	os.chdir(input_dir)


def RemoveFiles(
		method,
		numb_bead,
		temperature,
		rotor,
		rotor_type,
		preskip,
		postskip,
		numb_block,
		final_dir_in_work):

	col_block = genfromtxt(
		final_dir_in_work +
		"/results/output.eng",
		unpack=True,
		usecols=[0],
		skip_header=preskip,
		skip_footer=postskip)
	if (int(len(col_block)) == numb_block - (preskip + postskip)):

		temperature1 = "%8.6f" % temperature
		if (rotor_type == "LINEAR"):
			file_rotdens = rotor + "_T" + \
				str(temperature1) + "t" + str(numb_bead) + ".rot"
			call(["rm", final_dir_in_work + "/" + file_rotdens])
		else:
			if (method == 'PIMC'):
				numb_bead2 = int(numb_bead)
			else:
				numb_bead2 = int(numb_bead - 1)
			file_rotdens_mod = rotor + "_T" + \
				str(temperature1) + "t" + str(numb_bead)
			if (os.path.exists(final_dir_in_work + "/" + file_rotdens_mod + ".rho")):
				call(["rm", final_dir_in_work + "/" + file_rotdens_mod + ".rho"])
				call(["rm", final_dir_in_work + "/" + file_rotdens_mod + ".eng"])
				call(["rm", final_dir_in_work + "/" + file_rotdens_mod + ".esq"])
