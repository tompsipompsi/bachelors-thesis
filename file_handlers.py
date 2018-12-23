'''
Some functions to mess with the files.
'''
import numpy as np
import pandas as pd
from protein import Protein

def csv_writer(path, filename, array):
	'''
	Args:
		filename: Filename of the .csv being saved.
		array: Array of values to be saved.
	'''
	df = pd.DataFrame(np.array(array))
	df.to_csv(path+filename, sep=',', header=False, index=False)

	with open(path+filename) as f: #This part is to remove extra newline in the end
		lines = f.readlines()
		last = len(lines) - 1
		lines[last] = lines[last].replace('\r','').replace('\n','')
	with open(filename, 'w') as wr:
		wr.writelines(lines)

def csv_loader(path, filename):
	'''
	Args:
		filename: Name of the file being loaded.
		obj: Optional value indicating is the .csv being loaded to 
			protein objects.
	Returns:
		array: Array of protein objects or values.
	Raises:
		
	'''
	df = pd.read_csv(path+filename)
	array = []
	for val in df.values:
		protein = Protein(val[0], val[1])
		if len(protein.get_splitted_ec_number()) > 0: #this is to take care i.e. 'n2' -numbers, they will have the 'n' removed
			array.append(protein)
	return array

def npy_saver(path, protein_array, save_ca=True):
	'''
	Saves the curvature and torsion of each protein to a single file, and 
	different proteins to separate files. Curvature is index 0 and 
	torsion index 1.
	Alternatively saves the alpha carbon coordinates if the save_ca -parameter
	is set to true.
	Numpy array is used because of simplicity and efficiency.

	Args:
		protein_array: contains protein objects.
	'''
	for protein in protein_array:
		if save_ca:
			arr = np.array([protein.get_ca_coordinates()])
		else:
			if len(protein.get_curvature()) > 0 and len(protein.get_torsion()) > 0:
				arr = np.array([protein.get_curvature(), protein.get_torsion()])
			else:
				continue
		print(arr[0])
		file = str(protein.get_uniprot_id()) + '.npy'
		np.save(path+file, arr)

def npy_loader(path, connection_array, is_ca=True):
	'''
	Args:
		connection_array: All protein objects found in the uniprot_sport.dat
	Returns:
		protein_array: Array of protein objects, that have curvature and torsion.
	Raises:
		FileNotFoundError: An error occured when trying to read a file that does
			not exists. Error is expected to happen often and simply passed.
	'''
	protein_array = []
	for protein in connection_array:
		file = str(protein.get_uniprot_id()) + '.npy'
		try:
			arr = np.load(path+file)
			if is_ca:
				protein.set_ca_coordinates(arr[0])
			else:
				protein.set_curvature(arr[0])
				protein.set_torsion(arr[1])
			protein_array.append(protein)
		except FileNotFoundError:
			pass
	return protein_array

