'''
Some functions to mess with the files.
'''
import numpy as np
import pandas as pd
from protein import Protein

def csv_writer(filename, array):
	'''
	Args:
		filename: Filename of the .csv being saved.
		array: Array of values to be saved.
	'''
	path = 'general_files/'
	df = pd.DataFrame(np.array(array))
	df.to_csv(path+filename, sep=',', header=False, index=False)

	with open(path+filename) as f: #This part is to remove extra newline in the end
		lines = f.readlines()
		last = len(lines) - 1
		lines[last] = lines[last].replace('\r','').replace('\n','')
	with open(filename, 'w') as wr:
		wr.writelines(lines)

def csv_loader(filename, obj=False):
	'''
	Args:
		filename: Name of the file being loaded.
		obj: Optional value indicating is the .csv being loaded to 
			protein objects.
	Returns:
		array: Array of protein objects or values.
	Raises:
		
	'''
	path = 'general_files/'
	df = pd.read_csv(path+filename)
	array = []
	for val in df.values:
		if obj:
			protein = Protein(val[0], val[1])
			if len(protein.get_splitted_ec_number()) > 0: #this is to take care i.e. 'n2' -numbers, they will have the 'n' removed
				array.append(protein) #4942
			#else: #44
		else:
			array.append([val[0], val[1]])
	return array

def npy_saver(protein_array):
	'''
	Saves the curvature and torsion of each protein to a single file, and 
	different proteins to separate files. Curvature is index 0 and 
	torsion index 1.
	Numpy array is used because of simplicity and efficiency.

	Args:
		protein_array: contains protein objects.
	'''
	path = 'proteins_curvature_torsion/'
	for protein in protein_array:
		if len(protein.get_curvature()) > 0 and len(protein.get_torsion()) > 0:
			arr = np.array([protein.get_curvature(), protein.get_torsion()])
			print(arr[0])
			file = str(protein.get_uniprot_id()) + '.npy'
			np.save(path+file, arr)

def npy_loader(connection_array):
	'''
	Args:
		connection_array: All protein objects found in the uniprot_sport.dat
	Returns:
		protein_array: Array of protein objects, that have curvature and torsion.
	Raises:
		FileNotFoundError: An error occured when trying to read a file that does
			not exists. Error is expected to happen often and simply passed.
	'''
	path = 'proteins_curvature_torsion/'
	protein_array = []
	for protein in connection_array:
		file = str(protein.get_uniprot_id()) + '.npy'
		try:
			arr = np.load(path+file)
			protein.set_curvature(arr[0])
			protein.set_torsion(arr[1])
			protein_array.append(protein)
		except FileNotFoundError:
			pass
	return protein_array