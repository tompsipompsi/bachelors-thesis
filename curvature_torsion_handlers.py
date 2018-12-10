'''
Functions that play with curvature, torsion and pdb files.
'''
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from backbone_curve_cls import backbone_curve_cls
#import pylab

def read_id_connection(cleaned_ec):
	'''
	Args:
		cleaned_ec: EC number with missing values removed.
	Returns:
		PDB id in upper case if found, else None.
	Raises:
		
	'''
	path = 'proteins_ec/'
	file = 'download_pdbs_' + cleaned_ec[0] +'.csv'
	pdb_data_arr = []
	with open(path+file) as data:
		for line in data:
			if cleaned_ec in line: #excludes some of the wrong lines, but not all
				line = line.rstrip('\n')
				one_protein_arr = line.split('\t')
				pdb_data_arr.append(one_protein_arr)

	for line in pdb_data_arr: #find first pdb_id that has human 
		if line[0] == cleaned_ec and line[2] == 'Homo sapiens':
			return line[1].upper()
	for line in pdb_data_arr: #if no human pdbs found
		if line[0] == cleaned_ec:
			return line[1].upper()
	return None

def read_pdb_file(pdb_id):
	'''
	Args:
		pdb_id: String representing a PDB id.
	Returns:
		Numpy array of alpha carbon (CA) coordinates of one protein.
	Raises:
		FileNotFoundError: If no file for the PDB id is found. 
			In this case returns empty array.
	'''
	path = 'proteins_ec/pdb_files_' + pdb_id[0] + '/'
	file = 'pdb_' + pdb_id + '.pdb'

	parser = PDBParser(PERMISSIVE=1)
	ca_coordinates = []
	try:
		structure = parser.get_structure(pdb_id, path+file)
	except FileNotFoundError:
		return ca_coordinates

	for model in structure.get_list():
		for chain in model.get_list():
			for residue in chain.get_list():
				if residue.has_id("CA"):
					ca = residue["CA"]
					#print(ca.get_coord())
					ca_coordinates.append(np.array(ca.get_coord()))
	return np.array(ca_coordinates)

def interpolate_curve(ca_coordinates):
	'''
	Args:
		ca_coordinates: Numpy array of alpha carbon (CA) 
			coordinates of one protein.
	Returns:
		interpolated: Numpy array of CA coordinates with 
			missing values interpolated.
	'''
	m, n = ca_coordinates.shape
	interpolated = np.linspace(start=0, stop=1, num=m)
	return interpolated

def curvature_and_torsion(ca_coordinates):
	'''
	Args:
		ca_coordinates: Numpy array of alpha carbon (CA) 
			coordinates of one protein.
	Returns:
		curvature: Numpy array of curvature of one protein.
		torsion: Numpy array of torsion of one protein.
	Raises:
		ValueError: If alpha carbon (CA) coordinates have any
			problems. If happens, returns empty arrays.
	'''
	bcc = backbone_curve_cls()
	try:
		tck, u = bcc.splprep(ca_coordinates)
	except ValueError:
		return [], []

	tx = interpolate_curve(ca_coordinates)
	lder = bcc.splade(tx, tck)

	lder = np.array(lder) #convert the list to an array for curvature and torsion functions
	curvature = bcc.curvature(lder)
	torsion = bcc.torsion(lder)
	return curvature, torsion

def set_curvature_and_torsion(protein_array):
	'''
	Args:
		protein_array: Array of protein objects.
	Returns:
		curvature_torsion_arr: Array of protein objects that
			have curvature an torsion set.
	'''
	curvature_torsion_arr = []
	for protein in protein_array:
		cleaned_ec = protein.get_cleaned_ec_number()
		pdb_id = read_id_connection(cleaned_ec)
		if pdb_id != None:
			ca_coordinates = read_pdb_file(pdb_id)
			if len(ca_coordinates) > 0:
				curvature, torsion = curvature_and_torsion(ca_coordinates)
				protein.set_pdb_id(pdb_id) 
				protein.set_curvature(curvature)
				protein.set_torsion(torsion)
				curvature_torsion_arr.append(protein)
			else:
				print(cleaned_ec)
				#break #remove to read all
		else:
			print(pdb_id, cleaned_ec)
			#break #remove to read all
	return curvature_torsion_arr

def npy_saver(protein_array):
	'''
	Args:
		protein_array: contains protein objects.
	'''
	path = 'proteins_ca_coordinates/'
	for protein in protein_array:
		arr = np.array([protein.get_ca_coordinates()])
		#print(arr[0])
		file = str(protein.get_uniprot_id()) + '.npy'
		np.save(path+file, arr)

def npy_saver_for_single_protein(protein):
	'''
	Args:
		protein_array: contains protein objects.
	'''
	path = 'proteins_ca_coordinates/'
	arr = np.array([protein.get_ca_coordinates()])
	file = str(protein.get_uniprot_id()) + '.npy'
	np.save(path+file, arr)

def save_ca_coordinates(protein_array):
	'''
	Args:
		protein_array: Array of protein objects.
	Returns:
	'''
	#coordinate_array = []
	for protein in protein_array:
		cleaned_ec = protein.get_cleaned_ec_number()
		pdb_id = read_id_connection(cleaned_ec)
		if pdb_id != None:
			ca_coordinates = read_pdb_file(pdb_id)
			if len(ca_coordinates) > 0:
				protein.set_ca_coordinates(ca_coordinates)
				#coordinate_array.append(protein)
				npy_saver_for_single_protein(protein)
			else:
				print(cleaned_ec)
				#break #remove to read all
		else:
			print(pdb_id, cleaned_ec)
			#break #remove to read all
	#print(coordinate_array)
	#npy_saver(coordinate_array) 

def npy_loader(connection_array):
	'''
	Args:
		connection_array: All protein objects found in the uniprot_sport.dat
	Returns:
		protein_array: Array of protein objects, that have ca_coordinates.
	Raises:
		FileNotFoundError: An error occured when trying to read a file that does
			not exists. Error is expected to happen often and simply passed.
	'''
	path = 'proteins_ca_coordinates/'
	protein_array = []
	for protein in connection_array:
		file = str(protein.get_uniprot_id()) + '.npy'
		try:
			arr = np.load(path+file)
			print(arr[0])
			protein.set_ca_coordinates(arr[0])
			protein_array.append(protein)
		except FileNotFoundError:
			pass
	return protein_array

def set_curvature_and_torsion_from_ca(protein_array):
	'''
	Args:
		protein_array: Array of protein objects.
	Returns:
		curvature_torsion_arr: Array of protein objects that
			have curvature an torsion set.
	'''
	curvature_torsion_arr = []
	for protein in protein_array:
		print(protein.get_uniprot_id())
		ca_coordinates = protein.get_ca_coordinates()
		curvature, torsion = curvature_and_torsion(ca_coordinates)
		protein.set_curvature(curvature)
		protein.set_torsion(torsion)
		curvature_torsion_arr.append(protein)
	return curvature_torsion_arr
