import numpy as np
from datetime import datetime
from Bio import SwissProt
import matplotlib.pyplot as plt

from protein import Protein
import binary_vector as bv
import curvature_torsion_handlers as cth
import file_handlers as fh

def read_uniprot_sequence():
	'''
	Returns:
		connection_array:
		protein_array:
	'''
	file = 'uniprot_sprot.dat'
	#counter = 0 #temporary for testing, no need to read the whole file yet
	protein_array = []
	ec_array = []
	connection_array = []
	for record in SwissProt.parse(open(file)):
		if 'EC=' in record.description:
			#counter += 1
			#sequence is the string of the primary sequence
			#given by the markers of the residues
			print(record.sequence)
			#print(record.accessions) #holds the uniprot ids

			#description consists of ';' separated parts
			print(record.description)
			tokens = record.description.split(';')

			for token in tokens:
				if 'EC=' in token:
					parts = token.split('=')		#split header
					ec_parts = parts[1].split(' ')	#split additional content
					if ec_parts[0] not in ec_array:
						print('EC: ', ec_parts[0])		#print EC number as a string
						ec_array.append(ec_parts[0])
						connection_array.append([ec_parts[0], record.accessions[0]])
						protein_array.append(Protein(ec_parts[0], record.accessions[0]))
			#if counter >= 10000:
			#	break
	return connection_array, protein_array

def get_min_max_data(protein_array):
	'''
	Args:
		protein_array: Array containing protein objects.
	Returns:
		Minimum and maximum values of the curvature and torsion of all proteins.
	'''
	c_max, c_min, t_max, t_min = [], [], [], []
	for protein in protein_array:
		c_max.append(max(protein.get_curvature()))
		c_min.append(min(protein.get_curvature()))
		t_max.append(max(protein.get_torsion()))
		t_min.append(min(protein.get_torsion()))
	return min(c_min), max(c_max), min(t_min), max(t_max)

def normalize_data(protein_array, binary_vector):
	'''
	Args:
		protein_array: Array containing protein objects.
	Returns:
		
	'''
	c_min, c_max, t_min, t_max = get_min_max_data(protein_array)
	bins = 100
	def normalize_values(array, min_s, max_s, bins):
		'''
		Args:
			array: Contains array of values to be normalized.
			min_s: Smallest value of all the data (curvature or torsion).
			max_s: Largest value of all the data (curvature or torsion).
			bins: Amount of bins which normalized values are scaled.
		Returns:
			Numpy array of normalized values of the input array.
		'''
		norm_data = []
		for sample in array:
			norm_sample = int((sample-min_s)/(max_s-min_s)*bins)
			norm_data.append(norm_sample)
		return np.array(norm_data)
	arr = []
	vector = []
	for protein in protein_array:
		curvature = protein.get_curvature()
		torsion = protein.get_torsion()
		norm_curvature = normalize_values(curvature, c_min, c_max, bins)
		norm_torsion = normalize_values(torsion, t_min, t_max, bins)
		c_histog = np.histogram(norm_curvature, range(bins))[0]
		t_histog = np.histogram(norm_torsion, range(bins))[0]
		#print(norm_curvature, norm_torsion)
		print(c_histog+t_histog)
		#plt.hist(histog, bins)
		#plt.hist(norm_torsion, bins)
		#plt.show()
		#print(c_min,c_max)
		#protein.set_feature_vector(c_histog+t_histog)
		arr.append(c_histog+t_histog)
		vector.append(binary_vector)
	arr = np.array(arr)
	vec = np.array(vector)
	return arr, vec


def fit_data(X, Y_i):
	X = X - np.outer(np.sum(X,1), np.ones(X.shape[1]))
	x_lamba = 1
	#X_t = np.transpose(X)
	W = np.dot(np.dot(np.linalg.pinv(np.dot(X.T, X)+x_lamba*np.eye(X.shape[0])), X.T), Y_i)
	return W #*X-text

def linear_regression(protein_array, proteins_with_features, binary_vector):
	#for protein in protein_array:
	print(fit_data(proteins_with_features, binary_vector))
		#print(protein.get_feature_vector())
		#break

#rmse = np.sqrt(np.mean(y_i-y)**2)
#y_i = y_i - np.outer(np.ones(y.shape[0]), np.mean(y,0))
#=> same for prediction y_p

#np.corrcoef(vec(y_i), vec(y_p)) => 2x2 matrix
#Y = training data

def main():
	start_time = datetime.now() #Datetime for benchmarking

	#connection_array, protein_array = read_uniprot_sequence() #over 5min
	#connection_array = read_uniprot_sequence()
	#fh.csv_writer('new_connection_array.csv', connection_array)

	protein_array = fh.csv_loader('connection_array.csv', True)

	#binary_vector = bv.find_greatest_ec_values(protein_array)
	#binary_vector = bv.insert_ec_into_binary_vector(binary_vector, protein_array)
	#np.save('general_files/binary_vector', binary_vector) #later no need to do the vector every time again

	binary_vector = np.load('general_files/binary_vector.npy')

	#cth.save_ca_coordinates(protein_array) #16min
	#ca_array = cth.npy_loader(protein_array)
	#protein_curvature_torsion_arr = cth.set_curvature_and_torsion_from_ca(ca_array) #9min
	
	#protein_curvature_torsion_arr = cth.set_curvature_and_torsion(protein_array) #29min, 22.5min
	#fh.npy_saver(protein_curvature_torsion_arr)
	protein_curvature_torsion_arr = fh.npy_loader(protein_array) #protein_array = connection array

	#for protein in protein_curvature_torsion_arr:
	#	print(protein.get_ec_number(), protein.get_uniprot_id())
	#print(len(protein_array), len(protein_curvature_torsion_arr)) 
	
	proteins_with_features, vector = normalize_data(protein_curvature_torsion_arr, binary_vector)
	
	linear_regression(protein_curvature_torsion_arr, proteins_with_features, vector)

	end_time = datetime.now()
	print("Start time: ", start_time, " Finish time: ", end_time)
	print("Runtime: " , end_time - start_time)

if __name__ == '__main__':
	main()
