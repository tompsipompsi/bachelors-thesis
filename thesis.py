import numpy as np
from datetime import datetime
from Bio import SwissProt
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import pylab 	

from protein import Protein
#import binary_vector as bv
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

def get_mean_std_data(protein_array):
	'''
	Args:
		protein_array: Array containing protein objects.
	Returns:
		Mean and standard deviation of the curvature and torsion of all proteins.
	'''
	c_mean, c_std, t_mean, t_std = [], [], [], []
	for protein in protein_array:
		c_mean.append(np.mean(protein.get_curvature()))
		c_std.append(np.std(protein.get_curvature()))
		t_mean.append(np.mean(protein.get_torsion()))
		t_std.append(np.std(protein.get_torsion()))
	return np.mean(c_mean), np.mean(c_std), np.mean(t_mean), np.mean(t_std)

def normalize_data(protein_array, bins):
	'''
	Args:
		protein_array: Array containing protein objects.
	Returns:
		
	'''
	c_mean, c_std, t_mean, t_std = get_mean_std_data(protein_array)
	for protein in protein_array:
		val = protein.get_torsion()
		ival = np.where(val > t_mean + 3*t_std)[0]
		val[ival] = t_mean + 3*t_std

		ival = np.where(val < t_mean - 3*t_std)[0]
		val[ival] = t_mean - 3*t_std
		protein.set_torsion(val)

		c_val = protein.get_curvature()
		ival = np.where(c_val > c_mean + 3*c_std)[0]
		c_val[ival] = c_mean + 3*c_std
		protein.set_curvature(c_val)
	c_min = 0
	c_max = c_mean + 3*c_std
	t_min = t_mean - 3*t_std
	t_max = t_mean + 3*t_std

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
	for protein in protein_array:
		curvature = protein.get_curvature()
		torsion = protein.get_torsion()
		norm_curvature = normalize_values(curvature, c_min, c_max, bins)
		norm_torsion = normalize_values(torsion, t_min, t_max, bins)
		c_histog = np.histogram(norm_curvature, range(bins))[0]
		t_histog = np.histogram(norm_torsion, range(bins))[0]
		#protein.set_feature_vector(c_histog+t_histog)
		arr.append(np.concatenate([c_histog, t_histog]))
		print(protein.get_ec_number())
	arr = np.array(arr)
	return arr


def fit_data(X, Y_i):
	X = X - np.outer(np.sum(X,1), np.ones(X.shape[1]))
	x_lambda = 1
	#w_opt = np.dot(np.dot(np.linalg.pinv(np.dot(X.T, X)+x_lambda*np.eye(X.shape[1])), X.T), Y_i) 
	w_opt = np.dot(np.linalg.pinv(np.dot(X.T, X)+x_lambda*np.eye(X.shape[1])), np.dot(X.T, Y_i))
	return w_opt
'''
def predict(x_test, W): 
	y_pred = np.dot(x_test.T, W)
	return y_pred
'''
def ridge_regression(protein_array, proteins_with_features, binary_vector, nmax):
	x_train, x_test, y_train, y_test = train_test_split(proteins_with_features, binary_vector, test_size=0.5)
	w_opt = fit_data(x_train, y_train)

	#y_pred = predict(x_test, w_opt)

	y_test_i = np.dot(w_opt.T, x_test.T) 
	
	val = np.dot(y_train, y_test_i)
	
	index = np.argmax(val, 0)
	#pylab.plot(index)
	#pylab.show()
	#pylab.hist() tai jtn tonnepäin

	pred = y_train[index]

	def f_pred(pred, y_test):
		tp = np.sum((pred==1) * (y_test==1)) #true positives
		fp = np.sum((pred==1) * (y_test==0)) #false positives
		fn = np.sum((pred==0) * (y_test==1)) #false negatives
		precision = tp/(tp+fp)
		recall = tp/(tp+fn)
		f_1 = 2*precision*recall/(precision+recall) #f1-score

		print(precision, recall, f_1)

	for i in range(4): #4 is the maximum length of ec number
		pred_1 = pred[:,nmax[i]:nmax[i+1]]
		y_test_1 = y_test[:,nmax[i]:nmax[i+1]]
		f_pred(pred_1, y_test_1)

	#y_i = y_i - np.outer(np.ones(y_i.shape[0]), np.mean(y,0))
	#y_pred_w_index = y_pred_w_index - np.outer(np.ones(y_pred_w_index.shape[0]), np.mean(y,0))

	#np.corrcoef(vec(y_i), vec(y_pred_w_index)) #=> 2x2 matrix, variables not right

	#rmse = np.sqrt(np.mean(y_test-pred)**2) 
	#rmse = np.sqrt(np.sum(np.square(y_test - y_pred))/len(y_test))
	#print(rmse)

def create_binary_vector(protein_curvature_torsion_arr):
	ec_array = []
	for protein in protein_curvature_torsion_arr:
		ec = protein.get_splitted_ec_number()
		ec_number = [0 for _ in range(4)] #4 is the maximum length of ec number
		ec_number[0] = int(ec[0])
		if len(ec) > 1:
			ec_number[1] = int(ec[1])
		if len(ec) > 2:
			ec_number[2] = int(ec[2])
		if len(ec) > 3:
			ec_number[3] = int(ec[3])
		ec_array.append(ec_number)
	ec_array = np.array(ec_array)
	
	xmax = np.max(ec_array, 0)
	nmax = np.concatenate([[0], np.cumsum(xmax)])
	m, n = ec_array.shape
	binary_vector = np.zeros((m, nmax[-1]))
	for i in range(m):
		for j in range(n):
			if ec_array[i,j] > 0:
				binary_vector[i, ec_array[i,j] - 1 + nmax[j]] = 1
	return binary_vector, nmax

def main():
	start_time = datetime.now() #Datetime for benchmarking

	bins = 1000 #1500 works

	'''
	Connection array is used to store the EC number and UniProt id.
	It is saved as .csv, and loaded directly into protein objects.
	Array of protein objects can also be directly obtained from reader.
	'''
	#connection_array, protein_array = read_uniprot_sequence()
	#fh.csv_writer('general_files/', 'connection_array.csv', connection_array)
	protein_array = fh.csv_loader('general_files/', 'connection_array.csv')

	'''
	CA stuff. Optional with the next section.
	'''
	cth.save_ca_coordinates(protein_array) #14min
	ca_array = fh.npy_loader('proteins_ca_coordinates/', protein_array, is_ca=True)
	#protein_curvature_torsion_arr = cth.set_curvature_and_torsion_from_ca(ca_array) #9min
	
	'''
	Straight from protein array
	'''
	#protein_curvature_torsion_arr = cth.set_curvature_and_torsion(protein_array) #29min, 22.5min
	#fh.npy_saver('proteins_curvature_torsion/', protein_curvature_torsion_arr)
	protein_curvature_torsion_arr = fh.npy_loader('proteins_curvature_torsion/', protein_array)

	'''
	Binary vector
	'''
	#binary_vector, nmax = create_binary_vector(protein_curvature_torsion_arr)
	#np.save('general_files/binary_vector', binary_vector)
	binary_vector = np.load('general_files/binary_vector.npy')

	'''
	Machine learning
	'''
	#proteins_with_features = normalize_data(protein_curvature_torsion_arr, bins)
	#ridge_regression(protein_curvature_torsion_arr, proteins_with_features, binary_vector, nmax)
	
	end_time = datetime.now()
	print("Start time: ", start_time, " Finish time: ", end_time)
	print("Runtime: " , end_time - start_time)

if __name__ == '__main__':
	main()
