'''
Some functions to create the binary vector.
'''
import numpy as np

def find_greatest_ec_values(protein_array):
	'''
	Args:
		protein_array: Array of protein objects.
	Returns:
		binary_vector: Numpy array containing four subarrays with
			the lenght of largest EC numbers found.
	'''
	vector1, vector2, vector3, vector4 = [], [], [], []
	for protein in protein_array:
		i = 0
		for ec_number in protein.get_splitted_ec_number():
			number = int(ec_number)
			if i == 0:
				vector1.append(number)
			elif i == 1:
				vector2.append(number)
			elif i == 2:
				vector3.append(number)
			elif i == 3:
				vector4.append(number)
			i += 1
	binary_vector = create_binary_vector(max(vector1), max(vector2), max(vector3), max(vector4))
	return binary_vector

def create_binary_vector(len1, len2, len3, len4):
	'''
	Args:
		len1, len2, len3, len4: integers
	Returns:
		binary_vector: Numpy array with four subarrays of
			lengths: len1, len2, len3, len4.
	'''
	vec1 = np.zeros((len1,), dtype=int)
	vec2 = np.zeros((len2,), dtype=int)
	vec3 = np.zeros((len3,), dtype=int)
	vec4 = np.zeros((len4,), dtype=int)
	binary_vector = np.array([vec1, vec2, vec3, vec4])
	return binary_vector 

def insert_ec_into_binary_vector(binary_vector, protein_array):
	'''
	Args:
		binary_vector: Numpy array.
		protein_array: Array of protein objects.
	Returns:
		binary_vector: Numpy array with found EC numbers inserted into it.
	'''
	for protein in protein_array:
		ec_number = protein.get_splitted_ec_number()
		for i in range(len(ec_number)):
			ec = int(ec_number[i]) - 1 #to make ec compatible with indexes
			if i == 0:
				binary_vector[0][ec] = 1
			elif i == 1:
				binary_vector[1][ec] = 1
			elif i == 2:
				binary_vector[2][ec] = 1
			elif i == 3:
				binary_vector[3][ec] = 1
	return binary_vector