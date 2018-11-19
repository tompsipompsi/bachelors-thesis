class Protein:
	def __init__(self, ec_number=None, uniprot_id=None, pdb_id=None):
		self.uniprot_id = uniprot_id
		self.ec_number = ec_number
		self.pdb_id = pdb_id
		if self.ec_number != None:
			self.splitted_ec_number = self.split_ec_number()
		self.curvature = []
		self.torsion = []
		self.ca_coordinates = []

	def get_uniprot_id(self):
		return self.uniprot_id

	def get_ec_number(self):
		return self.ec_number

	def is_integer(self, num):
		try: 
			int(num)
			return True
		except ValueError:
			return False

	def split_ec_number(self):
		#removes non-integer values from the ec_number and returns it as a array of strings
		splitted_ec_number = []
		for value in self.ec_number.split("."):
			if self.is_integer(value):
				splitted_ec_number.append(value)
		return splitted_ec_number

	def get_splitted_ec_number(self):
		return self.splitted_ec_number

	def get_cleaned_ec_number(self):
		ec_str = ""
		i = 0
		length = len(self.splitted_ec_number)
		#print(self.splitted_ec_number)
		for char in self.splitted_ec_number:
			ec_str += char
			if i < length-1:
				ec_str += "."
			i += 1
		return ec_str

	def set_pdb_id(self, pdb_id):
		self.pdb_id = pdb_id

	def set_curvature(self, curvature):
		self.curvature = curvature

	def set_torsion(self, torsion):
		self.torsion = torsion

	def set_ca_coordinates(self, ca_coordinates):
		self.ca_coordinates = ca_coordinates

	def get_curvature(self):
		return self.curvature

	def get_torsion(self):
		return self.torsion

	def get_ca_coordinates(self):
		return self.ca_coordinates

	def print_protein(self):
		print(self.ec_number, self.uniprot_id, self.curvature, self.torsion)

	"""def get_as_array(self):
		return [self.ec_number, self.uniprot_id, self.curvature, self.torsion]

	def get_as_str(self):
		string = ""
		string += self.ec_number + ',' + self.uniprot_id + ','
		if len(self.curvature) > 0:
			string = string.join(str(i) for i in self.curvature)
			string += ','
		if len(self.torsion) > 0:
			string = string.join(str(i) for i in self.torsion)
		return string"""