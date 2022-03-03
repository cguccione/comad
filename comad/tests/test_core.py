import unittest
from skbio.util import get_data_path

from comad.neufit import neufit

class TestCore(unittest.TestCase):

	def setUp(self):
		self.biom_filename = get_data_path('sample_biom')
		self.output_filename = 'sample_name'
		self.output_folder_path = get_data_path('testing_output')
		
		self._data_filename = get_data_path('sample_data.csv')
		self._taxonomy_filename = get_data_path('sample_taxonomy.csv')

	def test_ignore_negative_rasises(self):
		arg_ignore_level = -1
		with self.assertRaises(Exception):
			neufit(self.output_filename, 
				self.output_folder_path, 
				self._data_filename, 
				self._taxonomy_filename,
				arg_ignore_level)

	def test_rarefraction_negative_rasises(self):
		arg_rarefraction_level = -1
		with self.assertRaises(Exception):
			neufit(self.output_filename, 
				self.output_folder_path, 
				self._data_filename, 
				self._taxonomy_filename,
				arg_rarefraction_level)
			


if __name__ == '__main__':
    unittest.main()