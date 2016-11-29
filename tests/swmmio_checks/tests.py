import unittest
import os
from swmmio.swmmio import Model

class CreateDictionaryTest(unittest.TestCase):
    def test_inp_create_dictionary(self):
        """
        make sure the createDictionary method is return the proper data. check
        against a sample model with 380 junctions.
        """
        dir_path = os.path.dirname(os.path.realpath(__file__))
        test_model = Model(os.path.join(dir_path, 'sample_subset.inp'))

        juncs_dict = test_model.inp.createDictionary('[JUNCTIONS]')
        subcats_dict = test_model.inp.createDictionary('[SUBCATCHMENTS]')

        self.assertEqual(len(juncs_dict), 380)
        self.assertEqual(len(subcats_dict), 355)


if __name__ == '__main__':
    unittest.main()
