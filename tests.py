from Output_Reader import read_output
import numpy as np


def test_get_excitations_Furan():
    A = read_output.Output_Reader("data/testfiles/Furan.out")
    A.get_excitations()
    
    assert A.dict_excitations["A1"]['excitations'][0] == 5.2122152
    assert A.dict_excitations["A1"]['excitations'][1] == 6.5134867
    assert A.dict_excitations["A1"]['excitations'][2] == 7.1588209
    assert A.dict_excitations["A1"]['type'][0] == 'triplet'
    assert A.dict_excitations["A1"]['type'][1] == 'triplet'
    assert A.dict_excitations["A1"]['type'][2] == 'singlet'
    assert A.dict_excitations['B2']['excitations'][0] == 3.6728915
    assert A.dict_excitations['B2']['excitations'][1] == 6.1429092
    assert A.dict_excitations['B2']['excitations'][2] == 7.2666012
    assert A.dict_excitations['B2']['type'][0] == 'triplet'
    assert A.dict_excitations['B2']['type'][1] == 'singlet'
    assert A.dict_excitations['B2']['type'][2] == 'triplet'
    assert len(A.dict_excitations['B2']['warnings']) == 1
    

def test_get_three_letter_symmetry():
    A = read_output.Output_Reader("data/testfiles/Ethene.out")
    A.get_excitations()
    
    assert 'B1u' in A.dict_excitations.keys()