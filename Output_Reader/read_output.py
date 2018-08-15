import numpy as np

class Output_Reader():
    def __init__(self, output_file):
        with open(output_file, "r", encoding="utf-8") as f:
            self.__load_file = list(f)
            
    def get_excitations(self):
        self.dict_excitations = {}
        found_rs_index_check = 0
        excitation_found_check = 0
        for line in self.__load_file:
            if "triplet =" and "Operator symmetry" in line:
                if "triplet =   T" in line:
                    default_type = "triplet"
                else:
                    default_type = "singlet"
                symmetry = line.split()[5]
                if ")" in symmetry:
                    symmetry = symmetry.split(")")[0]
                self.dict_excitations[symmetry] = {"excitations":[],
                                                   "type":[],
                                                   "warnings":[]}
                excitation_found_check = 1
            elif "WARNING Complex eigenvalue" in line:
                if "Complex eigenvalues in Linear Response equations" not in self.dict_excitations[symmetry]["warnings"]:
                    self.dict_excitations[symmetry]["warnings"].append("Complex eigenvalues in Linear Response equations")
            elif "SOLUTION VECTORS NOT CONVERGED" in line:
                self.dict_excitations[symmetry]["warnings"].append("Linear Response equations did not converge")
            elif "@ Excited state no:" in line:
                excitation_found_check = 1
            elif "kJ / mol" in line and excitation_found_check == 1:
                excitation_energy = float(line.split()[1])
                self.dict_excitations[symmetry]["excitations"].append(excitation_energy)
            elif "-----  -----" in line and excitation_found_check == 1:
                found_rs_index_check = 1
                total_rs = 0
            elif found_rs_index_check == 1 and line == "\n":
                found_rs_index_check = 0
            elif found_rs_index_check == 1:
                total_rs += abs(float(line.split()[3]))
            elif excitation_found_check == 1 and "The numbers in paren" in line:
                if total_rs < 10**-6 and default_type == "triplet":
                    self.dict_excitations[symmetry]["type"].append("singlet")
                else:
                    self.dict_excitations[symmetry]["type"].append(default_type)
                excitation_found_check = 0
                found_rs_index_check = 0
                    
    def get_natural_occupations(self):
        None
            