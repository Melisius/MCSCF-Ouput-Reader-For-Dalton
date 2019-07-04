import numpy as np

class Output_Reader():
    def __init__(self, output_file):
        with open(output_file, "r", encoding="utf-8") as f:
            self.__load_file = list(f)
            
    def get_excitations(self):
        self.dict_excitations = {}
        found_rs_index_check = 0
        excitation_found_check = 0
        excitation_found_check = 0
        cas_found_check = 0
        ci_found_check = 0
        ci_groundstate_check = 0
        groundstate_ci = [0]
        alpha = "None"
        beta = "None"
        for line in self.__load_file:
            if ".CAS" in line:
                cas_found_check = 1
            elif cas_found_check == 1:
                cas = line.split()
                cas = [int(i) for i in cas]
                cas_idx = []
                for j in range(0, len(cas)):
                    for i in range(cas[j]):
                        cas_idx.append(j+1)
                cas_found_check = 0
            elif "Printout of coefficients in interval" and ci_groundstate_check == 0:
                ci_groundstate_check = 1
            elif ci_groundstate_check == 1 and "alpha-string" in line:
                groundstate_ci = line.split()[1:]
                groundstate_ci = [int(i) for i in groundstate_ci]
                ci_groundstate_check = 2
            elif "triplet =" and "Operator symmetry" in line:
                if "triplet =   T" in line:
                    default_type = "triplet"
                else:
                    default_type = "singlet"
                symmetry = line.split()[5]
                if ")" in symmetry:
                    symmetry = symmetry.split(")")[0]
                self.dict_excitations[symmetry] = {"excitations":[],
                                                   "type":[],
                                                   "orb_classification":[],
                                                   "ci_classification":[],
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
                # sym_x --> sym_y
                rs_symmetry = np.zeros((9,9))
            elif found_rs_index_check == 1 and line == "\n":
                found_rs_index_check = 0
            elif found_rs_index_check == 1:
                total_rs += abs(float(line.split()[3]))
                r_idx = int(line.split()[1].split("(")[1].split(")")[0])
                s_idx = int(line.split()[2].split("(")[1].split(")")[0])
                rs_symmetry[r_idx,s_idx] += float(line.split()[5])**2
            elif excitation_found_check == 1 and "The numbers in paren" in line:
                rs_symmetry = rs_symmetry**0.5
                if total_rs < 10**-6 and default_type == "triplet":
                    self.dict_excitations[symmetry]["type"].append("singlet")
                else:
                    self.dict_excitations[symmetry]["type"].append(default_type)
                if np.sum(rs_symmetry) != 0.0:
                    temp = []
                    rs_symmetry_tmp = np.copy(rs_symmetry)
                    while np.sum(rs_symmetry_tmp) > 0.0:
                        temp.append([np.where(rs_symmetry_tmp==rs_symmetry_tmp.max())[0][0],np.where(rs_symmetry_tmp==rs_symmetry_tmp.max())[1][0], np.max(rs_symmetry_tmp)])
                        rs_symmetry_tmp[np.where(rs_symmetry_tmp==rs_symmetry_tmp.max())[0][0],np.where(rs_symmetry_tmp==rs_symmetry_tmp.max())[1][0]] = 0.0
                    self.dict_excitations[symmetry]["orb_classification"].append(temp)
                else:
                    self.dict_excitations[symmetry]["orb_classification"].append([[0,0,0.0]])
                excitation_found_check = 0
                found_rs_index_check = 0
                ci_found_check = 1 
                ci_coeff = np.zeros((9,9))
            elif ci_found_check == 1 and "Coefficient of determinant" in line:
                coeff = float(line.split()[5])
            elif ci_found_check == 1 and "alpha-string" in line:
                alpha = line.split()[1:]
            elif ci_found_check == 1 and "beta-string" in line:
                beta = line.split()[1:]
            elif beta != "None" and alpha != "None":
                alpha = [int(i) for i in alpha]
                beta = [int(i) for i in beta]
                det = []
                for i in alpha:
                    if i not in groundstate_ci:
                        det.append(i)
                for i in groundstate_ci:
                    if i not in alpha:
                        det.append(-i)
                for i in beta:
                    if i not in groundstate_ci:
                        det.append(i)
                for i in groundstate_ci:
                    if i not in beta:
                        det.append(-i)
                det.sort()
                if len(det) == 2:
                    det = [int(abs(i)) for i in det]
                    idx1 = abs(cas_idx[int(det[0])-1])
                    idx2 = abs(cas_idx[int(det[1])-1])
                    ci_coeff[idx1,idx2] += coeff**2
                alpha = "None"
                beta = "None"
            elif "Magnitude of CI coefficients" in line and ci_found_check == 1:
                ci_coeff = ci_coeff**0.5
                if np.sum(ci_coeff) != 0.0:
                    temp = []
                    ci_coeff_tmp = np.copy(ci_coeff)
                    while np.sum(ci_coeff_tmp) > 0.0:
                        temp.append([np.where(ci_coeff_tmp==ci_coeff_tmp.max())[0][0],np.where(ci_coeff_tmp==ci_coeff_tmp.max())[1][0], np.max(ci_coeff_tmp)])
                        ci_coeff_tmp[np.where(ci_coeff_tmp==ci_coeff_tmp.max())[0][0],np.where(ci_coeff_tmp==ci_coeff_tmp.max())[1][0]] = 0.0
                    self.dict_excitations[symmetry]["ci_classification"].append(temp)
                else:
                    self.dict_excitations[symmetry]["ci_classification"].append([[0,0,0.0]])
                ci_found_check = 0
            elif ">> NO ELEMENTS <<" in line and ci_found_check == 1 :
                self.dict_excitations[symmetry]["ci_classification"].append([[0,0,0.0]])
                ci_found_check = 0
                


    def get_natural_occupations(self):
        None
    

    def get_mp2srdft_contribution(self):
        self.mp2_energy = 0.0
        self.hf_energy = 0.0
        for line in self.__load_file:
            if "@   Short-range Hartree-Fock total energy        :" in line:
                self.hf_energy = float(line.split(":")[1])
            elif "@   + MP2 contribution from long-range integrals :" in line:
                self.mp2_energy = float(line.split(":")[1])
                break
            elif "Final results from SIRIUS" in line:
                break
    
                
    def get_mcsrdft_contribution(self):
        self.mp2_energy = 0.0
        self.hf_energy = 0.0
        for line in self.__load_file:
            if "@   Short-range Hartree-Fock total energy        :" in line:
                self.hf_energy = float(line.split(":")[1])
            elif "@    Final MC-SRDFT energy:" in line:
                self.mcscf_energy = float(line.split(":")[1])
                break
            elif "Final results from SIRIUS" in line:
                break
                
                
    def get_mcscf_contribution(self):
        self.mp2_energy = 0.0
        self.hf_energy = 0.0
        for line in self.__load_file:
            if "@   Hartree-Fock total energy   :" in line:
                self.hf_energy = float(line.split(":")[1])
            elif "Final MCSCF energy:" in line:
                self.mcscf_energy = float(line.split(":")[1])
                break
            elif "Final results from SIRIUS" in line:
                break
    

    def get_mp2_contribution(self):
        self.mp2_energy = 0.0
        self.hf_energy = 0.0
        for line in self.__load_file:
            if "@   Hartree-Fock total energy   :" in line:
                self.hf_energy = float(line.split(":")[1])
            elif "@   + MP2 contribution          :" in line:
                self.mp2_energy = float(line.split(":")[1])
                break
            elif "Final results from SIRIUS" in line:
                break
    

    def get_hfsrdft_contribution(self):
        self.hf_energy = 0.0
        self.srExc_energy = 0.0
        self.srEJ_energy = 0.0
        for line in self.__load_file:
            if "Final HF-SRDFT energy:" in line:
                self.hf_energy = float(line.split(":")[1])
            elif "Ex-sr + Ec-sr" in line:
                self.srExc_energy = float(line.split("Ec-sr")[1])
            elif "+ EJsr = sr" in line:
                self.srEJ_energy = float(line.split("energy")[1])
            elif "!!!Final results from SIRIUS" in line:
                break
            

    def get_dft_contribution(self):
        self.dft_energy = 0.0
        for line in self.__load_file:
            if "@    Final DFT energy:" in line:
                self.dft_energy = float(line.split(":")[1])
                break


    def get_spin_spin_coupling_constants(self):
        self.sscc = {}
        coupling_found = False
        for line in self.__load_file:
            if "Indirect spin-spin coupling between" in line:
                if coupling_found == True:
                    self.sscc[atom1+"//"+atom2] = values[:couplings_counter,:]
                atom1 = line.split("between")[1].split("and")[0].replace(" ","").replace("_","  _")
                atom2 = line.split("between")[1].split("and")[1].replace(" ","").replace("_","  _").replace(":\n","")
                values = np.zeros((10,11))
                couplings_counter = 0
                coupling_found = True
            elif coupling_found == True:
                if "Mass number atom 1:" in line:
                    values[couplings_counter,0] = float(line.split("Abundance:")[1].split("%")[0])
                if "Mass number atom 2:" in line:
                    values[couplings_counter,1] = float(line.split("Abundance:")[1].split("%")[0])
                elif "Isotropic coupling" in line:
                    try:
                        values[couplings_counter,2] = float(line.split(":")[1].split("H")[0])
                    except:
                        values[couplings_counter,2] = False
                elif "Anisotropic coupling" in line:
                    try:
                        values[couplings_counter,3] = float(line.split(":")[1].split("H")[0])
                    except:
                        values[couplings_counter,3] = False
                elif "Asymmetry" in line:
                    try:
                        values[couplings_counter,4] = float(line.split(":")[1])
                    except:
                        values[couplings_counter,4] = False
                elif "S parameter" in line:
                    try:
                        values[couplings_counter,5] = float(line.split(":")[1].split("H")[0])
                    except:
                        values[couplings_counter,5] = False
                elif "A parameter" in line:
                    try:
                        values[couplings_counter,6] = float(line.split(":")[1].split("H")[0])
                    except:
                        values[couplings_counter,6] = False
                elif "Isotropic DSO contribution" in line:
                    try:
                        values[couplings_counter,7] = float(line.split(":")[1].split("H")[0])
                    except:
                        values[couplings_counter,7] = False
                elif "Isotropic PSO contribution" in line:
                    try:
                        values[couplings_counter,8] = float(line.split(":")[1].split("H")[0])
                    except:
                        values[couplings_counter,8] = False
                elif "Isotropic SD contribution" in line:
                    try:
                        values[couplings_counter,9] = float(line.split(":")[1].split("H")[0])
                    except:
                        values[couplings_counter,9] = False
                elif "Isotropic FC contribution" in line:
                    try:
                        values[couplings_counter,10] = float(line.split(":")[1].split("H")[0])
                    except:
                        values[couplings_counter,10] = False
                    couplings_counter += 1
            if "End of Static Property Section" in line and coupling_found == True:
                self.sscc[atom1+"//"+atom2] = values[:couplings_counter,:]
                break
