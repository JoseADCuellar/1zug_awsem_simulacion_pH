try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except ModuleNotFoundError:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
import numpy as np

def read_fasta(fastaFile):
    with open(fastaFile) as input_data:
        data = ""
        for line in input_data:
            if(line[0] == ">"):
                print(line)
            elif(line == "\n"):
                pass
            else:
                data += line.strip("\n")
    return data

def debye_huckel_term(self, k_dh=15*4.184, forceGroup=30, screening_length=1.0, chargeFile=None,periodic=False):#kdh_lammps 3.357
        # screening_length (in the unit of nanometers)
        print("Debye Huckel term is ON,yes MAM")
        k_dh *= self.k_awsem*0.1
        print("k_dh",k_dh)
        k_screening = 1.0
        # screening_length = 1.0  # (in the unit of nanometers)
        min_seq_sep = 1
        dh = CustomBondForce(f"({k_dh}*charge_i*charge_j/r)*exp(-{k_screening}*r/{screening_length})") #orginal
        dh.addPerBondParameter("charge_i")
        dh.addPerBondParameter("charge_j")

         # Add Periodic Boundary Condition. 02082024 Rebekah Added --- Start
        if periodic:
            dh.setUsesPeriodicBoundaryConditions(True)
            is_periodic=dh.usesPeriodicBoundaryConditions()
            print("\ndebye_huckel_term is in PBC",is_periodic)
         # 02082024 Rebekah Added --- End
        
        
        structure_interactions_dh = []
        if chargeFile is None:
            for i in range(self.nres):
                for j in range(i+min_seq_sep,self.nres):
                    charge_i = 0.0
                    charge_j = 0.0
                    if self.seq[i] == "R" or self.seq[i]=="K":
                        charge_i = 1.0
                    if self.seq[i] == "D" or self.seq[i]=="E":
                        charge_i = -1.0
                    if self.seq[j] == "R" or self.seq[j]=="K":
                        charge_j = 1.0
                    if self.seq[j] == "D" or self.seq[j]=="E":
                        charge_j = -1.0
                    if charge_i*charge_j!=0.0:
                        cb_atom_i = self.cb[i]
                        if cb_atom_i == -1:
                            cb_atom_i = self.ca[i]  # if mutated, and CB isn't found, then use CA instead
                        cb_atom_j = self.cb[j]
                        if cb_atom_j == -1:
                            cb_atom_j = self.ca[j]  # if mutated, and CB isn't found, then use CA instead
                        structure_interactions_dh.append([cb_atom_i, cb_atom_j, [charge_i, charge_j]])
                        # print([self.seq[i], self.seq[j],self.cb[i], self.cb[j], [charge_i, charge_j]])
        else:
            chargeInfo = np.loadtxt(chargeFile, dtype=[('index', int), ('charge', float)])
            for i in range(self.nres+2):# mas 2 por los terminales extras si hay 105 residuso por el range va de 0 a 104
                #print(self.nres + 2) le suma los dos para le terminal que para la frata 4hs5A 107 los terminales son 105 106 self.nres=105
                charge_i = chargeInfo[i][1]
                for j in range(i+min_seq_sep,self.nres+2):               
                    charge_j = chargeInfo[j][1]
                    if charge_i*charge_j!=0.0:
                        if i == self.nres + 1: #cterminal
                            cb_atom_i = self.o[-1]
                        elif i == self.nres:# nterminal
                            cb_atom_i = self.ca[0]
                        else:
                            cb_atom_i = self.cb[i]
                            if cb_atom_i == -1:
                            	cb_atom_i = self.ca[i]    
                        if j == self.nres + 1: #cterminal
                            cb_atom_j = self.o[-1]
                        elif j == self.nres: #nterminal
                            cb_atom_j = self.ca[0]
                        else:
                            cb_atom_j = self.cb[j]
                            if cb_atom_j == -1:
                            	cb_atom_j = self.ca[j]  # if mutated, and CB isn't found, then use CA instead
                            
                        structure_interactions_dh.append([cb_atom_i, cb_atom_j, [charge_i, charge_j]])
        for structure_interaction_dh in structure_interactions_dh:
            dh.addBond(*structure_interaction_dh)
            #Change####
          #  charge_info = f"Aminoácido {i} (índice {chargeInfo[i][0]}), carga: {charge_i}, Aminoácido {j} (índice {chargeInfo[j][0]}), carga: {charge_j}"
           # print(charge_info)
            #charge_file.write(charge_info + "\n")
            #Change####

        dh.setForceGroup(forceGroup)
        return dh
