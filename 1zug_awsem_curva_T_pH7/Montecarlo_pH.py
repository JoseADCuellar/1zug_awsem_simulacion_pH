import numpy as np
import pandas as pd
import random
import math

def convert_sequence_to_three_letter(seq):
    # Diccionario de mapeo de una letra a tres letras
    one_to_three = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
        'C': 'CYS', 'E': 'GLU', 'Q': 'GLN', 'G': 'GLY',
        'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
        'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
        'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
        'X': 'NTR', 'Z': 'CTR',
    }
    # Convertir la secuencia de una letra a tres letras
    seq_list = [one_to_three[aa] for aa in seq]
    return seq_list

def calcular_distancia(punto1, punto2):
    return np.sqrt((punto1[0] - punto2[0])**2 + (punto1[1] - punto2[1])**2 + (punto1[2] - punto2[2])**2)

def charge_flip(aob, old_charge):
    new_charge = 0.0
    if aob == -1.0:  # Ácido
        if old_charge == -1.0:
            new_charge = 0.0
        elif old_charge == 0.0:
            new_charge = -1.0
    elif aob == 1:  # Base
        if old_charge == 1.0:
            new_charge = 0.0
        elif old_charge == 0.0:
            new_charge = 1.0
    return new_charge

class Montecarlo_change_charge:
    
    def __init__(self, matriz_posiciones, charged_residues, seq, enlaces_matrix, pH):
        self.matriz_posiciones = matriz_posiciones
        self.charged_residues = charged_residues
        self.seq_list = convert_sequence_to_three_letter(seq)
        self.type_dict = {
        'ALA': 'N', 'ARG': 'B', 'ASN': 'P', 'ASP': 'A', 'CYS': 'A', 'GLU': 'A', 'GLN': 'P', 'GLY': 'G', 'HIS': 'B',
        'ILE': 'N', 'LEU': 'N', 'LYS': 'B', 'MET': 'N', 'PHE': 'N', 'PRO': 'N', 'SER': 'P', 'THR': 'P', 'TRP': 'N',
        'TYR': 'A', 'VAL': 'N', 'NTR': 'B', 'CTR': 'A',
        }
        self.pH = pH
        self.enlaces_matrix = enlaces_matrix
        
    def choose_residue_mc(self):
        self.residue_mc, self.charge = random.choice(self.charged_residues)

    def acid_basic_residue_mc(self):
        self.resname_residue_mc = self.seq_list[self.residue_mc]     
        self.acid_basic_residue_mc = 1 if self.type_dict[self.resname_residue_mc] == 'B' else -1
    
    def distances_residue_mc_neighbors(self):
        # Obtener la posición del residuo de interés
        self.position_residue_mc = self.matriz_posiciones[self.matriz_posiciones['Residue_Number'] == self.residue_mc][['X', 'Y', 'Z']].values[0]
        # Calcular las distancias a todos los demás residuos
        self.neighbors = []
        self.distances = []
        for index, row in self.matriz_posiciones.iterrows():
            if row['Residue_Number'] != self.residue_mc:
                self.position_neighbor = [row['X'], row['Y'], row['Z']]
                self.distance = calcular_distancia(self.position_residue_mc, self.position_neighbor)
                self.neighbors.append(row['Residue_Number'])  # Cambio de nombre a 'vecinos'
                self.distances.append(self.distance)
        return self.distances, self.neighbors
        
    def get_resnames_neighbors(self):
        self.distances_residue_mc_neighbors()
        
        self.neighbors = [int(neighbor) for neighbor in self.neighbors]
        # Obtener los nombres de los residuos correspondientes a los índices en 'vecinos'
        self.resnames_neighbors = [self.seq_list[residue] for residue in self.neighbors]
        #return self.resnames_neighbors

    def get_charge_of_neighbors(self):
        self.distances_residue_mc_neighbors()
        
        self.charge_of_neighbors = []
        for residue in self.neighbors:
            charge_neighbor = next((charge for r, charge in self.charged_residues if r == residue), 0)
            self.charge_of_neighbors.append(charge_neighbor)
        return self.charge_of_neighbors

    def neighbors_acid_basic(self):
        self.acids = [-1 if self.type_dict[res] == 'A' else 0 for res in self.resnames_neighbors]
        self.basics = [1 if self.type_dict[res] == 'B' else 0 for res in self.resnames_neighbors]
        self.glicinas = [1 if res == 'GLY' else 0 for res in self.resnames_neighbors]
        return self.acids, self.basics, self.glicinas

    def count_polars_no_polars(self):
        self.get_resnames_neighbors()
        
        # Parametros Polares y No_polares, radio de la esfera para medir
        R_max = 5/10  # Umbral para distancia de polares
        R_max_no_polares = 7/10  # Umbral para distancia de no polares
        tau = 0.1*100  # Valor para tau en la fórmula
        
        self.polares = []
        self.no_polares = []
        for res, distancia in zip(self.resnames_neighbors, self.distances):
            if self.type_dict[res] in ['P', 'B', 'A']:
                valor_polar = 1 if distancia <= R_max else math.exp(-tau * (distancia - R_max)**2)
                self.polares.append(valor_polar)
            else:
                self.polares.append(0)
        
            if self.type_dict[res] == 'N':
                valor_no_polar = 1 if distancia <= R_max_no_polares else math.exp(-tau * (distancia - R_max_no_polares)**2)
                self.no_polares.append(valor_no_polar)
            else:
                self.no_polares.append(0)
    
        return self.polares, self.no_polares

    def sum_polares_no_polares(self):
        self.count_polars_no_polars()
        self.environment_polar = sum(self.polares)  # Sumar vecinos polares
        self.environment_no_polar = sum(self.no_polares)# Sumar vecinos no polares
        return self.environment_no_polar, self.environment_polar

    def calculate_term_polar(self): # todos los Bp Bnp estan multiplicados por Kb * T = 0.001987 * 300, por un tema de unidades
        self.sum_polares_no_polares()
        self.polar_no_polar_parameter ={'ASP':(0.72426,5.57075),'CTR':(0.72426,5.57075),'GLU':(0.380763,9.44846),'LYS':(0.0110311,5.65882),
                                       'ARG':(0.0110311,5.65882),'NTR':(0.0110311,5.65882),'HIS':(0.602896,9.36507),'CYS':(0.037 * 0.001987 * 300,85.91 * 0.001987 * 300),
                                       'TYR':(0.0, 0.0)}
        Bp, Bnp = self.polar_no_polar_parameter[self.resname_residue_mc]
    
        Npmax = 3.9
        Nnpmax = 17.78
        alfaP = 0.416683
        alfanp = 0.049
        
        # Calculate Up and Unp
        self.Up = math.exp(-alfaP * (self.environment_polar - Npmax)**2) if self.environment_polar <= Npmax else 1
        self.Unp = math.exp(-alfanp * (self.environment_no_polar - Nnpmax)**2) if self.environment_no_polar <= Nnpmax else 1
        
        # Calculate the term
        self.term_polar = self.acid_basic_residue_mc  * (Bnp * self.Unp - Bp * self.Up)
        return self.term_polar


    def calculate_term_elec(self):
        self.get_charge_of_neighbors()
        self.distances_residue_mc_neighbors()
        
        K_elec = 2.43232/10
        L = 10.0/10.0
        self.E_elect = []
        
        for charge_neighbor, distance in zip(self.charge_of_neighbors, self.distances):    
            if charge_neighbor == -1:  # Vecino ácido
                energy_electrostatic = self.charge * ((-1 / distance) * math.exp(-distance / L))
                self.E_elect.append(energy_electrostatic)
            elif charge_neighbor == 1:  # Vecino básico
                energy_electrostatic = self.charge * ((1 / distance) * math.exp(-distance / L))
                self.E_elect.append(energy_electrostatic)
            elif charge_neighbor == 0:  # Vecino sin carga
                self.E_elect.append(0)
        # Sumar la energia Electrostatica
        self.term_elec = sum(self.E_elect) * K_elec
        return self.term_elec

    def calculate_new_charge(self):
        self.choose_residue_mc()
        self.new_charge = charge_flip(self.acid_basic_residue_mc, self.charge)
        self.delta_carga = self.new_charge - self.charge
        return self.new_charge, self.delta_carga
    
    def calculate_delta_term_elec(self):
        self.calculate_new_charge()
        self.get_charge_of_neighbors()
        self.distances_residue_mc_neighbors()
        self.calculate_term_elec()
                
        K_elec = 2.43232/10
        L = 10.0/10.0
        self.new_E_elect = []
        
        for charge_neighbor, distance in zip(self.charge_of_neighbors, self.distances):    
            if charge_neighbor == -1:  # Vecino ácido
                energy_electrostatic = self.new_charge * ((-1 / distance) * math.exp(-distance / L))
                self.new_E_elect.append(energy_electrostatic)
            elif charge_neighbor == 1:  # Vecino básico
                energy_electrostatic = self.new_charge * ((1 / distance) * math.exp(-distance / L))
                self.new_E_elect.append(energy_electrostatic)
            elif charge_neighbor == 0:  # Vecino sin carga
                self.new_E_elect.append(0)
        self.new_term_elec = sum(self.new_E_elect) * K_elec
        self.delta_term_elec = self.new_term_elec - self.term_elec
        return self.delta_term_elec

    def calculate_delta_term_polar(self):
        self.calculate_new_charge()
        self.calculate_term_polar()
        
        self.delta_term_polar = self.delta_carga * self.term_polar
        return self.delta_term_polar

    def asignation_pKa_ref(self):
        self.acid_basic_residue_mc()
        # pKas libres en solución
        pKa_ref_dict = {
            'ASP': 4.0, 'GLU': 4.5, 'LYS': 10.6, 'ARG': 12.0,
            'HIS': 6.4, 'CYS': 8.3, 'TYR': 11.0, 'NTR': 7.5, 'CTR': 3.5
        }
        self.pKa_ref = pKa_ref_dict.get(self.resname_residue_mc)
        return self.pKa_ref
        
    def calculate_delta_term_pH(self):
        self.asignation_pKa_ref()
        kb = 0.001987
        T = 300
        self.delta_term_pH = self.delta_carga * (int(self.pH) - int(self.pKa_ref)) * kb * T * np.log(10)

    def accept_or_reject(self):
        self.calculate_delta_term_elec()
        self.calculate_delta_term_pH()
        self.calculate_delta_term_polar()
        
        kb = 0.001987
        T = 300
        self.lista_cambios = []
        self.change_energy_mc = self.delta_term_pH + self.delta_term_elec + self.delta_term_polar
        
        if self.change_energy_mc < 0:
        # Si la energía disminuye, aceptar la nueva carga
            for index, (residuo, _) in enumerate(self.charged_residues):
                if residuo == self.residue_mc:
                    self.charged_residues[index] = (residuo, self.new_charge)
                    self.lista_cambios.extend(self.charged_residues)  # Agregar elementos directamente
        else:
            # Si la energía aumenta, aceptar con una probabilidad determinada por el criterio de Metropolis
            random_prob = random.uniform(0, 1)
            if random_prob <= np.exp(-self.change_energy_mc/(kb*T)):
                for index, (residuo, _) in enumerate(self.charged_residues):
                    if residuo == self.residue_mc:
                        self.charged_residues[index] = (residuo, self.new_charge)
                        self.lista_cambios.extend(self.charged_residues)  # Agregar elementos directamente
            else:
                # Si no se acepta el cambio, agregar la configuración actual a la lista de cambios
                self.lista_cambios.extend(self.charged_residues)  # Agregar elementos directamente

        bond_matrix = []
        for _, row in self.enlaces_matrix.iterrows():
            if row['seq_i'] == self.residue_mc or row['seq_j'] == self.residue_mc:
                bond_index = int(row['bond_index'])  # Asegurar que bond_index es un entero
                # Obtener las cargas correspondientes de charged_residues
                carga_i = next(carga for residuo, carga in self.lista_cambios if residuo == row['seq_i'])
                carga_j = next(carga for residuo, carga in self.lista_cambios if residuo == row['seq_j'])
                bond_matrix.append([bond_index, carga_i, carga_j])
        
        # Convertir bond_matrix a DataFrame
        self.bond_matrix_df = pd.DataFrame(bond_matrix, columns=['bond_index', 'carga_i', 'carga_j'])
        return self.lista_cambios, self.bond_matrix_df
