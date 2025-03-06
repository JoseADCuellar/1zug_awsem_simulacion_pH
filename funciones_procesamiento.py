
import numpy as np
import pandas as pd
import random
import math

#############################################
##### MAS FUNCIONES #########################
#############################################

def procesador_de_archivo_con_residuos_cargados(archivo):
    with open(archivo, 'r') as file:
        residues = []
        for line in file:
            parts = line.split()
            residue = int(parts[0])
            charge = float(parts[1])
            residues.append((residue, charge))

    # Filtrar los residuos con carga diferente de cero
    charged_residues = [(residue, charge) for residue, charge in residues if charge != 0.0]
    
    return charged_residues


def constructor_df_de_posiciones_de_CB(cb_atoms_matrix, positions): # Funciona con los terminales
    cb_positions = []

    # Iterar sobre los índices de átomos CB y números de residuo en cb_atoms_matrix
    for atom_index, residue_number in cb_atoms_matrix:
        cb_position = positions[atom_index]
        cb_positions.append([residue_number-1, atom_index, cb_position.x, cb_position.y, cb_position.z])
    
    # Definir las columnas del DataFrame
    columns = ['Residue_Number', 'Atom_Index', 'X', 'Y', 'Z']
    cb_positions_df = pd.DataFrame(cb_positions, columns=columns)
    # Ordenar por la columna Residue_Number
    cb_positions_df = cb_positions_df.sort_values(by='Residue_Number').reset_index(drop=True)
    
    return cb_positions_df
    
def constructor_df_de_los_objetos_force(last_force, num_bonds, atom_list, oa):
    enlaces_info = []

    for bond_index in range(num_bonds):
        bond_parameters = last_force.getBondParameters(bond_index)
        atom_index_i = bond_parameters[0]
        atom_index_j = bond_parameters[1]
        charge_i, charge_j = bond_parameters[2]
        residue_index_i = atom_list[atom_index_i].residue.index
        residue_index_j = atom_list[atom_index_j].residue.index
        num_seq_i = oa.residues[residue_index_i].index
        num_seq_j = oa.residues[residue_index_j].index
        atom_name_i = atom_list[atom_index_i].name
        atom_name_j = atom_list[atom_index_j].name
        #### poner los terminales ficticios
        if atom_name_i == 'O':
            num_seq_i = oa.nres + 1
        if atom_name_j == 'CA' and num_seq_j == 0: # si es ionizalbe y le falta el CB puede fallar
            num_seq_j = oa.nres
        if atom_name_j == 'O':
            num_seq_j = oa.nres + 1
        if atom_name_i == 'CA'and num_seq_i == 0: # si es ionizalbe y le falta el CB puede fallar
            num_seq_i = oa.nres
        ####
        enlace_info = [bond_index, num_seq_i, num_seq_j,atom_name_i,atom_name_j, charge_i, charge_j]
        enlaces_info.append(enlace_info)

    columns = ['bond_index', 'seq_i', 'seq_j','atom_i', 'atom_j', 'carga_i', 'carga_j']
    enlaces_matrix = pd.DataFrame(enlaces_info, columns=columns)
    
    return enlaces_matrix

def id_and_residues_df_oasistem(oa):
    cb_atoms_matrix = []
    residues = list(oa.pdb.topology.residues())  # Convertir los residuos a una lista

    for idx, residue in enumerate(residues):
        has_cb = False  # Bandera para verificar si el residuo tiene un átomo CB

        for atom in residue.atoms():
            # Caso general: Carbono Beta (CB)
            if atom.name == 'CB':
                cb_atoms_matrix.append((atom.index, int(residue.id)))
                has_cb = True

            # Caso N-terminal: Carbono Alfa (CA)
            if idx == 0 and atom.name == 'CA':  # Residuo N-terminal
                cb_atoms_matrix.append((atom.index, int(oa.nres + 1)))

            # Caso C-terminal: Oxígeno (O)
            if idx == len(residues) - 1 and atom.name == 'O':  # Residuo C-terminal
                cb_atoms_matrix.append((atom.index, int(oa.nres + 2)))

        # Si el residuo no tiene CB, usar CA
        if not has_cb:
            for atom in residue.atoms():
                if atom.name == 'CA':
                    cb_atoms_matrix.append((atom.index, int(residue.id)))
                    break  # No necesitamos seguir buscando después de encontrar CA

    return cb_atoms_matrix





