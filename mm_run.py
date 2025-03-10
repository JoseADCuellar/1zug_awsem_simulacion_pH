#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
from time import sleep
import fileinput
import importlib.util

from openawsem import *
from openawsem.helperFunctions.myFunctions import *
###change###
from openmm.app import PDBFile, PDBReporter
from openmm import Vec3
import json
from Montecarlo_pH import *
from funciones_procesamiento import *
###change###

# simulation_platform = "CPU"  # OpenCL, CUDA, CPU, or Reference
# simulation_platform = "OpenCL"
do = os.system
cd = os.chdir


def run(args):
    simulation_platform = args.platform
    platform = Platform.getPlatformByName(simulation_platform)
    if simulation_platform == "CPU":
        if args.thread != -1:
            platform.setPropertyDefaultValue("Threads", str(args.thread))
        print(f"{simulation_platform}: {platform.getPropertyDefaultValue('Threads')} threads")

    # if mm_run.py is not at the same location of your setup folder.
    setupFolderPath = os.path.dirname(args.protein)
    setupFolderPath = "." if setupFolderPath == "" else setupFolderPath
    proteinName = pdb_id = os.path.basename(args.protein)


    pwd = os.getcwd()
    toPath = os.path.abspath(args.to)
    checkPointPath = None if args.fromCheckPoint is None else os.path.abspath(args.fromCheckPoint)
    forceSetupFile = None if args.forces is None else os.path.abspath(args.forces)
    parametersLocation = "." if args.parameters is None else os.path.abspath(args.parameters)
    os.chdir(setupFolderPath)



    # chain=args.chain.upper()
    chain=args.chain
    pdb = f"{pdb_id}.pdb"

    if chain == "-1":
        chain = getAllChains("crystal_structure.pdb")
        print("Chains to simulate: ", chain)


    if args.to != "./":
        # os.system(f"mkdir -p {args.to}")
        os.makedirs(toPath, exist_ok=True)
        os.system(f"cp {forceSetupFile} {toPath}/forces_setup.py")
        os.system(f"cp crystal_structure.fasta {toPath}/")
        os.system(f"cp crystal_structure.pdb {toPath}/")
        # os.system(f"cp {pdb} {args.to}/{pdb}")
        # pdb = os.path.join(args.to, pdb)

    if args.fromOpenMMPDB:
        input_pdb_filename = proteinName
        seq=read_fasta("crystal_structure.fasta")
        print(f"Using Seq:\n{seq}")
    else:
        suffix = '-openmmawsem.pdb'
        if pdb_id[-len(suffix):] == suffix:
            input_pdb_filename = pdb_id
        else:
            input_pdb_filename = f"{pdb_id}-openmmawsem.pdb"
        seq=None

    if args.fasta == "":
        seq = None
    else:
        seq = seq=read_fasta(args.fasta)
        print(f"Using Seq:\n{seq}")
    # start simulation
    collision_rate = 5.0 / picoseconds
    checkpoint_file = "checkpnt.chk"
    checkpoint_reporter_frequency = 10000



    snapShotCount = 400
    stepsPerT = int(args.steps/snapShotCount)
    Tstart = args.tempStart
    Tend = args.tempEnd
    if args.reportFrequency == -1:
        if stepsPerT == 0:
            reporter_frequency = 4000
        else:
            reporter_frequency = stepsPerT
    else:
        reporter_frequency = args.reportFrequency
    # reporter_frequency = 4000

    print(f"using force setup file from {forceSetupFile}")
    spec = importlib.util.spec_from_file_location("forces", forceSetupFile)
    # print(spec)
    forces = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(forces)


    oa = OpenMMAWSEMSystem(input_pdb_filename, k_awsem=1.0, chains=chain, xml_filename=openawsem.xml, seqFromPdb=seq, includeLigands=args.includeLigands)  # k_awsem is an overall scaling factor that will affect the relevant temperature scales
    myForces = forces.set_up_forces(oa, submode=args.subMode, contactParameterLocation=parametersLocation)
    # print(forces)
    # oa.addForces(myForces)
    oa.addForcesWithDefaultForceGroup(myForces)

    if args.fromCheckPoint:
        integrator = LangevinIntegrator(Tstart*kelvin, 1/picosecond, args.timeStep*femtoseconds)
        simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
        simulation.loadCheckpoint(checkPointPath)
    else:
        # output the native and the structure after minimization
        integrator = CustomIntegrator(0.001)
        simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
        simulation.context.setPositions(oa.pdb.positions)  # set the initial positions of the atoms
        simulation.reporters.append(PDBReporter(os.path.join(toPath, "native.pdb"), 1))
        simulation.reporters.append(DCDReporter(os.path.join(toPath, "movie.dcd"), 1))
        simulation.step(int(1))
        simulation.minimizeEnergy()  # first, minimize the energy to a local minimum to reduce any large forces that might be present
        simulation.step(int(1))


        # print("------------------Folding-------------------")
        # oa = OpenMMAWSEMSystem(input_pdb_filename, k_awsem=1.0, chains=chain, xml_filename=OPENAWSEM_LOCATION+"awsem.xml")  # k_awsem is an overall scaling factor that will affect the relevant temperature scales
        # myForces = forces.set_up_forces(oa, submode=args.subMode, contactParameterLocation=parametersLocation)
        # oa.addForces(myForces)

        integrator = LangevinIntegrator(Tstart*kelvin, 1/picosecond, args.timeStep*femtoseconds)
        # integrator.setRandomNumberSeed(A_NUMBER_AS_RANDOM_SEED)
        # integrator = CustomIntegrator(0.001)
        simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
        # simulation.loadState(os.path.join(toPath, 'output.xml'))
        simulation.context.setPositions(oa.pdb.positions)  # set the initial positions of the atoms
        simulation.context.setVelocitiesToTemperature(Tstart*kelvin)  # set the initial velocities of the atoms according to the desired starting temperature
        # simulation.context.setVelocitiesToTemperature(Tstart*kelvin, A_RANDOM_SEED_NUMBER)
        #simulation.minimizeEnergy()  # first, minimize the energy to a local minimum to reduce any large forces that might be present


    print("reporter_frequency", reporter_frequency)
    simulation.reporters.append(StateDataReporter(stdout, reporter_frequency, step=True, potentialEnergy=True, temperature=True))  # output energy and temperature during simulation
    simulation.reporters.append(StateDataReporter(os.path.join(toPath, "output.log"), reporter_frequency, step=True, potentialEnergy=True, temperature=True)) # output energy and temperature to a file
    simulation.reporters.append(PDBReporter(os.path.join(toPath, "movie.pdb"), reportInterval=reporter_frequency))  # output PDBs of simulated structures
    simulation.reporters.append(DCDReporter(os.path.join(toPath, "movie.dcd"), reportInterval=reporter_frequency, append=True))  # output PDBs of simulated structures
    # simulation.reporters.append(DCDReporter(os.path.join(args.to, "movie.dcd"), 1))  # output PDBs of simulated structures
    # simulation.reporters.append(PDBReporter(os.path.join(args.to, "movie.pdb"), 1))  # output PDBs of simulated structures
    #####Change#####
    #checkpoint_file = os.path.join(toPath, "checkpnt.chk")
    #simulation.reporters.append(CheckpointReporter(checkpoint_file, checkpoint_reporter_frequency))
    #####Change#####
    #simulation.reporters.append(CheckpointReporter(os.path.join(toPath, checkpoint_file), checkpoint_reporter_frequency))  # save progress during the simulation
    print("Simulation Starts")
    start_time = time.time()
#####Change####       
    interrupt_frequency = args.interruptFrequency
    cb_atoms_matrix = id_and_residues_df_oasistem(oa)#crea un df con los id y numero de residuo que me da el obejto oa (sistema.topology)
    charged_residues = procesador_de_archivo_con_residuos_cargados('charge_on_residues.dat')#lee el archivo y los pasa a una lista de tuplas (Residuo,carga)
    num_bonds = myForces[-1].getNumBonds()
    #print(num_bonds)
    last_force = myForces[-1]
    atom_list = list(oa.pdb.topology.atoms())
    seq = read_fasta("crystal_structure.fasta.terminals")
    print(seq)
    Hawsem_state=[]
    # Antes de empezar la simulacion, se abre el archivo en modo 'write' para crear/limpiar el archivo output
    with open('Hawsem.state', 'w') as f:
        pass

    if args.simulation_mode == 0:
        for step in range(0, int(args.steps), interrupt_frequency):
            remaining_steps = min(interrupt_frequency, int(args.steps) - step)
            simulation.step(remaining_steps)
        
        # Obtener las posiciones actuales de los átomos
            state = simulation.context.getState(getPositions=True)
            positions = state.getPositions()

    
            cb_positions_df = constructor_df_de_posiciones_de_CB(cb_atoms_matrix, positions) #loop for para iterar sobre cb_atom_matrix y obtener las posiciones de cada CB. Me devuelve un df ['Residue_Number', 'Atom_Index', 'X', 'Y', 'Z']
                
            enlaces_matrix = constructor_df_de_los_objetos_force(last_force, num_bonds, atom_list, oa) #accede al objeto force de Debye Huckel y a sus parametros para armar un df ['bond_index', 'seq_i', 'seq_j', 'carga_i', 'carga_j'] usa un loop for num_bonds
                
            pH = 3
                
            charged_residues, bond_matrix_df = Montecarlo(pH, charged_residues, cb_positions_df, seq, enlaces_matrix)

            for _, row in bond_matrix_df.iterrows():
                bond_index = int(row['bond_index'])
                new_charge_i = row['carga_i']
                new_charge_j = row['carga_j']
                bond_parameters = last_force.getBondParameters(bond_index)
                atom_index_i, atom_index_j, _ = bond_parameters
                last_force.setBondParameters(bond_index, atom_index_i, atom_index_j, (new_charge_i, new_charge_j))

            last_force.updateParametersInContext(simulation.context)

            # Guardar enlaces_matrix en un archivo CSV
            #enlaces_matrix.to_csv('enlaces_matrix.csv', index=False)
            
            # Guardar el estado de charged_residues en Hawsem.state
            if total_steps % reporter_frequency == 0:
                Hawsem_state.append(json.dumps(charged_residues) + str(pH))#+ str(step)) puedo borrar el if y guardar cada cambio pero es mejor asi, para que sea igual al reporter_frecuency que lo defino en los argumentos
                #with open('Hawsem.state', 'a') as f:
                #    f.write(json.dumps(charged_residues) + str(pH) + "\n")
            
        with open('Hawsem.state', 'w') as f:
            for state in Hawsem_state:
                f.write(state + "\n")

####Change####
    elif args.simulation_mode == 1:
        # Inicializar variables necesarias
        pH=7
        interrupt_frequency = 250
        snapShotCount = 500
        stepsPerT = int(args.steps/snapShotCount)
        deltaT = (Tend - Tstart) / snapShotCount
        print("stepsperT", stepsPerT)
        lista_steps =generar_lista(int(interrupt_frequency),int(stepsPerT))
        #print(lista_steps)
        total_steps = 0  # Inicializar contador de steps totales
        Hawsem_state = []  # Lista para almacenar los estados de charged_residues
        for i in range(snapShotCount):
            integrator.setTemperature((Tstart + deltaT * i) * kelvin)
            for steps in lista_steps:
                if steps == int(interrupt_frequency):
                    simulation.step(steps)
                    total_steps += steps #por ahora solo usa esta forma, seria mejor no usar un contador, pero siempre que de sin resto genial
                    state = simulation.context.getState(getPositions=True)
                    positions = state.getPositions()
                    cb_positions_df = constructor_df_de_posiciones_de_CB(cb_atoms_matrix, positions)
                    enlaces_matrix = constructor_df_de_los_objetos_force(last_force, num_bonds, atom_list, oa)
                    # Llamar a Montecarlo
                    mc_pH = Montecarlo_change_charge(cb_positions_df, charged_residues,  seq, enlaces_matrix,pH)
                    charged_residues, bond_matrix_df = mc_pH.accept_or_reject()
            # Actualizar parámetros de enlace según los resultados de Monte Carlo
                    for _, row in bond_matrix_df.iterrows():
                        bond_index = int(row['bond_index'])
                        new_charge_i = row['carga_i']
                        new_charge_j = row['carga_j']
                        bond_parameters = last_force.getBondParameters(bond_index)
                        atom_index_i, atom_index_j, _ = bond_parameters
                        last_force.setBondParameters(bond_index, atom_index_i, atom_index_j, (new_charge_i, new_charge_j))

                    last_force.updateParametersInContext(simulation.context)
                #else:
                #    continue
            # Guardar el estado de charged_residues en Hawsem.state
                if total_steps % reporter_frequency == 0:
                    Hawsem_state.append(json.dumps(charged_residues) + str(pH))

    # Guardar todos los estados al final
        with open('Hawsem.state', 'w') as f:
            for state in Hawsem_state:
                f.write(state + "\n")
    print("Simulation completed successfully for mode 1.")
    
            # simulation.saveCheckpoint('step_%d.chk' % i)
            # simulation.context.setParameter("k_membrane", 0)
            # if i < snapShotCount/2:
            #     simulation.context.setParameter("k_membrane", (i % 2) * k_mem)
            #     simulation.context.setParameter("k_single_helix_orientation_bias", (i % 2) * k_single_helix_orientation_bias)
            # else:
            #     simulation.context.setParameter("k_membrane", k_mem)
            #     simulation.context.setParameter("k_single_helix_orientation_bias", k_single_helix_orientation_bias)

            # simulation.context.setParameter("k_membrane", (i)*(k_mem/snapShotCount))
            # simulation.context.setParameter("k_single_helix_orientation_bias", (i)*(k_single_helix_orientation_bias/snapShotCount))
            # print(simulation.context.getParameter("k_membrane"))


    # simulation.step(int(1e6))

    time_taken = time.time() - start_time  # time_taken is in seconds
    hours, rest = divmod(time_taken,3600)
    minutes, seconds = divmod(rest, 60)
    print(f"---{hours} hours {minutes} minutes {seconds} seconds ---")

    timeFile = os.path.join(toPath, "time.dat")
    with open(timeFile, "w") as out:
        out.write(str(time_taken)+"\n")

    # accompany with analysis run
    simulation = None
    time.sleep(10)
    os.chdir(pwd)
    print(os.getcwd())
    if args.fasta == "":
        analysis_fasta = ""
    else:
        analysis_fasta = f"--fasta {args.fasta}"
    if args.includeLigands:
        additional_cmd = "--includeLigands"
    else:
        additional_cmd = ""
    os.system(f"{sys.executable} mm_analyze.py {args.protein} -t {os.path.join(toPath, 'movie.dcd')} --subMode {args.subMode} -f {args.forces} {analysis_fasta} {additional_cmd} -c {chain}")

def main():
    # from run_parameter import *
    parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatic copy the template file, \
        run simulations")

    parser.add_argument("protein", help="The name of the protein")
    parser.add_argument("--name", default="simulation", help="Name of the simulation")
    parser.add_argument("--to", default="./", help="location of movie file")
    parser.add_argument("-c", "--chain", type=str, default="-1")
    parser.add_argument("-t", "--thread", type=int, default=-1, help="default is using all that is available")
    parser.add_argument("-p", "--platform", type=str, default="OpenCL")
    parser.add_argument("-s", "--steps", type=float, default=2e4, help="step size, default 1e5")
    parser.add_argument("--tempStart", type=float, default=800, help="Starting temperature")
    parser.add_argument("--tempEnd", type=float, default=200, help="Ending temperature")
    parser.add_argument("--fromCheckPoint", type=str, default=None, help="The checkpoint file you want to start from")
    parser.add_argument("-m", "--simulation_mode", type=int, default=1,
                    help="default 1,\
                            0: constant temperature,\
                            1: temperature annealing")
    parser.add_argument("--subMode", type=int, default=-1)
    parser.add_argument("-f", "--forces", default="forces_setup.py")
    parser.add_argument("--parameters", default=None)
    parser.add_argument("-r", "--reportFrequency", type=int, default=-1, help="default value step/400")
    parser.add_argument("--fromOpenMMPDB", action="store_true", default=False)
    parser.add_argument("--fasta", type=str, default="crystal_structure.fasta")
    parser.add_argument("--timeStep", type=float, default=2)
    parser.add_argument("--includeLigands", action="store_true", default=False)
    parser.add_argument("--interruptFrequency", type=int, default=5, help="Frequency of interruptions during simulation")
    args = parser.parse_args()


    with open('commandline_args.txt', 'a') as f:
        f.write(' '.join(sys.argv))
        f.write('\n')
    print(' '.join(sys.argv))

    run(args)

if __name__=="__main__":
    main()
