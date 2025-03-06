# Ejemplo de simulación a pH constante con openAWSEM


Para poder correr el ejemplo es necesario instalar openAWSEM. En la página openawsem.org se encuentra detallado todo el proceso y comandos para su instalación.

## Correr el ejemplo

- Copiar la carpeta 1zug_awsem_curva_T_pH7 a tu escritorio.
- Entrar a la carpeta usando la terminal de linux: 
```
cd /home/user/Escritorio/1zug_awsem_curva_T_pH7/
```
- Una vez dentro, activar el entorno virtual de openAWSEM
- correr el siguiente comando para inciar la simulación: 
```
./mm_run.py 1zug --platform CPU --steps 5000 --tempStart 700 --tempEnd 900 -m 1 -f forces_setup.py --r 100
```

Lo que estamos corriendo es una simulacion corta de 5000 pasos de 700 a 900 K de temperatura. "-m 1" = se usa para aclarar que es una curva de temperatura, tambien se puede correr con "-m 0" para correr a temperatura constante, pero el pH no esta incorporado ahi correctamente. 


--r 100 es el reportero de la simulación, la informacion se guardara cada 100 pasos

1zug es el nombre del archivo PDB 

# Archivos que se agregaron

- crystal_structure.fasta.terminals  (la secuencia de la 1zug agregando terminales como "X" n-terminal, "Z" c-terminal) 
- charge_on_residues.dat (es igual al charge.txt pero tiene agregado la carga de los terminales)  
- debyeHuckelTerms_with_terminals.py (Este archivo se usa para agregar las cargas a los atomos en la simulación)  
- funciones_procesamiento.py (funciones complementarias para antes de correr la simulación, se usan en el mm_run.py) 
- Montecarlo_pH.py (Aca se encuentra el modelo de pH que se usa para decidir el cambio de cargas de cada aminoacido para simular el pH)

# Archivos que se modificaron
- mm_run.py (Este scrip es el que corre la simulacion, dentro de sub mode 1 se hicieron todas las modificaciones)
- forces_setup.py (aca modifique un poco el forces orginal para que lea bien el archivo de debyeHuckelTerms_with_terminals.py y no el que tiene por defecto)
