Now, you are in folder AssignmentMPI, you have the following options: 

-Run on Astral with default parameters*: execute scriptAutoLaptop.sh (no need to load modules)
-Run on your laptop with default parameters*: execute scriptAutoLaptop.sh (no need to load modules)

-Modify the number of processor to run on Astral : modify the last line of bashScript<Scheme> (and do not forget to modify the name of the job)
-Modify the number of processor to run on your laptop: modify the number after mpirun on the scriptAutoLaptop.sh

-Modify the number of points: At the begining of each scheme file (in ParallelC and in the main in SerieCPP), you simply modify the constant variable NB_POINTS.

If you want to check the code, you will find all parallel codes in the ParallelC folder and all the serial code in the SerieCPP folder.
If you want to check the outputs of the simulations, everything will be written in the appropriate CSV file.

At the end of a simulation, if you want to clean all useless files ".o", ".e", CSV files etc... Execute clean.sh


*default parameters: 1000 points and the optimal number of processors for each scheme, see the report.
