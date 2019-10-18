CC=gcc
FILE_SET_1=simulation_clib.c pcg_basic.c 
OBJ_1=simulation_clib.so

#FILE_SET_2=analyzer_clib.c pcg_basic.c 
#OBJ_2=analyzer_clib.so

FILE_SET_3=ensembleSimulator.c pcg_basic.c 
OBJ_3=ensembleSimulator.so

FILE_SET_4=pertParamSim.c pcg_basic.c 
OBJ_4=pertParamSim.so




build: $(FILE_SET_1)  $(FILE_SET_2)
	$(CC) -fPIC -shared -o $(OBJ_1) $(FILE_SET_1) 
	#$(CC) -fPIC -shared -o $(OBJ_2) $(FILE_SET_2) 
	$(CC) -fPIC -shared -o $(OBJ_3) $(FILE_SET_3) 
	$(CC) -fPIC -shared -o $(OBJ_4) $(FILE_SET_4) 

