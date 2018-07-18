# make rules for C code

c: $(OBJ)/karl.so

$(OBJ)/%.o: $(SRC)/physics/%.c $(SRC)/include/*.h
	gcc -I$(SRC)/include -c -Wall -Werror -fpic $< -o $@

$(OBJ)/lambert_w.o: $(SRC)/physics/lambert_w.cpp
	gcc -c -Wall -Werror -fpic $< -o $@

veberic_lambert_w/LambertW.o: veberic_lambert_w/LambertW.cc
	cd veberic_lambert_w ; make ; cd -

$(OBJ)/karl.so: $(OBJ)/apply_christoffel.o $(OBJ)/lambert_w.o veberic_lambert_w/LambertW.o
	gcc -shared -o $(OBJ)/karl.so $(OBJ)/apply_christoffel.o $(OBJ)/lambert_w.o \
              veberic_lambert_w/LambertW.o -lm
