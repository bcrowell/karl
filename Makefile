PYTHON = python3
OBJ = obj
SRC = src

PY = \
    $(OBJ)/io_util.py        $(OBJ)/test.py \
    $(OBJ)/lambert_w.py      $(OBJ)/test_lambert_w.py  \
    $(OBJ)/schwarzschild.py  $(OBJ)/test_schwarzschild.py   \
    $(OBJ)/kruskal.py        $(OBJ)/test_kruskal.py   \
    $(OBJ)/transform.py      $(OBJ)/test_transform.py \
    $(OBJ)/runge_kutta.py    $(OBJ)/test_runge_kutta.py   \
    $(OBJ)/angular.py        $(OBJ)/test_angular.py \
    $(OBJ)/math_util.py      $(OBJ)/test_math_util.py \
    $(OBJ)/vector.py         $(OBJ)/test_vector.py

TESTS = math_util lambert_w angular schwarzschild runge_kutta kruskal transform vector

VPATH = src

test: $(PY)
	@for test in $(TESTS); do \
	  $(PYTHON3) $(OBJ)/test_$${test}.py; \
	done

test_kruskal: $(PY)
	$(PYTHON3) $(OBJ)/test_kruskal.py

test_runge_kutta: $(PY)
	$(PYTHON3) $(OBJ)/test_runge_kutta.py

$(PY): $(OBJ)/%.py: $(SRC)/%.pp
	filepp -DLANG=python $< -o $@
	@chmod +x $@

clean:
	rm -f *~
	cd doc && make clean && cd -
