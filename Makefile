PYTHON = python3
OBJ = obj
SRC = src
JS = js

PY = \
    $(OBJ)/io_util.py        $(OBJ)/test.py \
    $(OBJ)/lambert_w_stuff.py $(OBJ)/test_lambert_w.py  \
    $(OBJ)/schwarzschild.py  $(OBJ)/test_schwarzschild.py   \
    $(OBJ)/kruskal.py        $(OBJ)/test_kruskal.py   \
    $(OBJ)/transform.py      $(OBJ)/test_transform.py \
    $(OBJ)/runge_kutta.py    $(OBJ)/test_runge_kutta.py   \
    $(OBJ)/angular.py        $(OBJ)/test_angular.py \
    $(OBJ)/math_util.py      $(OBJ)/test_math_util.py \
    $(OBJ)/vector.py         $(OBJ)/test_vector.py \
    $(OBJ)/test_math.py

JS_FILES = \
    $(JS)/io_util.js        $(JS)/test.js \
    $(JS)/lambert_w_stuff.js $(JS)/test_lambert_w.js  \
    $(JS)/schwarzschild.js  $(JS)/test_schwarzschild.js   \
    $(JS)/kruskal.js        $(JS)/test_kruskal.js   \
    $(JS)/transform.js      $(JS)/test_transform.js \
    $(JS)/runge_kutta.js    $(JS)/test_runge_kutta.js   \
    $(JS)/angular.js        $(JS)/test_angular.js \
    $(JS)/math_util.js      $(JS)/test_math_util.js \
    $(JS)/vector.js         $(JS)/test_vector.js \
    $(JS)/test_math.js

TESTS = math math_util lambert_w angular schwarzschild runge_kutta kruskal transform vector

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

js: $(JS_FILES)
	@#

clean_js:
	rm -f js/*.js js/*.jsi

test_js:
	cd js ; rhino -opt -1 test_math.js ; cd -
	cd js ; rhino -opt -1 test_math_util.js ; cd -
	cd js ; rhino -opt -1 test_lambert_w.js ; cd -
	cd js ; rhino -opt -1 test_schwarzschild.js ; cd -
	cd js ; rhino -opt -1 test_kruskal.js ; cd -

$(JS)/%.js: $(SRC)/%.pp
	@filepp -DLANG=js $< -o $@i
	pj/pj.rb $@i <$@i >$@
	@-js-beautify --replace -n $@

clean:
	rm -f *~ src/*~ obj/*~ pj/*~ js/*~
	cd doc && make clean && cd -
