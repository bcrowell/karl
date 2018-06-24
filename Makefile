PYTHON = python3
OBJ = obj
SRC = src
JS = js
BROWSER_PHYSICS = browser/physics

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
    $(JS)/io_util.js          $(JS)/test.js \
    $(JS)/lambert_w_stuff.js  $(JS)/test_lambert_w.js  \
    $(JS)/schwarzschild.js    $(JS)/test_schwarzschild.js   \
    $(JS)/kruskal.js          $(JS)/test_kruskal.js   \
    $(JS)/transform.js        $(JS)/test_transform.js \
    $(JS)/runge_kutta.js      $(JS)/test_runge_kutta.js   \
    $(JS)/angular.js          $(JS)/test_angular.js \
    $(JS)/math_util.js        $(JS)/test_math_util.js \
    $(JS)/vector.js           $(JS)/test_vector.js \
                              $(JS)/test_math.js

BROWSER_PHYSICS_FILES = \
    $(BROWSER_PHYSICS)/io_util.js        \
    $(BROWSER_PHYSICS)/lambert_w_stuff.js \
    $(BROWSER_PHYSICS)/schwarzschild.js  \
    $(BROWSER_PHYSICS)/kruskal.js        \
    $(BROWSER_PHYSICS)/transform.js      \
    $(BROWSER_PHYSICS)/runge_kutta.js    \
    $(BROWSER_PHYSICS)/angular.js        \
    $(BROWSER_PHYSICS)/math_util.js      \
    $(BROWSER_PHYSICS)/vector.js         \

TESTS = math math_util lambert_w angular schwarzschild kruskal transform vector runge_kutta

VPATH = src

test: $(PY)
	@for test in $(TESTS); do \
	  $(PYTHON3) $(OBJ)/test_$${test}.py; \
	done
	@echo "To make this a real test, do: rm obj/*.py ; make test . Otherwise is doesn't test for regressions in PJ."

$(PY): $(OBJ)/%.py: $(SRC)/%.pp $(SRC)/*.h
	filepp -DLANG=python $< -o $@
	@chmod +x $@

$(BROWSER_PHYSICS_FILES): $(BROWSER_PHYSICS)/%.js: $(JS)/%.js
	@mkdir -p $(BROWSER_PHYSICS)
	cp $< $@

js: $(JS_FILES) $(BROWSER_PHYSICS_FILES)
	@#

clean_js:
	rm -f js/*.js js/*.jsi

test_js:
	@for test in $(TESTS); do \
	  cd js ; rhino -opt -1 test_$${test}.js ; cd -; \
	done

$(JS)/%.js: $(SRC)/%.pp $(SRC)/*.h
	@filepp -DLANG=js $< -o $@i
	pj/pj.rb $@i karl <$@i >$@
	@rm $@i
	@-js-beautify --replace -n -s 2 $@

clean:
	rm -f *~ src/*~ obj/*~ pj/*~ js/*~
	cd doc && make clean && cd -
