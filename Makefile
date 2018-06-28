PYTHON = python3

SRC = src
OBJ = obj
JS = js

include depend
include c.mk
# ... rules for C code

VPATH = src

.PHONY: clean clean_js js py c all

depend: gen_depends.py
	python3 gen_depends.py >depend

all:
	make depend
	make py
	make js
	make c

clean_js:
	rm -f js/*.js js/*.jsi

clean:
	rm -f *~ src/*~ obj/*~ pj/*~ js/*~ js/*.jsi obj/*.pyc
	cd doc && make clean && cd -
