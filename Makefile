MAKEFLAGS += --no-print-directory

PYTHON = python3

SRC = src
OBJ = obj
JS = js

FILEPP_OPTIONS = -Isrc/include

VPATH = src

.PHONY: clean clean_js js js_all py c all doc

all:
	make depend
	make py
	make js_all
	make c
	make test

depend: gen_depends.py
	@python3 gen_depends.py >depend

js_all:
	make js
	mkdir -p js/lib
	cp src/lib/*.js js/lib
	cp src/lib/*.js browser/lib # except  ...
	rm -f browser/lib/loader.js # ... loader.js
	cp js/*.js browser/physics # except ...
	rm -f browser/physics/test* # ... tests
	cp src/browser/*.js browser
	cp src/browser/*.html browser
	cp src/browser/util/*.js browser/util

clean_js:
	rm -f js/*.js js/*.jsi

clean_py:
	rm -f obj/*.py

clean:
	rm -f *~ src/*/*~ obj/*~ pj/*~ js/*~ js/*.jsi obj/*.pyc
	cd doc && make clean && cd -

doc: doc/doc.pdf
	#

doc/doc.pdf: doc/doc.tex
	cd doc ; make ; cd -

include depend
include c.mk
# ... rules for C code

