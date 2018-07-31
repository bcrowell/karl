MAKEFLAGS += --no-print-directory

PYTHON = python3

SRC = src
OBJ = obj
JS = js

FILEPP_OPTIONS = -Isrc/include

VPATH = src

KARL_DATA = /usr/share/karl
STAR_CATALOG = $(KARL_DATA)/star_catalog.sqlite

.PHONY: clean clean_js js js_all py c all doc

all:
	make depend
	make py
	make js_all
	make c
	make test

optics: obj/optics.py
	make py
	@chmod +x obj/optics.py
	obj/optics.py && display stars.png

animation:
	convert -delay 20 -loop 0 anim*.png animation.gif
	#display animation.gif
	# ... doesn't work, instead do gimp : filters : animation : playback

depend: gen_depends.py
	@python3 gen_depends.py >depend

js_all:
	make js
	mkdir -p js/lib
	mkdir -p browser/physics
	mkdir -p browser/lib
	mkdir -p browser/util
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

uninstall:
	rm -f $(STAR_CATALOG)
	rmdir $(KARL_DATA)

doc: doc/doc.pdf
	#

doc/doc.pdf: doc/doc.tex
	cd doc ; make ; cd -

include depend
include c.mk
# ... rules for C code

$(STAR_CATALOG): /usr/share/kstars/stars.dat
	ruby star_catalog/build_star_catalog.rb
	mkdir -p $(KARL_DATA)
	mv mag7.sqlite $(STAR_CATALOG)


