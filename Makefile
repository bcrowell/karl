PYTHON = python3

SOURCES = io_util.py test.py \
    lambert_w.py  test_lambert_w.py  \
    schwarzschild.py  test_schwarzschild.py   \
    runge_kutta.py  test_runge_kutta.py   \
    util.h math.h init.h

TESTS = lambert_w schwarzschild runge_kutta


%.py: %.pp
	filepp $< -o $@
	@chmod +x $@

test: $(SOURCES)
	@for test in $(TESTS); do \
	  $(PYTHON3) test_$${test}.py; \
	done

clean:
	rm -f *~
	cd doc && make clean && cd -
