PYTHON = python3

SOURCES = io_util.py test.py \
    lambert_w.py  test_lambert_w.py  \
    schwarzschild.py  test_schwarzschild.py   \
    kruskal.py  test_kruskal.py   \
    runge_kutta.py  test_runge_kutta.py   \
    angular.py test_angular.py \
    math_util.py test_math_util.py \
    util.h math.h init.h 

TESTS = math_util lambert_w angular schwarzschild runge_kutta kruskal 


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
