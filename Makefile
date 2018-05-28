test:
	./karl.py

clean:
	rm -f *~
	cd doc && make clean && cd -
