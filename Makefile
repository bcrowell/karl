test:
	./karl.py

doc:
	pdflatex doc

invariants:
	maxima -b invariants.mac

christoffel:
	maxima -b christoffel.mac
	# Now cut and paste output into clean_up_christoffel.rb and run it.

