test:
	./karl.py

doc:
	pdflatex doc

invariants:
	maxima -b invariants.mac

christoffel:
	maxima -b christoffel.mac
	# Now cut and paste output into christoffel.rb and run it.

