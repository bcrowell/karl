star catalog
============

This script builds a sqlite database of stars, used by the optics app.

Requires ruby's sqlite bindings: apt-get install ruby-sqlite3 .
To run this script, do a
sudo make /usr/share/karl/star_catalog.sqlite .

The script build_star_catalog.rb reads the file /usr/share/kstars/stars.dat, which
can be obtained on a Debian system by doing an apt-get install kstars-data.
It writes the file star_catalog.sqlite. 

