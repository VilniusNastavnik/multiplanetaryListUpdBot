# multiplanetaryListUpdBot
Python program to update https://it.wikipedia.org/wiki/Sistemi_multiplanetari with data from http://exoplanet.eu/, https://exoplanetarchive.ipac.caltech.edu/ and https://simbad.u-strasbg.fr/simbad/.

The program creates a file 'tabella.wiki' with the template to put in https://it.wikipedia.org/wiki/Sistemi_multiplanetari.

Before publishing the changes check them carefully!

It uses astroquery.simbad (python -m pip install -U --pre astroquery) 

Tested with astroquery.simbad 0.4.8 and python 3.11

Exclude objects with mass > 13 Jupiter's mass as probable brown dwarfs
