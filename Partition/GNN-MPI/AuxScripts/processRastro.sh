aky_converter --no-links rastro*.rst > rastro.paje
pj_dump rastro.paje | grep ^State > rastro.csv
AuxPrograms/mpitimecalc rastro.csv mpitime.csv

rm rastro.csv
rm rastro.paje
rm rastro*.rst

