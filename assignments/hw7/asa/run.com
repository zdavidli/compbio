#Accessible Surface Area Calculation
#mpprun -map [Group #] -np 1 /usr/local/bin/python2.2 -OO asa.py A.pdb > A.asa
#mpprun -map [Group #] -np 1 /usr/local/bin/python2.2 -OO asa.py B.pdb > B.asa
