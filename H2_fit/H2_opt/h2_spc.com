%MEM=10000MB
%nproc=8
%Chk=h2_spc.chk
#P CCSD(T)/aug-cc-pVQZ SCF=(CONVER=8,maxcycle=200,XQC) 

HEME + H2 2D Scan

0 1
 h   0.0 0.0  0.371193
 h   0.0 0.0 -0.371193 

