%MEM=10000MB
%nproc=8
%Chk=h2.chk
#P CCSD(T)/aug-cc-pVQZ SCF=(CONVER=8,maxcycle=200,XQC) scan

HEME + H2 2D Scan

0 1
 h
 h   1 dh2

dh2          0.5  S 100 0.01

