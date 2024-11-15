#!bin/bash

#first, compile: "bash compilefdhelix.sh"

#if sequence is short: "./fdhelix.exe abdna AAAAAAAAAATTTTTTTTTTATATATATATAAAAAAAAAATTTTTTTTTT | grep C1 > ATAT.pdb"
#if sequence is in file: "bash getfdhelix.sh seqfile.txt ATAT.pdb"

var=$(cat $1)

./fdhelix.exe abdna $var > "${2}AA_fdh.pdb"
grep "C1'" "${2}AA_fdh.pdb" > "${2}CG.pdb"

echo "source leaprc.DNA.bsc1
mol = loadpdb ${2}AA_fdh.pdb
savepdb mol ${2}AA_leap.pdb
quit" > leap${2}.in

tleap -f leap${2}.in > leap${2}.out

echo "parm ${2}AA_leap.pdb
trajin ${2}AA_leap.pdb
strip @H=
trajout ${2}AA.pdb
run" > cpptraj${2}.in

cpptraj -i cpptraj${2}.in > cpptraj${2}.out

rm leap${2}.in leap${2}.out cpptraj${2}.in cpptraj${2}.out ${2}AA_fdh.pdb ${2}AA_leap.pdb