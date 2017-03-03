#Prime Factorization

This projects aims to find two prime factors if they exist for any number given. It utilizes gmp for large number manipluation, calculating primes, and verifying prime. It also utilizes mpi for distributed processing. They primary target is Unix-like systems, however only minor porting would be required for a Windows system.

To complile make sure you having a working mpi environment and gmp installation. The most basic compilation can be done with:

    mpicc -o factor factor.c -lgmp -lmpi

To later run the program use mpirun

    mpirun -np 4 factor 14

np is the number of processes desired

factor is the path to your binary

The final value is the number to be factored

The program will output the results to stdout so you have a variety of options for collecting the output such as: piping to a file, running within a tmux/screen session, nohup, etc.
