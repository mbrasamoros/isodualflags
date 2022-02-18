# seeds-algorithm
This repository hosts the implementation of the seeds-algorithm to compute the number of numerical semigroups of each given genus.

For a quick survey on the sequence of the number of numerical semigroups of each given genus and the conjectures related to it see the entry http://oeis.org/A007323 in the On-line Encyclopedia of Integer Sequences.

One can compile the non-parallelized version as follows.

    gcc -Wall -o seeds-algorithm-noparallelization.out seeds-algorithm-noparallelization.cpp

and then execute it to compute the number of numerical semigroups of a given genus. For instance, to obtain the number of semigroups of genus 35 the command is

    ./seeds-algorithm-noparallelization.out 35

and the output will be

    n35=66687201 (without parallelization) time taken 1
    
Alternatively, one can compile the Cilk++ parallelized version as follows.

    g++ -fcilkplus -flto -Ofast -march=native -mtune=native -fwhole-program -o seeds-algorithm.out seeds-algorithm.cpp

and then execute it to compute the number of numerical semigroups of a given genus. For instance, to obtain the number of semigroups of genus 40 the command is

    ./seeds-algorithm.out 40

and the output will be

    n40=774614284 (24 workers) time taken 514 miliseconds
