
# LocatingArrayGenerator
### Authors: Dr. Ryan Dougherty, CDT Dylan Green, CDT Hyunmook Kang, CDT Grace Kim
## Overview
The purpose of this program is to return the number of rows needed to locate given inputs of covering arrays using a multistage genetic algorithm.

The program will append rows for non locating pairs in a Covering Array using a two-stage genetic algorithm until a complete Locating Array is created. The total number of rows needed to make it locating is returned at the end of the program.
## Usage
#### Compiling

To compile the program use the command navigate to the location of :
```
g++ -std=c++20 -O3 LocatingDetectingGeneticAlg.cpp
```
#### Running the Program
Run the ```a.out``` after compilation of ```LocatingDetectingGeneticAlg.cpp``` by running the command: ```./a.out```

Example Output:
```
------------d=1 ./evaluation/2^10-t2_l3------------
Read file with 20 rows.
Computing rows of d-sets...
Ready to look at pairs...
largest_rows_num_map size=6
(3, 21) (4, 48) (5, 50) (6, 39) (7, 15) (8, 7)
There were 99 non-locating pairs
Total N=25, time=285ms.
------------d=1 ./evaluation/2^15-t2_l3------------
Read file with 20 rows.
Computing rows of d-sets...
Ready to look at pairs...
largest_rows_num_map size=6
(3, 45) (4, 98) (5, 136) (6, 95) (7, 45) (8, 1)
There were 378 non-locating pairs
Total N=28, time=2173ms.
```
- The first output explains locating an input Covering Array with parameters: k = 10, t = 2, d = 1, v - 2, lambda = 3.
- The second output explains locating an input Covering Array with parameters: k = 15, t = 2, d = 1, v - 2, lambda = 3.

#### Covering Array Options
1. To use the <strong>Lovász Local Lemma</strong>, set:
```LLL_instead_of_file = ;``` on line 768 to ```LLL_instead_of_file = false;```
2. To use your own custom-generated Covering Arrays, set
```LLL_instead_of_file = ;``` on line 768 to ```LLL_instead_of_file = true;```


## Details and Definitions
<strong>Combinatorial Interaction Testing (CIT)</strong> refers to the creation of test suites that detect or locate desired interactions. Past research has shown success in faster location of specific interactions using a two-stage creation of uniform <strong>Locating Arrays (LAs)</strong>. Given that LAs guarantee unique location for every such set of interactions, we contribute to past research by implementing the ability to generate non-uniform arrays using a multistage genetic algorithm that splits up the second stage of generating locating pairs.

<strong>Covering Arrays (CAs)</strong> A formal model that guarantees all interactions of size at most a specified number,  $t$, appear at least a specified number, &lambda; times.

<strong>Locating Arrays (LAs)</strong> - Test suites that can detect all covered interactions between test variables.