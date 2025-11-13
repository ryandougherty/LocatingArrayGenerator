
# LocAG: a fast Locating Array Generator
### Authors: Dr. Ryan Dougherty, CDT Dylan Green, CDT Hyunmook Kang, CDT Grace Kim, Dr. Stephanie Forrest
## Overview
The purpose of this program is to generate a locating array. We use a two-stage approach, where the second uses a genetic algorithm.

The program will append rows for non locating pairs from a given covering array using a two-stage genetic algorithm until a complete Locating Array is created. The total number of rows needed to make it locating is returned at the end of the program.
## Usage
#### Compiling

On Windows, it may be ne

To compile the program use the command navigate to the location of :
```
<!-- g++ -std=c++20 -O3 LocAG.cpp -o ./LocAG -->
g++ -std=c++20 -O3 LocAG.cpp phase1/phase1.cpp phase2/phase2.cpp utils/utils.cpp -o ./LocAG -ltbb
```
#### Running the Program
Run the ```./LoCAG``` program after compilation by running the command: ```./LocAG <name of config>```. The names of real-world systems we have modeled are in the ```LocAG.cpp``` file. 

## Details and Definitions
<strong>Combinatorial Interaction Testing (CIT)</strong> refers to the creation of test suites that detect or locate desired interactions. Past research has shown success in faster location of specific interactions using a two-stage creation of uniform <strong>Locating Arrays (LAs)</strong>. Given that LAs guarantee unique location for every such set of interactions, we contribute to past research by implementing the ability to generate non-uniform arrays using a multistage genetic algorithm that splits up the second stage of generating locating pairs.

<strong>Covering Arrays (CAs)</strong> A formal model that guarantees all interactions of size at most a specified number,  $t$, appear at least a specified number, &lambda; times.

<strong>Locating Arrays (LAs)</strong> - Test suites that can detect all covered interactions between test variables of size at most $d$.
