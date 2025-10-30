#include <vector>
#include <map>
#include <string>

// A type to represent a t-way interaction. 
// You could use a string "col1:val1,col2:val2" or a more complex struct.
typedef std::string Interaction;

// A map to store the number of times each interaction has been covered so far.
// This is the main data structure for tracking coverage.
std::map<Interaction, int> coverageCounts;

// We also need a master list of all possible t-way interactions.
std::vector<Interaction> allInteractions; 

// The final array
std::vector<std::vector<int>> array;