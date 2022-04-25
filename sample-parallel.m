(*Limits Mathematica to requested resources*)
Unprotect[$ProcessorCount];$ProcessorCount = 8;
 
(*Prints the machine name that each kernel is running on*)
Print[ParallelEvaluate[$MachineName]];

(*Prints all Mersenne PRime numbers less than 2000*)
Print[Parallelize[Select[Range[2000],PrimeQ[2^#-1]&]]];
