# Catmull-Rom interpolation
The code tabulates the splines for the 4D mechanism computed directly using Sundials IDA.
The code can be reused to evaluate splines in Fluent.

This code proves that an efficient implementation of a lookup table is ~8x faster than an efficient solution using IDA from sundials.
On average, building the table took 15.819 s for breaks [20,20,20,20]. 

Average time taken per inference iteration is 0.012010 ms 
Average time taken by sundials is 0.093376 ms 
RMSE is 0.000620 
This is a 7.77x speed up factor for a very low RMSE.

See the output file for evaluation output. The columns are defined as:

interpolated|exact|T|xNH3|xO2|XNO

Anton Fadic
