There are several steps to build and run this file

Name = Assignment_2_Kunj/Jarrow

1). Load Module
<pre>
module load intel/2022.0
</pre>

2). Build
<pre>
icc -qopenmp -O3 -xHost -qopt-report=5 -qopt-report-phase=vec -std=c++11 main.cpp -o main
</pre>


3). Run
<pre>q
./main
</pre>

4). Tree information
- I implemented a Jarrow-Rudd Rree
- I used backward propagation to price options