# DataSample-and-Test-Suite
This project is for create the sequence that pass all the NIST test suite and TestU01 batteries and Hamano et al.'s new test 
but will be rejected by our proposed DM-1 and DM-2 test.
/* generate the sample */
The code that generate the sample data is in file GenerateSample.c, which uses the file rdseed.data as input, and generate the 
output file named biggap10_4_1. Compile this c source file and run, then the generated file biggap10_4_1 will be the same as 
sample.data.

/* NIST and TestU01 test suite*/
Unzip these two compression packs， there are README files in each of them. Follow the README files and compile them.
Then you can test the generated file biggap10_4_1(or sample.data directly). Notice that the size of this test data is
1*10^9 bits. And when do the NIST test, the input sequence size should be 10^6, and the number of sequences should be 
1000. When do the TestU01 battery test, just transform the sample size with 10^9 in the correspondding function.

/*test binary file by NIST linear complexity test, HSY test and our proposed DM-1 and DM-2 test*/
These source files are put in the directory ./LinearComplexityRandomTest. Just compile the c source files in this directory.
And then use the command like:
                ./a.out inputfile output.txt
Notice that in the source file LC_Random_Test.c, there are some lines for paralleling by using omp instructions, so when
compiling, you can use parameter -openmp to parallel the codes for high speed.
Then the test result of these 4 test methods will be written in the file output.txt, which included both P and U results.
