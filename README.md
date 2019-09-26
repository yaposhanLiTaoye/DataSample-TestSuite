# DataSample-and-Test-Suite
This project is for generating a sequence that passs all the NIST test suite and TestU01 batteries and Hamano et al.'s new test 
but will be rejected by our proposed DM-1 and DM-2 test.

##1. Generate the Sample Sequence 

   The code for generating the sample data is in the file "GenerateSample.c", in which uses the data file "rdseed.data" as input. The 
output file is named "biggap10_4_1". Compile this single c source file and run it, then you will get the generated file biggap10_4_1,
which is the same as "sample.data".

##2. Test the Generated Sample Data by NIST and TestU01 Test Suite

   Unzip these two compression packs，there are README or INSTALL files in both of them. Follow the README or INSTALL files to compile them,and you will get the corresponding executable files. Then test the generated file "biggap10_4_1". Notice that the size of this data file is10^9 bits. When do the NIST test, the input sequence size should be 10^6, and the number of sequences should be 1000. While for TestU01 batteries, the input size is simply 10^9.

##3. Test the Generated Sample Data by NIST Linear Complexity Test, HSY Test and Our Proposed DM-1 and DM-2 Test

   These source files are placed in the directory ./LinearComplexityRandomTest. Enter this directory and compile these c source files there with a command like:

            gcc -o FileTest *.c -lm -fopenmp 
   
   Notice that the input parameter "-fopenmp" in this command is for paralleling, because in the file "LC_Random_Test.c" there are some codes written for paralleling by the instruction of omp. Then there will be a executable file named "FileTest", and you can now test the sample data by
           
           ./FileTest ../biggap10_4_1 TestResult.txt
   
   The test result of the binary file "biggap_10_4_1" by these 4 test methods then will be written in the file TestResult.txt, which includes  both P and U results and the number of the final rejected sequences according to the criteria of U for each test.
