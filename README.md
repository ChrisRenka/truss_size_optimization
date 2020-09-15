# truss_size_optimization
Codes for truss size optimization

Folder "EB-A-SHADE" contains the codes for the algorithm.
Folder "fit" contains the codes for setting up the structures and Finite Element Analysis.

To select the structre, modify the variable "func", at line 42 of the file EB-A-SHADE/mainLSHADE.cpp.
Supported values are:

* 10: 10 bar truss;
* 17: 17 bar truss;
* 25: 25 bar truss;
* 721: 72 bar truss, case 1;
* 722: 72 bar truss, case 2;
* 120: 120 bar truss;
* 200: 200 bar truss.

File input.txt contains the values for the algorithm parameters. These are already set up according to the experimental settings. If you wish to reproduce the results from the folder "tests_article", keep the values as is.

**You need to create a folder called "tests" inside the folder EB-A-SHADE. Otherwise, there will be a segmentation fault when you try to run the program.**

After setting the desired value for "func" and creating the "tests" folder, to run the program:

```
cd EB-A-SHADE
make
./main <runs>
```
The value "runs" is the number of independent runs you wish to execute. For 30 runs, for example, use:
  
```
./main 30  
```
