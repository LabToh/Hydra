# Example of Cross Validation of qHydra and eqHydra

2018.3.27
Created by Shogo Satoyama

## Usage for Cross Validation
Please type these commands.   
***Please change paths line 22-24 in eqHydraCrossValidation.cpp and qHydraCrossValidation.cpp***
### Compile Option
***qHydra Compile Option***
```
g++ -lafopencl qHydraCrossValidation.cpp -O3 -std=c++11 \
  -I/path/to/bayesopt/include \
  -I/path/to/bayesopt/utils \
  -I/path/to/bayesopt/matplotpp \
  -Wl,-rpath,/path/to/bayesopt/lib \
  /path/to/bayesopt/lib/libbayesopt.a \
  /path/to/bayesopt/lib/libnlopt.a \
  -o qHydra
```
***eqHydra Compile Option***
```
g++ -lafopencl eqHydraCrossValidation.cpp -O3 -std=c++11 \
  -I/path/to/bayesopt/include \
  -I/path/to/bayesopt/utils \
  -I/path/to/bayesopt/matplotpp \
  -Wl,-rpath,/path/to/bayesopt/lib \
  /path/to/bayesopt/lib/libbayesopt.a \
  /path/to/bayesopt/lib/libnlopt.a \
  -o eqHydra
```
### Run

***qHydra***    
`./qHydra`    

***eqHydra***    
`./eqHydra`  
