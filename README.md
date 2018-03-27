# qHydra and eqHydra

2018.3.27
Created by **Shogo Satoyama**

## About
This program recomends genes to desease based on the quadripartite graph.

### Data Set
Following three adjacency matrices are used as input data of our programs qHydra and eqHydra.

The first matrix is designed as a adjacency matrix of connectors (C) × genes (G).
Each row of the matrix corresponds to a connector (a Gene Ontology term or a pathway).Each Column corresponds to a gens.
***Note that the first parts of the columns are occupied by disease related genes (DG) (see InputData.png).***
Then, genes where association with any diseases are located of the following columns.That is, the matrix consists of two matrices , A (C × DG) and B (C × G\DG), where G\DG means a set of genes without DG.

The second matrix is designed as a transpose of A of the first matrix.
***That is, the order of the disease related gene of the rows of the second matrix is identical to that of the columns of the first matrix.***

The third matrix is designed as an adjacency matrix of diseases (D) × disease related genes (DG).
That is, each row of the matrix corresponds to a disease, where as column of the matrix corresponds to a disease related gene.

## Dependencies
**C++ Libraries**    
ArrayFire https://arrayfire.com/   
BayesOpt https://rmcantin.bitbucket.io/html/index.html

**R Packages**   
org.Hs.eg.db https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html   
annotate https://www.bioconductor.org/packages/release/bioc/html/annotate.html    
igraph https://cran.r-project.org/web/packages/igraph/index.html    

Please install these tools.



## Usage for Cross Validation
Please type these commands.   
***Please change paths where is line 22-24 in eqHydraCrossValidation.cpp and qHydraCrossValidation.cpp***

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



***Please change line 853,854 and 872 in qHydraCrossValidation.cpp and line 970,971 and 990 in eqHydraCrossValidation.cpp if you use Gene Ontolog as connector.***   

**qHydraCrossValidation.cpp**   
line 853: Matrix<int> first_bipartite("Pathway_Gene.txt"); -> Matrix<int> first_bipartite("GeneOntology_Gene.txt");   
line 854: Matrix<int> second_bipartite("DiseaseRelatedGene_Pathway.txt"); -> Matrix<int> second_bipartite("DiseaseRelatedGene_GeneOntology.txt");   
line 872: qHydrabo.OutFileNameSetting("qHydra(pathway)"); -> qHydrabo.OutFileNameSetting("qHydra(GO)");   

**eqHydraCrossValidation.cpp**   
line 970: Matrix<int> first_bipartite("Pathway_Gene.txt"); -> Matrix<int> first_bipartite("GeneOntology_Gene.txt");   
line 971: Matrix<int> second_bipartite("DiseaseRelatedGene_Pathway.txt"); -> Matrix<int> second_bipartite("DiseaseRelatedGene_GeneOntology.txt");   
line 990: eqhydrabo.OutFileNameSetting("eqHydra(pathway)"); -> eqhydrabo.OutFileNameSetting("eqHydra(GO)");   


## Usage for Computing Precision and Recall Enhancement
Please type these commands.

### Compile Option

***Computing Precision and Recall Enhancement***
```
g++ -O3 -std=c++11 eqHydraRecallPrecision.cpp -o PrecisonAndRecallEnhancement
```   
Then, please run PrecisonAndRecallEnhancement program.   
`./PrecisonAndRecallEnhancement file_name length`   

If you compute precision and recall enhancement which length of recommendation list is 5, you should run the PrecisonAndRecallEnhancement program like this.   
`./PrecisonAndRecallEnhancement "eqHydra(Pathway)FirstLambdaXXXXXSecondLambdaYYYYYThirdLambdaZZZZZ.txt" 5`



## Usage for Prediction
Please type these commands.

### Compile Option
```
g++ -lafopencl eqHydraAllScore.cpp -O3 -std=c++11 -o AllScore
```
### Run
`./AllScore`   
and Run R Progarams   
`TopN.R`   
`Plot.R`

***Please change paths where is line 18 in TopN.R and line 41-42 in Plot.R***
