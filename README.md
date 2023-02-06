# Empirical Evaluation of Frequency Based Statistical Models for Estimating Killable Mutants

This replication package contains the following:
- results of the manual classification of live mutants (`manual-classification/`)
- automatically generated test suites (`chaos-test-suites.tar.gz`)
- data about test execution in the form of kill matrices (`killmatrix/`)
- scripts to compute the estimators
- plots and scripts to generate them
- A extended version of the original paper


## Automatically Generated Test Suites
The generated test suites (along with the original ones) are stored in the `chaos-test-suites.tar.gz` file.

We generated the test suites using EvoSuite (1.7.0-SNAPSHOT).

We generated the `random` test suites using the following parameters:

```
     -Dsearch_budget=60
     -Duse_separate_classloader=false
     -generateRandom
     -Dassertion_strategy=all
     -Dsandbox=TRUE
     -Dassertions=TRUE
     -Dassertion_timeout=120
     -Dminimization_timeout=300 
```
We generated the `dynamosa` test suites using the following parameters:

```
    -Dsearch_budget=60 \
    -generateMOSuite \
    -Duse_separate_classloader=false \
    -Dcriterion=BRANCH:LINE:WEAKMUTATION:INPUT:OUTPUT:METHOD:EXCEPTION \
    -Dpopulation=50 \
    -Dranking_type=PREFERENCE_SORTING \
    -Dalgorithm=DYNAMOSA \
    -Dassertion_timeout=120 \
    -Dcoverage=TRUE \
    -Dsandbox=TRUE \
    -Dassertions=TRUE \
    -Dminimization_timeout=300 \
```


## Data about Test Execution in the Form of Kill Matrices
The generated kill matrices are stored in the `killmatrix/killmatrix.tar.gz` file.
They are represented as a compresses sparse row (csr) matrix and stored in npz format.
The matrices are generated from the parsed JUGE output 
(since it is very large, we don't provide it in this bundle, an example is in 
`./raw/commons-csv.random.03-killmatrix.tar.gz`).  

## Scripts to Compute the Estimators
The scripts are stored in the `scripts` folder.  
To compute estimations run

```
Rscript scripts/estimation.R -i ./killmatrix -o ./data/summary.csv --cores 8 --total ./data/total.csv --ci
```

Other parameters can be printed with `Rscript scripts/estimation.R -h`.  
All required R packages will be installed automatically.
Note, the `reticulate` package require Python compiled with shared library support.
For instance, it can be installed with [pyenv](https://github.com/pyenv/pyenv) by executing `env PYTHON_CONFIGURE_OPTS="--enable-shared" pyenv 
install 3.8.2`.
`--ci` parameter controls the computation of confidence intervals for `PCG`, `PNPMLE`, and `UNPMLE` estimators from the `SPECIES` R package.
It is slow and unstable, restart it if it hangs or exits with the timeout. In that case partial results can be found in `./logs` folder.
After the estimators are computed to add manual estimations use  
`python scripts/postprocessing.py manual -s ./data/summary.csv --classified ./manual-classification`  
Next, one can compute mean difference across subjects with  
`python scripts/postprocessing.py mean -s ./data/summary.csv -o ./data/mean.csv --mean-type classes -t organic`

## Plots and Scripts to Generate them
The plots are stored in the `plots` folder.  
To regenerate them use
`python scripts/postprocessing.py plot -c estimators -s ./data/summary.csv 
-o ./plots/estimators` and
`python scripts/postprocessing.py plot -c subjects -s ./data/summary.csv
-o ./plots/subjects`
