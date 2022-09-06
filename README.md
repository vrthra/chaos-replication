# Estimating Equivalent Mutants: Are we there yet?

This replication package contains the following:

- results of the manual classification of live mutants (`manual-classification/`)
- automatically generated test suites (`chaos-test-suites.tar.gz`)
- data about test execution in the form of kill matrices
- scripts to compute the estimators
- plots and scripts to generate them


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
