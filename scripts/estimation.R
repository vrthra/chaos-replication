##########
packages = c("stringr", "argparse", "SpadeR", "SPECIES", "tools", "glue",
             "foreach", "doParallel", "data.table", "rlist", "remotes",
             "wiqid", "reticulate", "Matrix", "jsonlite")

library("wiqid")
library("jsonlite")

#renv::hydrate()
load_install = function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg, repos = "http://cran.us.r-project.org")
}

for (pkg in packages) {
  load_install(pkg)
}
# reticulate::use_python("~/.pyenv/shims/python3", required = TRUE)
reticulate::py_config()
if (!require("pathlibr")) remotes::install_github("mstr3336/pathlibr")
##########
library(stringr)
library(argparse)
library(SpadeR)
library(SPECIES)
library(tools)
library(glue)
library(pathlibr)
library(foreach)
library(doParallel)
library(data.table)
library(rlist)
library(Matrix)

all_subjects = c("commons-beanutils",
                 "commons-cli",
                 "commons-codec",
                 "commons-collections",
                 "commons-compress",
                 "commons-configuration",
                 "commons-csv",
                 "commons-dbcp",
                 "commons-dbutils",
                 "commons-digester",
                 "commons-email",
                 "commons-exec",
                 "commons-fileupload",
                 "commons-functor",
                 "commons-imaging",
                 "commons-io",
                 "commons-lang",
                 "commons-math",
                 "commons-net",
                 "commons-pool",
                 "commons-validator")

default_subjects = c("commons-csv",
                     "commons-collections",
                     "commons-compress",
                     "commons-configuration",
                     "commons-dbcp",
                     "commons-imaging",
                     "commons-io",
                     "commons-lang",
                     "commons-net",
                     "commons-math")

all_testsuites = c("organic", "evosuite", "randoop", "dynamosa", "random")
default_testsuites = c("organic", "dynamosa", "random")

estimators_fun = list(chao_lower = 'get_lowerbound',
                      bootstrap = 'estimator_bootstrap', # should be fixed
                      sample_coverage = 'get_sample_coverage_estimate', # sample coverage
                      chao_upperbound = 'get_upperbound',
                      jack_our = 'get_custom_jackknife_estimate',
                      zelterman = 'estimator_zelterman',
                      chao2 = 'get_chao2',
                      chao2_corrected = 'get_chao2_bias_corrected',
                      species_pcg = 'estimator_pcg',
                      species_unpmle = 'estimator_unpmle',
                      species_pnpmle = 'estimator_pnpmle',
                      species_chaobunge = 'estimator_ChaoBunge',
                      species_jackknife = 'estimator_jackknife',
                      spader = 'estimators_spader')

# Encode total number of mutants per program
all_estimators = names(estimators_fun)

default_estimators = c('spader',
                       'bootstrap',
                       'zelterman',
                       'species_pcg',
                       'species_unpmle',
                       'species_pnpmle',
                       'species_chaobunge',
                       'species_jackknife',
                       'sample_coverage'
)

counts_estimators = c('chao_lower',
                      'sample_coverage',
                      'chao_upperbound',
                      'jack_our',
                      'zelterman',
                      'chao2',
                      'chao2_corrected',
                      'species_pcg',
                      'species_unpmle',
                      'species_pnpmle',
                      'species_chaobunge',
                      'species_jackknife',
                      'spader'
)
matrix_estimators = list('bootstrap') # necessary to perform bootstrapping for CI
estimators_additional_param = list(chao_upperbound = c('upper_bound'),
                                   species_pcg = c('ci'),
                                   species_unpmle = c('ci'),
                                   species_pnpmle = c('ci')
)

all_catch_types = c('classes', 'methods', 'assertions', 'delimiters')
default_catch_types = c('classes', 'methods')

item_or_default = function(list_, item, default = NULL) {
  value = list_[[item]]
  if (is.null(value)) return(default) else return(value)
}

run_estimator <- function(estimator, params) {
  tryCatch({
    res = do.call(estimator, params) # one of predict_population_size estimator
    return(res)
  },
    error = function(cond) {
      cat(glue_collapse(cond, sep = "\n"), file = "mylog.txt", append = TRUE)
      print(cond)
      message = return(list(N = -1, lowCI = -1, uppCI = -1))
    }
  )
}

# load sparse matrix produced by scipy.sparse
load_npz <- function(f_name) {
  np = import("scipy.sparse")
  smatrix = np$load_npz(f_name)
  n_mutants = dim(smatrix)[1]
  smatrix = smatrix[, colSums(smatrix) != n_mutants] # remove tests that always fail
  return(smatrix)
}

load_csv_data <- function(f_name) {
  matrix = fread(f_name)
  if (max(matrix) > 1) { # remove index column
    matrix = matrix[, -1]
  }
  return(matrix)
}

load_data <- function(datafile, data_type) {
  if (!file.exists(datafile)) {
    print(glue("Err: no data file found at {datafile}"))
    return(NULL)
  }
  if (data_type == "abundance" || data_type == "freq") {
    d = read.table(datafile, sep = " ", header = FALSE)
  }else if (data_type == "csv") {
    d = load_csv_data(datafile)
  }else if (data_type == "npz") {
    d = load_npz(datafile)
  }
  if (sum(d) == 0) {
    print("WARN: data contains only zeros")
    return(NULL)
  }
  if (max(d) > 1 && data_type == "npz") {
    d = sign(d)
  }
  return(d)
}

get_vector_counts <- function(data, v_size = length(data)) {
  # data - matrix of 0/1
  # return vector of counts including not-found mutants i.e. f[0]
  counts = as.data.frame(table(factor(data, levels = c(0:v_size))))
  return(counts$Freq)
}

get_matrix_counts <- function(data) {
  # data - matrix of 0/1
  # return vector of counts including not-found mutants i.e. f[0]
  n = ncol(data)
  counts = as.data.frame(table(factor(rowSums(data), levels = 0:n)))
  return(counts$Freq)
}

load_counts <- function(data_type, d) {
  if (data_type == "abundance") {
    # deprecated
    counts = get_vector_counts(d$V1)
    counts_no0 = counts[-1]
    n_samples = sum(d$V1) #length(d) #stub for abundance
    Sn = length(d[d > 0])
    total_mutants = length(d$V1)
  }else if (data_type == "freq") {
    # deprecated
    counts = get_vector_counts(d$V1[-1])
    counts_no0 = counts[-1]
    n_samples = dd$V1[1]
    Sn = length(dd$V1[dd$V1 > 0]) - 1 # minus d[1]
    total_mutants = 0
  }else if (data_type == "csv" || data_type == "npz") {
    counts = get_matrix_counts(d)
    counts_no0 = counts[-1]
    n_samples = dim(d)[2]
    rs = rowSums(d)
    total_mutants = nrow(d)
    Sn = length(rs[rs > 0])
  }
  return(list(counts = counts_no0, n_samples = n_samples, Sn = Sn, total_mutants = total_mutants))
}


get_npz_file <- function(data_dir, catch_type, subject, testsuite, run = "01") {
  return((data_dir %//%
    subject %//%
    glue("{subject}_{testsuite}_{catch_type}-killmatrix.npz"))$show)
}

get_abundance_file <- function(data_dir, catch_type, subject, testsuite, run = "01") {
  return((data_dir %//%
    subject %//%
    glue("{subject}_{testsuite}_{catch_type}-killmatrix.abundance"))$show)
}

get_npz_file_old <- function(data_dir, catch_type, subject, testsuite, run = "01") {
  return((data_dir %//%
    catch_type %//%
    subject %//%
    glue("{testsuite}-{run}-killmatrix.npz"))$show)
}

compute_estimator <- function(counts_no0, identity_matrix, estimator, n_samples, total_mutants, ci) {
  est_type = get_est_type(estimator)
  if (est_type == 'counts')
    estimation = estimate_from_counts(estimator, counts_no0, n_samples, upper_bound = total_mutants, ci = ci)
  else if (est_type == 'matrix')
    estimation = estimate_from_matrix(estimator, identity_matrix)
  else
    estimation = NULL
  if (is.data.frame(estimation)) {
    return(estimation) # if we compute many estimations at once (SpadeR) just return the results
  }
  estN = if (!is.null(estimation$N)) estimation$N else -1
  # if no confidence interval is given, return point estimation
  lowCI = if (!is.null(estimation$lowCI)) estimation$lowCI else estN
  uppCI = if (!is.null(estimation$uppCI)) estimation$uppCI else estN
  res_df = list(estimator = estimator,
                estimation = estN,
                CILower = lowCI,
                CIUpper = uppCI)
  return(res_df)
}

compute_estimators <- function(subject, testsuite, estimators, data_dir, catch_type, ci, data_type, log_dir, total_mutants_df = NULL) {
      if (data_type == 'npz') {
    data_file = get_npz_file(data_dir, catch_type, subject, testsuite)
  }else if (data_type == 'abundance') {
    data_file = get_abundance_file(data_dir, catch_type, subject, testsuite)
  }
  identity_matrix = load_data(data_file, data_type)
  if (is.null(identity_matrix)) {
    return(NULL)
  }
  counts_obj = load_counts(data_type, identity_matrix)
  counts = counts_obj$counts
  n_samples = counts_obj$n_samples
  Sn = counts_obj$Sn
  total_mutants = counts_obj$total_mutants # if total mutants is not available, take form killmatrix
  if (!is.null(total_mutants_df)) {
    total_mutants = total_mutants_df[which(total_mutants_df$subject == subject & total_mutants_df$testsuite == testsuite), 'total_mutants']
  }
  res_df = list()
  for (estimator in estimators) {
    file.create(as.character(log_dir %//% glue("{estimator}_{subject}_{testsuite}_{catch_type}_log.txt")))
    res = compute_estimator(counts, identity_matrix, estimator, n_samples, total_mutants = total_mutants, ci) # error
    rr = do.call("cbind", res)
    rr = cbind(subject = subject,
               testsuite = testsuite,
               type = catch_type,
               total_mutants = total_mutants,
               killed = Sn,
               rr)
    fwrite(as.data.table(rr), as.character(log_dir %//% glue("{estimator}_{subject}_{testsuite}_{catch_type}_log.txt")))
    res_df = append(res_df, list(res))
  }
  all_data = do.call("rbind", res_df)
  all_data = cbind(subject = subject,
                   testsuite = testsuite,
                   type = catch_type,
                   total_mutants = total_mutants,
                   killed = Sn,
                   all_data)
  return(all_data)
}

get_est_type <- function(estimator) {
  if (estimator %in% counts_estimators) {
    return('counts')
  }
  if (estimator %in% matrix_estimators) {
    return('matrix')
  }
  return(NULL)
}

estimate_from_matrix <- function(estimator, m) {
  fun_name = estimators_fun[[estimator]]
  additional_params = item_or_default(estimators_additional_param, estimator, list())
  params = append(list(data = m), additional_params)
  return(run_estimator(fun_name, params))
}

get_additional_params <- function(estimator, ...) {
  args = list(...)
  params = list()
  params_for_est = item_or_default(estimators_additional_param, estimator, c())
  for (a in names(args)) {
    if (a %in% params_for_est) {
      params[a] = args[a]
    }
  }
  return(params)
}

estimate_from_counts <- function(estimator, counts, n_samples, ...) {
  fun_name = estimators_fun[[estimator]]
  additional_params = get_additional_params(estimator, ...)
  params = append(list(n_samples = n_samples, counts = counts), additional_params)
  return(run_estimator(fun_name, params))
}

estimators_spader <- function(n_samples, counts, ...) {
  est = list(
    "Chao2 (Chao, 1987)" = "spader_chao",
    "Chao2-bc" = "spader_chao2bc",
    "iChao2 (Chiu et al. 2014)" = "spader_ichao2",
    "ICE (Lee & Chao, 1994)" = "spader_ice",
    "ICE-1 (Lee & Chao, 1994)" = "spader_ice1")
  est_values = names(est)
  incidence_freq_count = c(n_samples, unlist(Map(c, seq_along(counts), counts)))
  # k - the cut-off point (default = 10)
  res = SpadeR::ChaoSpecies(incidence_freq_count, datatype = "incidence_freq_count", k = 10, conf = 0.95)$Species_table
  #SpadeR::SpecInciiChao2()
  res = res[trimws(rownames(res)) %in% est_values,]
  res$estimator = sapply(trimws(rownames(res)), function(x)est[[x]])
  rownames(res) = NULL
  res = res[, c("estimator", "Estimate", "95%Lower", "95%Upper")]
  setnames(res, old = c("Estimate", "95%Lower", "95%Upper"), new = c("estimation", "CILower", "CIUpper"))
  return(res)
}

# penalized conditional NPML estimator of the species number by Wang and Lindsay 2005
estimator_pnpmle <- function(counts, ci, ...) {
  # t is the cutoff value to define the relatively less abundant species
  # C 1 to calculate confidence interval
  # b the number of bootstrap samples
  c = if (ci) 1 else 0
  freq_of_freq = data.frame(freq = seq_along(counts), nj = counts)
  res = SPECIES::pnpmle(freq_of_freq, t = 15, C = c, b = 100, seed = NULL, conf = 0.95, dis = 0)
  return(list(N = res$Nhat, lowCI = res$CI[1], uppCI = res$CI[2]))
}


# unconditional NPML estimator of the species number by Norris and Pol- lock 1996, 1998
estimator_unpmle <- function(counts, ci, ...) {
  # t is the cutoff value to define the relatively less abundant species
  # C 1 to calculate confidence interval
  # b the number of bootstrap samples
  c = if (ci) 1 else 0
  freq_of_freq = data.frame(freq = seq_along(counts), nj = counts)
  res = SPECIES::unpmle(freq_of_freq, t = 15, C = c, method = "W-L", b = 100, conf = .95, seed = NULL, dis = 1)
  return(list(N = res$Nhat, lowCI = res$CI[1], uppCI = res$CI[2]))
}

# jackknife estimator of the species number by Burnham and Overton 1978 and 1979
estimator_jackknife <- function(counts, ...) {
  freq_of_freq = data.frame(freq = seq_along(counts), nj = counts)
  res = SPECIES::jackknife(freq_of_freq, k = 5, conf = 0.95)
  return(list(N = res$Nhat, lowCI = res$CI[1], uppCI = res$CI[2]))
}

# coverage-duplication based estimator from a Poisson-Gamma model by Chao and Bunge 2002
estimator_ChaoBunge <- function(counts, ...) {
  # t is the cutoff value to define the relatively less abundant species
  # C 1 to calculate confidence interval
  # b the number of bootstrap samples
  freq_of_freq = data.frame(freq = seq_along(counts), nj = counts)
  res = SPECIES::ChaoBunge(freq_of_freq, t = 10, conf = 0.95)
  return(list(N = res$Nhat, lowCI = res$CI[1], uppCI = res$CI[2]))
}

# Poisson-compound Gamma estimators of the species number by Wang 2010
estimator_pcg <- function(counts, ci, ...) {
  # t is the cutoff value to define the relatively less abundant species
  # C 1 to calculate confidence interval
  # b the number of bootstrap samples
  c = if (ci) 1 else 0
  freq_of_freq = data.frame(freq = seq_along(counts), nj = counts)
  res = SPECIES::pcg(freq_of_freq, t = 35, C = c, b = 100, alpha = 1:10, seed = NULL, conf = 0.95, dis = 0)
  return(list(N = res$Nhat, lowCI = res$CI[1], uppCI = res$CI[2]))
}

estimator_bootstrap <- function(data, t = 100) {

  bootstrap <- function(data) {
    n = ncol(data)
    # k = nrow(data)
    idx = sample(n, size = n, replace = TRUE)
    S0 = length(which(rowSums(data) > 0)) # number of caought mutants
    B = data[, idx]
    Y = rowSums(B)
    Bs = sum((1 - Y / n)^n)
    Bn = S0 + Bs
    # return(list(Bn, B))
    return(Bn)
  }

  # (Smith and van Belle 1984)
  # it works with incidence data(ID), not with counts
  # we sample n columns from ID with replacement from n samples
  data = data[rowSums(data) > 0,] # reduce matrix size by getting only caught mutants
  bres = replicate(t, bootstrap(data), simplify = TRUE)
  bn = mean(bres)
  varbn = 1 / (t - 1) * sum((bres - bn)^2)
  se = sqrt(varbn)
  return(list(N = bn, lowCI = bn - 1.96 * se, uppCI = bn + 1.96 * se))
}

# alternative zelterman from Bohning (poisson model)
# http://www.personal.soton.ac.uk/dab1f10/Boehning_Dublin_2008.pdf
# variance from https://doi.org/10.1016/j.stamet.2007.10.001
estimator_zelterman <- function(n_samples, counts) {
  n = sum(counts) # killed mutants
  f1 = counts[1]
  f2 = counts[2]
  lambda = 2 * f2 / f1
  N_est = n + n / (exp(lambda) - 1) # which is the same as = n/(1-exp(-lambda))
  G = exp(-lambda) / (1 - exp(-lambda))^2
  var = n *
    G *
    (1 + n * G * lambda^2 * (1 / f1 + 1 / f2))
  # var = n*G*(1+(lambda+2)/(1-exp(-lambda))) # from zelterman, can be wrong according to Bohning
  se = sqrt(var)
  lowCI = N_est - 1.96 * se
  uppCI = N_est + 1.96 * se
  return(list(N = N_est, se = se, lowCI = lowCI, uppCI = uppCI))
}

get_lowerbound <- function(n_samples, counts) {
  # deprecated
  # chao1987estimating.pdf Estimating the Population Size for Capture-Recapture Data with Unequal Catchability
  # An Improved Nonparametric Lower Bound of Species Richness via a Modified Good–Turing Frequency Formula 2014

  S_obs <- sum(counts)
  n <- S_obs
  f1 <- counts[1]
  f2 <- counts[2]
  f3 <- counts[3]
  f4 <- counts[4]
  var_Schao1 <- f2 * (1 / 4 * ((n - 1) / n)^2 * (f1 / f2)^4 +
    ((n - 1) / n)^2 * (f1 / f2)^3 +
    1 / 2 * (n - 1) / n * (f1 / f2)^2)
  S_chao1 <- S_obs + (n - 1) / n * f1^2 / (2 * f2)
  R <- exp(1.96 * (1 + var_Schao1 / (S_chao1 - S_obs)^2)^0.5) # FIXME: error in their formula
  #R <- exp(1.96*(log(1+var_Schao1/(S_chao1 - S_obs)^2))^0.5) #<- from reviewer
  S_lower <- S_obs + (S_chao1 - S_obs) / R
  S_upper <- S_obs + (S_chao1 - S_obs) * R

  return(list(N = S_chao1, se = -1, lowCI = S_lower, uppCI = S_upper))
}

get_sample_coverage_estimate <- function(n_samples, counts) {
  k = length(counts)
  f = counts
  r = sum(counts)
  sm = sum(f * c(1:k))
  smi = sum(f * c(1:k) * c(0:(k - 1)))
  C = 1 - f[1] / (sm)
  gamma2 = max(0, r / C * smi / ((k - 1) * sm^2) * k - 1)
  Nsc = r / C + f[1] * gamma2 / C
  H = array(, k)
  if (0 <= 1 - r / (C * (k - 1) * sm^2) * k * smi) {
    H[1] = gamma2 / C + ((gamma2 * f[1] + r) * sm -
      f[1]^2 * gamma2 -
      r * f[1]) / (C^2 * sm^2)
  }else {
    H[1] = -f[1] * r * (2 * C * sm - sm + f[1]) / (C^3 * (k - 1) * sm^4) *
      k *
      smi +
      gamma2 / C +
      ((gamma2 * f[1] + r) * sm -
        f[1]^2 * gamma2 -
        r * f[1]) / (C^2 * sm^2)
  }
  for (j in c(2:k)) {
    if (0 <= 1 - r / (C * (k - 1) * sm^2) * k * smi) {
      H[j] = f[1] / (C^2 * sm^2) * j * (gamma2 * f[1] + r)
    }else {
      H[j] = -f[1] * r * j / (C^2 * sm^2) - f[1]^2 * gamma2 * j / (C^2 * sm^2) + f[1] / C * k / (k - 1) *
        r *
        (j^2 / C / sm^2 -
          j / (C * sm^2) -
          2 * j * smi / (C * sm^3) -
          j * f[1] / (C^2 * sm^4) * smi)
    }
  }
  cov = outer(f, f) * (-1 / Nsc)
  diag(cov) = f * (1 - f / Nsc)
  varNsc = sum(H %*% cov %*% H) # matrix of 1 element, converting to scalar
  se = sqrt(varNsc)
  lowCI = Nsc - se * 1.96
  uppCI = Nsc + se * 1.96
  return(list(N = Nsc, se = se, lowCI = lowCI, uppCI = uppCI))
}

get_upperbound <- function(n_samples, counts, upper_bound, alpha = 0.05) {
  # Estimating the Richness of a Population When the Maximum Number of Classes
  # Is Fixed: A Nonparametric Solution to an Archaeological Problem
  # 2012
  U = upper_bound
  chao_est = get_chao2(n_samples, counts)
  s_hat = chao_est$N
  se = chao_est$se
  var_s = se * se
  s_obs = sum(counts)
  mu_y = log(s_hat - s_obs)
  sigma2 = log(1 + var_s / (s_hat - s_obs)^2)
  sigma = sqrt(sigma2)
  p = pnorm((log(U - s_obs) - mu_y) / sigma)
  z_p_alpha = qnorm(p * alpha / 2)
  z_p_1alpha = qnorm(p * (1 - alpha / 2))
  s_lower = s_obs + (s_hat - s_obs) * exp(sigma * z_p_alpha)
  s_upper = s_obs + (s_hat - s_obs) * exp(sigma * z_p_1alpha)
  N = s_hat
  # if (s_hat>upper_bound)
  N = (s_upper+s_lower)/2 # sometimes CI doesn't contain point estimate, what if we take a middle point as a new estimate?
  return(list(N = N, se = -1, lowCI = s_lower, uppCI = s_upper))
}

get_custom_jackknife_estimate <- function(n_samples, counts) {
  # not used any more
  k = length(counts)
  f = counts
  S = sum(counts)
  a = matrix(1, k, 5)
  a[1, 1] = 1 + (k - 1) / k
  a[1, 2] = 1 + (2 * k - 3) / k
  a[2, 2] = 1 - (k - 2)^2 / (k * (k - 1))
  a[1, 3] = 1 + (3 * k - 6) / k
  a[2, 3] = 1 - (3 * k^2 - 15 * k + 19) / (k * (k - 1))
  a[3, 3] = 1 + (k - 3)^3 / (k * (k - 1) * (k - 2))
  a[1, 4] = 1 + (4 * k - 10) / k
  a[2, 4] = 1 - (6 * k^2 - 36 * k + 55) / (k * (k - 1))
  a[3, 4] = 1 + (4 * k^3 - 42 * k^2 + 148 * k - 175) / (k * (k - 1) * (k - 2))
  a[4, 4] = 1 - (k - 4)^4 / (k * (k - 1) * (k - 2) * (k - 3))
  a[1, 5] = 1 + (5 * k - 15) / k
  a[2, 5] = 1 - (10 * k^2 - 70 * k + 125) / (k * (k - 1))
  a[3, 5] = 1 + (10 * k^3 - 120 * k^2 + 485 * k - 660) / (k * (k - 1) * (k - 2))
  a[4, 5] = 1 - ((k - 4)^5 - (k - 5)^5) / (k * (k - 1) * (k - 2) * (k - 3))
  a[5, 5] = 1 + (k - 5)^5 / (k * (k - 1) * (k - 2) * (k - 3) * (k - 4))
  N = sapply(c(1:5), function(x) S + sum(f * (a[, x] - 1)))
  T = sapply(c(1:4), function(x)Tk(N, f, a, S, x))
  p = sapply(c(1:4), function(x)p_value(T[x]))
  se = sapply(c(1:5), function(k)seN(N, f, a, k))
  null_hyp = as.data.table(cbind(T, p))
  null_hyp = rbind(null_hyp, list(NaN, NaN))
  optimal = rep('', 5)
  opt_index = which(null_hyp$p > 0.05)[1]
  optimal[opt_index] = '*'
  if (is.na(opt_index)) {
    opt_index = 5
  }
  lowCI = N - se * 1.96
  uppCI = N + se * 1.96
  return(list(N = N[opt_index], sn = S, se = se[opt_index], lowCI = lowCI[opt_index], uppCI = uppCI[opt_index]))
}

get_chao2 <- function(n_samples, counts) {
  m = n_samples
  f1 = counts[1]
  f2 = counts[2]
  s_obs = sum(counts)
  coeff = (m - 1) / m
  if (f1 > 0 && f2 > 0) {
    est = s_obs + coeff * (f1^2 / (2 * f2))
    var = f2 * (1 / 2 * coeff * (f1 / f2)^2 +
      coeff^2 * (f1 / f2)^3 +
      1 / 4 * coeff^2 * (f1 / f2)^4)
    T = est - s_obs
    K = exp(1.96 * log(1 + var / T^2)^(1 / 2))
    lowCI = s_obs + T / K
    uppCI = s_obs + T * K
  }
  if (f1 > 0 && f2 == 0) {
    est = s_obs + coeff * (f1 * (f1 - 1) / (2 * (f2 + 1)))
    var = 1 / 2 * coeff * (f1 * (f1 - 1) / 2) + coeff^2 * f1 * ((2 * f1 - 1) / 2)^2 - 1 / 4 * coeff^2 * f1^4 / est
    T = est - s_obs
    K = exp(1.96 * log(1 + var / T^2)^(1 / 2))
    lowCI = s_obs + T / K
    uppCI = s_obs + T * K
  }
  if ((f1 == 1 && f2 == 0) ||
    (f1 == 0 && f2 > 0) ||
    (f1 == 0 && f2 == 0)) {
    est = s_obs
    e_1 = sapply(seq_along(counts), function(i) exp(-i) - exp(-2 * i))
    e_2 = sapply(seq_along(counts), function(i) i * exp(-i))
    var = sum(counts * e_1) - 1 / m * sum(counts * e_2)^2
    e_3 = sapply(seq_along(counts), function(i) exp(-i))
    P = sum(counts * e_3) / s_obs
    lowCI = max(s_obs, s_obs / (1 - P) - 1.96 * sqrt(var) / (1 - P))
    uppCI = s_obs / (1 - P) + 1.96 * sqrt(var) / (1 - P)
  }
  return(list(N = est, se = sqrt(var), lowCI = lowCI, uppCI = uppCI))
}

get_chao2_bias_corrected <- function(n_samples, counts) {
  # under "no hetegogenity" assumption
  m = n_samples
  f1 = counts[1]
  f2 = counts[2]
  s_obs = sum(counts)
  coeff = (m - 1) / m
  if (f1 > 0 && f2 > 0) {
    est = s_obs + coeff * (f1 * (f1 - 1) / (2 * (f2 + 1)))
    var = 1 / 2 *
      coeff *
      (f1 * (f1 - 1) / (f2 + 1)) +
      coeff^2 *
        f1 *
        ((2 * f1 - 1) / (f2 + 1))^2 +
      1 / 4 * coeff^2 * f1^2 * f2 * (f1 - 1)^2 / (f2 + 1)^4
    T = est - s_obs
    K = exp(1.96 * log(1 + var / T^2)^(1 / 2))
    lowCI = s_obs + T / K
    uppCI = s_obs + T * K
  }
  if (f1 > 0 && f2 == 0) {
    est = s_obs + coeff * (f1 * (f1 - 1) / (2 * (f2 + 1)))
    var = 1 / 2 * coeff * (f1 * (f1 - 1) / 2) + coeff^2 * f1 * ((2 * f1 - 1) / 2)^2 - 1 / 4 * coeff^2 * f1^4 / est
    T = est - s_obs
    K = exp(1.96 * log(1 + var / T^2)^(1 / 2))
    lowCI = s_obs + T / K
    uppCI = s_obs + T * K
  }
  if ((f1 == 1 && f2 == 0) ||
    (f1 == 0 && f2 > 0) ||
    (f1 == 0 && f2 == 0)) {
    est = s_obs
    e_1 = sapply(seq_along(counts), function(i) exp(-i) - exp(-2 * i))
    e_2 = sapply(seq_along(counts), function(i) i * exp(-i))
    var = sum(counts * e_1) - 1 / m * sum(counts * e_2)^2
    e_3 = sapply(seq_along(counts), function(i) exp(-i))
    P = sum(counts * e_3) / s_obs
    lowCI = max(s_obs, s_obs / (1 - P) - 1.96 * sqrt(var) / (1 - P))
    uppCI = s_obs / (1 - P) + 1.96 * sqrt(var) / (1 - P)
  }
  return(list(N = est, se = sqrt(var), lowCI = lowCI, uppCI = uppCI))
}

get_counts <- function(subject, testsuite, data_dir, catch_type) {
  data_type = 'npz'
  data_file = get_npz_file(data_dir, catch_type, subject, testsuite)
  identity_matrix = load_data(data_file, data_type)
  if (is.null(identity_matrix)) {
    return(NULL)
  }
  counts_obj = load_counts(data_type, identity_matrix)
  counts = counts_obj$counts
  n_samples = counts_obj$n_samples
  Sn = counts_obj$Sn
  total_mutants = counts_obj$total_mutants
  rr = list(subject = subject,
            testsuite = testsuite,
            catch_type = catch_type,
            n_samples = n_samples,
            f1 = counts[1],
            f2 = counts[2])
  return(rr)
}

get_all_counts <- function(subject, testsuite, data_dir, catch_type) {
  data_type = 'npz'
  data_file = get_npz_file(data_dir, catch_type, subject, testsuite)
  identity_matrix = load_data(data_file, data_type)
  if (is.null(identity_matrix)) {
    return(NULL)
  }
  counts_obj = load_counts(data_type, identity_matrix)
  counts = counts_obj$counts
  n_samples = counts_obj$n_samples
  Sn = counts_obj$Sn
  total_mutants = counts_obj$total_mutants
  rr = list(subject = subject,
            testsuite = testsuite,
            catch_type = catch_type,
            n_samples = n_samples,
            counts = counts)
  return(rr)
}

print_all_counts <- function(subjects, testsuites, catch_types, input_dir, results_dir) {
  for (testsuite in testsuites) {
    for (subject in subjects) {
      for (catch_type in catch_types) { # all catch-types shoudl produce the same # of mutants
        data = get_all_counts(subject, testsuite, input_dir, catch_type)
        if (is.null(data))
          next
        counts = data$counts
        counts = counts[0:max(which(counts != 0))]
        fname = results_dir %//% glue("{subject}_{testsuite}_{catch_type}.txt")
        fwrite(list(counts), file = as.character(fname))
      } }
  }
}

print_counts <- function(subjects, testsuites, catch_types, input_dir, results_file) {
  all_data = list()
  for (testsuite in testsuites) {
    for (subject in subjects) {
      for (catch_type in catch_types) {
        data = get_counts(subject, testsuite, input_dir, catch_type)
        all_data = append(all_data, list(data))
      } } }
  all_data = as.data.table(do.call("rbind", all_data))
  fwrite(all_data, file = results_file)
}

get_total_mutants <- function(killed_mutants_path) {
  if (is.na(killed_mutants_path))
    return(NA)
  return(read.csv(file = killed_mutants_path))
}

estimate_all <- function(subjects, testsuites, estimators, catch_types, ci, input_dir, results_file, data_type, killed_mutants_path, log_dir, time_limit = 3600) {
  dir.create(log_dir, showWarnings = FALSE)
  cl <- parallel::makeCluster(cores, outfile = "", type = "FORK")
  doParallel::registerDoParallel(cl)
  total_mutants_df = get_total_mutants(killed_mutants_path)
  all_data =
    foreach(testsuite = testsuites, .combine = "rbind") %:%
      foreach(subject = subjects, .combine = "rbind") %:%
      foreach(catch_type = catch_types, .combine = "rbind", .packages = c("data.table", "pathlibr", "glue", "SPECIES"),
              .export = c("compute_estimators", "get_npz_file", "load_data", "load_counts", "compute_estimator", "get_est_type", "counts_estimators", "matrix_estimators", "estimate_from_counts", "get_abundance_file",
                          "estimators_fun", "estimate_from_counts", "item_or_default", "estimators_additional_param", "run_estimator", "estimator_bootstrap", "estimator_zelterman", "estimate_from_matrix", "get_additional_params",
                          "estimator_jackknife", "estimator_pcg", "estimator_ChaoBunge", "estimator_unpmle", "estimator_pnpmle")) %dopar% {
      withCallingHandlers({
        setTimeLimit(time_limit, transient = TRUE)
        compute_estimators(subject, testsuite, estimators, input_dir, catch_type, ci, data_type, log_dir, total_mutants_df)
      },
        error = function(e) {
          message = return(data.frame(subject = subject, testsuite = testsuite, type = catch_type, total_mutants = -1,
                                      killed = -1, estimation = -1, CILower = -1, CIUpper = 1))
          # timeout hit
          # do stuff to capture error messages here
        }
      )
    }
  parallel::stopCluster(cl)
  fwrite(as.data.table(all_data), file = results_file)
}

if (sys.nframe() == 0 || TRUE) {
  parser = ArgumentParser()
  parser$add_argument("-i", "--input-dir", help = "directory with data files")
  parser$add_argument("-s", "--subjects", nargs = '+', choices = all_subjects, default = default_subjects, help = "subjects list, divided with space")
  parser$add_argument("-t", "--testsuites", nargs = '+', choices = all_testsuites, default = default_testsuites, help = "testsuite to draw")
  parser$add_argument("-e", "--estimators", nargs = '+', choices = all_estimators, default = default_estimators, help = "estimators to calculate")
  parser$add_argument("-c", "--catch-types", nargs = '+', choices = all_catch_types, default = default_catch_types, help = "mutation catch type")
  parser$add_argument("-o", "--output", help = "Filename for the resulting comparison table file.csv")
  parser$add_argument("--ci", action = 'store_true', help = "Compute confidence intervals for estimators from species package, is slow")
  parser$add_argument("--data-type", choices = c('npz', 'abundance'), default = 'npz', help = "deprecated, use npz")
  parser$add_argument("--counts", choices = c('all', 'f1f2', 'none'), default = 'none', help = "don't estimate, just print counts")
  parser$add_argument("--log", default = './logs', help = "logs folder")
  parser$add_argument("--total", help = "csv file with total mutants data")
  parser$add_argument("--cores", default = 10, help = "cores")
  args = parser$parse_args()
  cores = args$cores
  if (args$counts == 'f1f2') {
    print_counts(args$subjects, args$testsuites, args$catch_types, args$input_dir, args$output)
  } else if (args$counts == 'all') {
    print_all_counts(args$subjects, args$testsuites, args$catch_types, args$input_dir, args$output)
  } else {
    estimate_all(args$subjects, args$testsuites, args$estimators, args$catch_types, args$ci, args$input_dir, args$output, args$data_type, args$total, args$log)
  }
}
