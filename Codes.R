rm(list = ls())
pkgs <- c("pbmcapply", "foreign", "withr", "DescTools")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0)
  install.packages(to_install, dependencies = TRUE)
if (!requireNamespace("parTreat", quietly = TRUE)) {
  install.packages("remotes")
  remotes::install_github("michaelpollmann/parTreat",
                          dependencies = FALSE,    # avoid pulling heavy deps
                          upgrade = "never",
                          build_vignettes = FALSE)
}
library(pbmcapply)
library(withr)
library(DescTools)
library(parTreat)
n.cores <- parallel::detectCores()

trimmed_mean <- function(x, trim = 0.1) {
  n <- length(x); k <- floor(trim * n); if (k <= 0L) return(mean(x))
  mean(sort(x)[(k+1):(n-k)])
}

winsorized_mean <- function(x, wins = 0.1) {
  n <- length(x); k <- floor(wins * n); if (k <= 0) return(mean(x))
  xs <- sort(x)
  xs[seq_len(k)] <- xs[k+1]; xs[(n-k+1):n] <- xs[n-k]
  mean(xs)
}

stat_unadj <- function(y, z, trim = 0.1, wins = 0.1) {
  y1 <- y[z==1]; y0 <- y[z==0]
  c(
    DM = mean(y1) - mean(y0),
    DiffMed = stats::median(y1) - stats::median(y0),
    HL      = DescTools::HodgesLehmann(y1, y0),
    TrimDM  = trimmed_mean(y1, trim) - trimmed_mean(y0, trim),
    WinsDM  = winsorized_mean(y1, wins) - winsorized_mean(y0, wins)
  )
}

stat_adj <- function(y, z, x) {
  n <- length(y)
  if (length(z) != n) stop("y and z must have equal length.")
  X <- if (is.null(dim(x))) data.frame(x = x) else as.data.frame(x)
  if (nrow(X) != n) stop("nrow(x) must equal length(y).")
  
  # coerce character -> factor; leave existing factors as is
  X[] <- lapply(X, function(col) if (is.character(col)) factor(col) else col)
  
  # split controls
  is_num <- vapply(X, is.numeric, TRUE)
  num_names <- names(X)[is_num]
  fac_names <- names(X)[!is_num]
  
  # drop NAs
  df_un <- cbind.data.frame(y = y, z = z, X)
  df_un <- df_un[stats::complete.cases(df_un), , drop = FALSE]
  N <- nrow(df_un); m <- sum(df_un$z == 1L)
  
  # ---------- tau.reg: lm(y ~ z + numeric + as.factor(factors)) ----------
  num_terms_reg <- if (length(num_names)) paste(num_names, collapse = " + ") else ""
  fac_terms_reg <- if (length(fac_names)) paste0("as.factor(", fac_names, ")", collapse = " + ") else ""
  rhs_reg <- paste(Filter(nzchar, c("z", num_terms_reg, fac_terms_reg)), collapse = " + ")
  fit1 <- stats::lm(stats::as.formula(paste0("y ~ ", rhs_reg)), data = df_un)
  tau.reg <- unname(coef(fit1)["z"])
  
  # ---------- tau.interact: center ONLY numeric X, interact z with centered numerics; add factors as as.factor(...) ----------
  if (length(num_names)) {
    Xc <- X
    Xc[num_names] <- lapply(Xc[num_names], function(col) col - mean(col, na.rm = TRUE))
    cn <- paste0(num_names, "_c")
    names(Xc)[match(num_names, names(Xc))] <- cn
    df_int <- cbind.data.frame(y = y, z = z, Xc)
    
    num_terms_int <- paste(cn, collapse = " + ")
    fac_terms_int <- if (length(fac_names)) paste0(" + ", paste0("as.factor(", fac_names, ")", collapse = " + ")) else ""
    f_int <- stats::as.formula(paste0("y ~ z * (", num_terms_int, ")", fac_terms_int))
    fit2 <- stats::lm(f_int, data = df_int)
    tau.interact <- unname(coef(fit2)["z"])
  } else {
    # no numeric covariates to interact with
    rhs_int <- paste(Filter(nzchar, c("z", fac_terms_reg)), collapse = " + ")
    fit2 <- stats::lm(stats::as.formula(paste0("y ~ ", rhs_int)), data = df_un)
    tau.interact <- unname(coef(fit2)["z"])
  }
  
  # ---------- tau.Radj: root of rank balance after adjusting (y - tau z) on UNcentered X with as.factor for factors ----------
  ctrl_rhs <- paste(
    Filter(nzchar, c(num_terms_reg, fac_terms_reg)),
    collapse = " + "
  )
  aux <- function(tau) {
    y_shift <- df_un$y - tau * df_un$z
    fit <- stats::lm(stats::as.formula(
      if (nzchar(ctrl_rhs)) paste0("y_shift ~ ", ctrl_rhs) else "y_shift ~ 1"
    ), data = transform(df_un, y_shift = y_shift))
    r <- resid(fit)
    sum(rank(r)[df_un$z == 1L]) - m * (N + 1) / 2
  }
  tau.Radj <- uniroot(aux, c(0, 1), extendInt = "yes")$root
  
  c(tau.reg = tau.reg, tau.interact = tau.interact, tau.Radj = tau.Radj)
}

# ------------------------------------------------
# Real data example: Progresa data
# ------------------------------------------------
set.seed(42)
Progressa = read.csv("./ag/PROGRESA.csv")
Y = Progressa$pri2000s
Z = Progressa$treatment
X = with(Progressa, data.frame(villages, pri1994, pan1994, prd1994,
                               votos1994, avgpoverty, pobtot1994))
X$villages = as.factor(X$villages)
B <- 1000; alpha <- 0.05; zc <- qnorm(1 - alpha/2)

# EIF and WAQ from Athey et al. (2023)
y0 <- Y[Z == 0]; y1 <- Y[Z == 1]
eif <- as.numeric(parTreat::eif_additive(y0, y1)$tau)
waq <- as.numeric(parTreat::waq(y0, y1)$tau)

That <- c(stat_unadj(Y, Z), stat_adj(Y, Z, X), eif, waq)

n  <- length(Z); n1 <- sum(Z == 1)
Tperm <- do.call(rbind, 
                 pbmclapply(1:B, 
                            function(b){
                              set.seed(b)
                              zpi <- integer(n); zpi[sample.int(n, n1, replace = FALSE)] <- 1
                              # EIF and WAQ from Athey et al. (2023)
                              y0.zpi <- Y[zpi == 0]; y1.zpi <- Y[zpi == 1]
                              eif.zpi <- as.numeric(parTreat::eif_additive(y0.zpi, y1.zpi)$tau)
                              waq.zpi <- as.numeric(parTreat::waq(y0.zpi, y1.zpi)$tau)
                              c(stat_unadj(Y, zpi), stat_adj(Y, zpi, X), eif.zpi, waq.zpi)
                            },
                            mc.cores = n.cores
                 )
)
dim(Tperm)

# Two-sided (1 - alpha) CI via permutation-based critical values
c_abs <- apply(abs(na.omit(Tperm)), 
               2, 
               quantile, 
               probs = 1 - alpha, 
               names = FALSE)
CI_lo <- That - c_abs; CI_hi <- That + c_abs

# output
out <- data.frame(
  methods = c('Difference-in-Means',
              'Difference-in-Medians', 
              'Rosenbaum unadj', 
              '0.1-trimmed-Difference-in-Means', 
              '0.1-Winsorized-Difference-in-Means', 
              'OLS adjusted', 
              'Lin estimator', 
              'Rosenbaum adj',
              'eif (Athey et al., 2023)', 
              'waq (Athey et al., 2023)'),
  estimate = as.numeric(That),
  std.err = as.numeric(CI_hi - CI_lo)/(2*zc),
  ci_lo = as.numeric(CI_lo),
  ci_hi = as.numeric(CI_hi),
  row.names = NULL
)

print(out[c(1:2,4:5,9:10,3,6:8),], digits = 4)

# Following are the results for B = 1e5
# methods                          estimate std.err ci_lo  ci_hi
# Difference-in-Means                3.622   1.922 -0.1446 7.390
# Difference-in-Medians              0.690   1.562 -2.3720 3.752
# 0.1-trimmed-Difference-in-Means    1.998   1.676 -1.2866 5.282
# 0.1-Winsorized-Difference-in-Means 2.590   1.718 -0.7775 5.958
# eif (Athey et al., 2023)           1.953   1.721 -1.4203 5.326
# waq (Athey et al., 2023)           1.306   1.708 -2.0418 4.653
# Rosenbaum unadj                    1.834   1.666 -1.4322 5.100
# OLS adjusted                       3.671   1.701  0.3378 7.004
# Lin estimator                      4.214   1.985  0.3224 8.105
# Rosenbaum adj                      2.185   1.384 -0.5289 4.898


# ------------------------------------------------
# Synthetic experiments
# ------------------------------------------------

# --- single replication ---
sim_once <- function(setting = "1a", N = 1e3, p = 0.5, 
                     tau0 = 2, alpha = 0.05, B = 1e3, seed = 42) {
  m = round(p * N); zc  <- qnorm(1 - alpha/2)
  withr::with_seed(seed, {
    if (setting == "1a"){
      X <- runif(N, min = -4, max = 4)
      a <- 3*X + rnorm(N)
    } else if (setting == "1b"){
      X <- runif(N, min = -4, max = 4)
      a <- 3*X + rt(N, df = 1)
    } else if (setting == "1c"){
      X <- runif(N, min = -4, max = 4)
      a <- 3*X + rt(N, df = 3)
    } else if (setting == "2a"){
      X <- runif(N, min = -4, max = 4)
      a <- 3*X + rnorm(N)
    } else if (setting == "2b"){
      X <- runif(N, min = -4, max = 4)
      a <- 3*X + rt(N, df = 1)
    } else if (setting == "2c"){
      X <- runif(N, min = -4, max = 4)
      a <- 3*X + rt(N, df = 3)
    } else if(setting == "3a"){
      X <- runif(N, min = -4, max = 4)
      a <- rexp(N, rate = 1/10) + rnorm(N)
    } else if (setting == "3b"){
      X <- runif(N, min = -4, max = 4)
      a <- rexp(N, rate = 1/10) + rt(N, df = 1)
    } else if (setting == "3c"){
      X <- runif(N, min = -4, max = 4)
      a <- rexp(N, rate = 1/10) + rt(N, df = 3)
    } else if (setting == "4a"){
      X <- exp(runif(N, min = -4, max = 4))
      a <- (X + sqrt(X))/4 + rnorm(N)
    } else if (setting == "4b"){
      X <- exp(runif(N, min = -4, max = 4))
      a <- (X + sqrt(X))/4 + rt(N, df = 1)
    } else if (setting == "4c"){
      X <- exp(runif(N, min = -4, max = 4))
      a <- (X + sqrt(X))/4 + rt(N, df = 3)
    }
    b <- a - tau0
    Z <- as.numeric(sample.int(N) <= m)
    Y <- a * Z + b * (1 - Z)
    if (setting %in% c("2a", "2b", "2c")){
      contam.indices = c(which(Z==1)[sample(m, round(0.05*m), replace = F)],
                         which(Z==0)[sample(N-m, round(0.05*(N-m)), replace = F)])
      Y[contam.indices] = 500
    }
    
    # EIF and WAQ from Athey et al. (2023)
    y0 <- Y[Z == 0]; y1 <- Y[Z == 1]
    eif <- as.numeric(parTreat::eif_additive(y0, y1)$tau)
    waq <- as.numeric(parTreat::waq(y0, y1)$tau)
    
    That <- c(stat_unadj(Y, Z), stat_adj(Y, Z, X), eif, waq)
    
    n  <- length(Z); n1 <- sum(Z == 1)
    
    Tperm <- replicate(B, {
      zpi <- integer(n); zpi[sample.int(n, n1, replace = FALSE)] <- 1
      # EIF and WAQ from Athey et al. (2023)
      y0.zpi <- Y[zpi == 0]; y1.zpi <- Y[zpi == 1]
      eif.zpi <- as.numeric(parTreat::eif_additive(y0.zpi, y1.zpi)$tau)
      waq.zpi <- as.numeric(parTreat::waq(y0.zpi, y1.zpi)$tau)
      c(stat_unadj(Y, zpi), stat_adj(Y, zpi, X), eif.zpi, waq.zpi)
    })
    
    c_abs <- apply(abs(Tperm),
                   1,
                   quantile,
                   probs = 1 - alpha,
                   na.rm = TRUE,
                   names = FALSE)
  })
  data.frame(
    method = c('Difference-in-Means',
               'Difference-in-Medians', 
               'Rosenbaum unadj', 
               '0.1-trimmed-Difference-in-Means', 
               '0.1-Winsorized-Difference-in-Means', 
               'OLS adjusted', 
               'Lin estimator', 
               'Rosenbaum adj',
               'eif (Athey et al., 2023)', 
               'waq (Athey et al., 2023)'),
    est    = as.numeric(That),
    lo     = as.numeric(That - c_abs),
    hi     = as.numeric(That + c_abs),
    width  = as.numeric(2 * c_abs),
    row.names = NULL
  )
}

# --- run many reps in parallel, summarize ---
run_permutation_study <- function(R = 100,
                                  setting = "1a", N = 1e3, p = 0.5, 
                                  tau0 = 2, alpha = 0.05, B = 1e3,
                                  cores = n.cores) {
  res_list <- pbmcapply::pbmclapply(
    seq_len(R),
    function(i){
      sim_once(setting, N, p, tau0, alpha, B, i)
      },
    mc.cores = cores
  )
  
  methods <- res_list[[1]]$method
  ci_arr  <- array(NA_real_, dim = c(R, length(methods), 2),
                   dimnames = list(NULL, methods, c("lo", "hi")))
  wid_mat <- matrix(NA_real_, nrow = R, ncol = length(methods),
                    dimnames = list(NULL, methods))
  
  for (i in seq_len(R)) {
    ci_arr[i, , "lo"] <- res_list[[i]]$lo
    ci_arr[i, , "hi"] <- res_list[[i]]$hi
    wid_mat[i, ]  <- res_list[[i]]$width
  }
  
  na_methods <- names(which(colSums(is.na(ci_arr[, , "lo"]) | is.na(ci_arr[, , "hi"])) > 0))
  if (length(na_methods)) {
    message(paste(na_methods, collapse = ", "), " produced NA estimate in some replications")
  }
  na_width_methods <- names(which(colSums(is.na(wid_mat)) > 0))
  if (length(na_width_methods)) {
    message(paste(na_width_methods, collapse = ", "), " produced NA se in some replications")
  }
  coverage  <- colMeans(ci_arr[, , "lo"] <= tau0 & tau0 <= ci_arr[, , "hi"], na.rm = TRUE)
  avg_width <- colMeans(wid_mat, na.rm = TRUE)
  median_width<- apply(wid_mat, 2, median, na.rm = TRUE)
  
  summary_mat <- cbind(coverage = coverage,
                       avg_width = avg_width,
                       median_width = median_width)
  
  list(ci_array = ci_arr, width_matrix = wid_mat, 
       reps = res_list, summary = summary_mat)
}

R <- 1000; B = 1000; N <- 1000; tau0 = 2; alpha = 0.05
set.seed(42)
p = 0.25
all_out <- list()
for (setting in paste0(rep(1:4, each = 3), letters[1:3])) {
  cat("Setting", setting, "\n")
  ans <- run_permutation_study(setting = setting, R = R, N = N, p = p, 
                               tau0 = tau0, alpha = alpha, B = B)
  print(round(ans$summary[c(1:2,4:5,9:10,3,6:8),], 3))
  
  tmp <- as.data.frame(ans$summary[c(1:2,4:5,9:10,3,6:8),])
  tmp$method  <- rownames(tmp)
  tmp$setting <- setting
  all_out[[setting]] <- tmp[, c("setting","method","coverage","avg_width","median_width")]
}

final_df <- do.call(rbind, all_out)
write.csv(final_df, 
          paste0("./ag/perm_study_summary_p",p,"_R",R,"_B",B,".csv"), 
          row.names = FALSE)

df <- read.csv(paste0("./ag/perm_study_summary_p",p,"_R",R,"_B",B,".csv"), stringsAsFactors = FALSE)
for (s in unique(df$setting)) {
  cat("Setting", s, "\n")
  tab <- df[df$setting == s, c("method","coverage","avg_width","median_width")]
  rownames(tab) <- tab$method
  print(round(tab[, c("coverage","avg_width","median_width")], 3))
  cat("\n")
}

set.seed(42)
p = 0.5
all_out <- list()
for (setting in paste0(rep(1:4, each = 3), letters[1:3])) {
  cat("Setting", setting, "\n")
  ans <- run_permutation_study(setting = setting, R = R, N = N, p = p, 
                               tau0 = tau0, alpha = alpha, B = B)
  print(round(ans$summary[c(1:2,4:5,9:10,3,6:8),], 3))
  
  tmp <- as.data.frame(ans$summary[c(1:2,4:5,9:10,3,6:8),])
  tmp$method  <- rownames(tmp)
  tmp$setting <- setting
  all_out[[setting]] <- tmp[, c("setting","method","coverage","avg_width","median_width")]
}

final_df <- do.call(rbind, all_out)
write.csv(final_df, 
          paste0("./ag/perm_study_summary_p",p,"_R",R,"_B",B,".csv"), 
          row.names = FALSE)

df <- read.csv(paste0("./ag/perm_study_summary_p",p,"_R",R,"_B",B,".csv"), stringsAsFactors = FALSE)
for (s in unique(df$setting)) {
  cat("Setting", s, "\n")
  tab <- df[df$setting == s, c("method","coverage","avg_width","median_width")]
  rownames(tab) <- tab$method
  print(round(tab[, c("coverage","avg_width","median_width")], 3))
  cat("\n")
}