# =============================================================================
# Title:
# GSE256335 Exact-Matched Synthetic Benchmark with Corrected Onset Semantics
#
# Purpose:
#   Generate and test a synthetic benchmark matched to the real mESC experiment
#   GSE256335, using the corrected ode.r semantics:
#
#       - gene-specific transcriptional onset;
#       - one common pharmacological shutoff time;
#       - no post-hoc observation-time shift.
#
# Exact matched design:
#       relative times after shutoff = 0, 30, 60, 120, 240 min
#       biological replicates        = 3
#       platform                     = RNA-seq
#       perturbation                 = complete SHUTOFF
#
# For each simulated gene:
#       1. sample kinetic parameters once;
#       2. sample ONE transcriptional onset shift once;
#       3. simulate a sufficiently long pre-shutoff history;
#       4. apply the common shutoff;
#       5. sample exactly the GSE256335 time points;
#       6. apply RNA-seq observation noise;
#       7. run test_sigma_nested().
#
# Null and alternative genes are simulated so the same run provides:
#       - empirical Type-I calibration;
#       - power;
#       - sigma_c recovery;
#       - numerical diagnostics.
#
# Outputs:
#       GSE256335_matched_corrected_onset_raw.tsv
#       GSE256335_matched_corrected_onset_summary.tsv
#       GSE256335_matched_corrected_onset.rdata
#       GSE256335_matched_corrected_onset_progress.tsv
#       GSE256335_matched_corrected_onset_sessionInfo.txt
#
# =============================================================================


# =============================================================================
# 0. Environment
# =============================================================================

setwd("~/postexport-kinetics/synthetic_dataset")

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

library(data.table)
library(parallel)
library(nnls)
library(MASS)
library(deSolve)

data.table::setDTthreads(1)

source("../ode_model/ode.r")
source("../commons/nested_test2.r")
source("../commons/platforms.r")


# =============================================================================
# 1. Settings
# =============================================================================

N_NULL <- 1000L
N_ALT  <- 1000L

N_BOOT <- 1999L

N_WORKERS <- 100L
TASK_BATCH_SIZE <- 100L

PSOCK_TIMEOUT_SECONDS <- 12L * 60L * 60L

BENCHMARK_VERSION <- "GSE256335_corrected_onset_v1"

RELATIVE_TIMES <- c(
  0,
  30,
  60,
  120,
  240
)

N_REPLICATES <- 3L

# Absolute latent time of drug addition.
# This is arbitrary internally; only relative post-shutoff times are analyzed.
T_STAR_ABSOLUTE <- 300

# Gene-specific onset heterogeneity.
MAX_ONSET_SHIFT <- 54L

POST_R_FRACTION <- 0

RANGE_NOISE <- c(
  "Very low",
  "Low",
  "Medium",
  "High"
)

SCALING_A <- TRUE
LAMBDA_TIME <- 0.5
LAMBDA_DIAG <- 0.1
REL_FLOOR <- 1e-8
TRUNCATE_NONNEGATIVE_BOOT <- FALSE
MAX_FAILURE_RATE <- 0.05

T_STAR_FIT <- 0

CHECKPOINT_DIR <- "GSE256335_matched_corrected_onset_checkpoints"
PROGRESS_FILE <- "GSE256335_matched_corrected_onset_progress.tsv"
RAW_FILE <- "GSE256335_matched_corrected_onset_raw.tsv"
SUMMARY_FILE <- "GSE256335_matched_corrected_onset_summary.tsv"
RDATA_FILE <- "GSE256335_matched_corrected_onset.rdata"
SESSION_FILE <- "GSE256335_matched_corrected_onset_sessionInfo.txt"

EXPECTED_TESTS_PER_GENE <- length(RANGE_NOISE)


# =============================================================================
# 2. Check corrected ode.r API
# =============================================================================

required_ode_objects <- c(
  "random_params",
  "generate_ODE_states",
  "rna_kinetics"
)

missing_ode <- required_ode_objects[
  !vapply(
    required_ode_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (
  length(missing_ode) > 0L
) {
  stop(
    paste(
      "Missing corrected ode.r objects:",
      paste(
        missing_ode,
        collapse = ", "
      )
    )
  )
}

ode_formals <- names(
  formals(
    generate_ODE_states
  )
)

required_ode_args <- c(
  "use_onset_shift",
  "nominal_onset_time",
  "post_R_fraction"
)

missing_ode_args <- setdiff(
  required_ode_args,
  ode_formals
)

if (
  length(missing_ode_args) > 0L
) {
  stop(
    paste0(
      "generate_ODE_states() does not look like the corrected ode.r. ",
      "Missing arguments: ",
      paste(
        missing_ode_args,
        collapse = ", "
      )
    )
  )
}


# =============================================================================
# 3. Output/checkpoint directory
# =============================================================================

if (!dir.exists(CHECKPOINT_DIR)) {
  dir.create(
    CHECKPOINT_DIR,
    recursive = TRUE
  )
}

CHECKPOINT_DIR <- normalizePath(
  CHECKPOINT_DIR,
  mustWork = TRUE
)


# =============================================================================
# 4. Stable seeds
# =============================================================================

stable_seed_from_string <- function(
  x
) {

  ints <- utf8ToInt(
    enc2utf8(x)
  )

  h <- 104729

  for (ii in ints) {
    h <- (
      h * 1000003 +
        ii
    ) %%
      2000000000
  }

  out <- as.integer(h)

  if (
    !is.finite(out) ||
      out <= 0L
  ) {
    out <- 1L
  }

  out
}


make_seed <- function(
  gene_i,
  truth,
  noise = "none",
  stream
) {

  stable_seed_from_string(
    paste(
      BENCHMARK_VERSION,
      gene_i,
      truth,
      noise,
      stream,
      sep = "|"
    )
  )
}


# =============================================================================
# 5. Observation model
# =============================================================================

add_rnaseq_noise <- function(
  dt,
  noise
) {

  scale_counts <- switch(
    noise,
    "Very low" = 10000,
    "Low"      = 5000,
    "Medium"   = 1000,
    "High"     = 200,
    stop(
      paste(
        "Unknown noise regime:",
        noise
      )
    )
  )

  mean_disp <- switch(
    noise,
    "Very low" = 0.01,
    "Low"      = 0.05,
    "Medium"   = 0.10,
    "High"     = 0.25,
    stop(
      paste(
        "Unknown noise regime:",
        noise
      )
    )
  ) * 0.25

  cv_disp <- switch(
    noise,
    "Very low" = 0.5,
    "Low"      = 0.7,
    "Medium"   = 0.8,
    "High"     = 1.0,
    stop(
      paste(
        "Unknown noise regime:",
        noise
      )
    )
  )

  simulate_rnaseq(
    data.table::copy(
      as.data.table(dt)
    ),
    targets = c(
      "N",
      "C",
      "C_s",
      "N_s"
    ),
    scale_counts = scale_counts,
    mean_disp = mean_disp,
    cv_disp = cv_disp
  )
}


# =============================================================================
# 6. Safe helpers
# =============================================================================

safe_value <- function(
  x,
  field,
  default = NA_real_
) {

  if (
    is.null(x) ||
      is.null(x[[field]]) ||
      length(x[[field]]) == 0L
  ) {
    return(default)
  }

  suppressWarnings(
    as.numeric(
      x[[field]][1]
    )
  )
}


safe_text <- function(
  x,
  field,
  default = NA_character_
) {

  if (
    is.null(x) ||
      is.null(x[[field]]) ||
      length(x[[field]]) == 0L
  ) {
    return(default)
  }

  as.character(
    x[[field]][1]
  )
}


safe_coef <- function(
  x,
  name
) {

  if (
    is.null(x) ||
      is.null(names(x)) ||
      !(name %in% names(x))
  ) {
    return(NA_real_)
  }

  unname(
    x[name]
  )
}


# =============================================================================
# 7. Base kinetic parameters
# =============================================================================

make_base_parameters <- function(
  gene_i,
  truth
) {

  set.seed(
    make_seed(
      gene_i = gene_i,
      truth = truth,
      stream = "kinetics"
    )
  )

  p <- random_params()

  p$R <- rgamma(
    1,
    shape = 5,
    scale = 20
  )

  upper_alpha_s <- min(
    p$alpha,
    r_alpha_s_max
  )

  if (
    upper_alpha_s <=
      r_alpha_s_min
  ) {

    p$alpha_s <- upper_alpha_s

  } else {

    p$alpha_s <- runif(
      1,
      r_alpha_s_min,
      upper_alpha_s
    )
  }

  if (
    truth == "null"
  ) {

    p$sigma_c <- 0

  } else {

    # Alternative sigma_c is drawn from the same range used by the simulator.
    p$sigma_c <- runif(
      1,
      r_sigma_c_min,
      r_sigma_c_max
    )
  }

  p
}


# =============================================================================
# 8. One onset per synthetic gene
# =============================================================================

make_onset_time <- function(
  gene_i,
  truth
) {

  set.seed(
    make_seed(
      gene_i = gene_i,
      truth = truth,
      stream = "onset"
    )
  )

  shift <- sample(
    -MAX_ONSET_SHIFT:
      MAX_ONSET_SHIFT,
    1L
  )

  list(
    onset_shift = shift,
    onset_time = shift
  )
}


# =============================================================================
# 9. Simulate one exact GSE256335 latent gene
# =============================================================================

simulate_latent_gene <- function(
  gene_i,
  truth
) {

  p <- make_base_parameters(
    gene_i = gene_i,
    truth = truth
  )

  onset <- make_onset_time(
    gene_i = gene_i,
    truth = truth
  )

  abs_sample_times <-
    T_STAR_ABSOLUTE +
    RELATIVE_TIMES

  sim_times <- seq(
    0,
    max(
      abs_sample_times
    ),
    by = 1
  )

  y0 <- c(
    N = 0,
    N_s = 0,
    C = 0,
    C_s = 0
  )

  dd <- generate_ODE_states(
    base_params = p,
    y0 = y0,
    times = sim_times,
    n_replicates = N_REPLICATES,
    model_kinetics = rna_kinetics,
    stimes = abs_sample_times,
    shutofftimes = abs_sample_times,
    t_star = T_STAR_ABSOLUTE,
    post_R_fraction = POST_R_FRACTION,

    # IMPORTANT:
    # one onset is generated above and passed explicitly.
    use_onset_shift = FALSE,
    max_shift = 0,
    nominal_onset_time = onset$onset_time
  )

  ts <- as.data.table(
    dd$intervention_tsampled_data
  )

  ts[
    ,
    time :=
      time -
      T_STAR_ABSOLUTE
  ]

  setorder(
    ts,
    time,
    replicate
  )

  if (
    !identical(
      as.numeric(
        sort(
          unique(
            ts$time
          )
        )
      ),
      as.numeric(
        RELATIVE_TIMES
      )
    )
  ) {
    stop(
      paste(
        "Wrong matched time points for gene",
        gene_i,
        truth
      )
    )
  }

  rep_check <- ts[
    ,
    uniqueN(
      replicate
    ),
    by = time
  ]

  if (
    nrow(rep_check) !=
      length(
        RELATIVE_TIMES
      ) ||
      any(
        rep_check$V1 !=
          N_REPLICATES
      )
  ) {
    stop(
      paste(
        "Wrong replicate structure for gene",
        gene_i,
        truth
      )
    )
  }

  # The corrected semantics require the intervention to be exact.
  post_check <- as.data.table(
    dd$intervention_data
  )[
    time >= T_STAR_ABSOLUTE
  ]

  if (
    any(
      abs(
        post_check$R
      ) >
        1e-12
    )
  ) {
    stop(
      paste(
        "Residual R found after complete shutoff for gene",
        gene_i,
        truth
      )
    )
  }

  list(
    data = ts,
    params = p,
    onset_shift = onset$onset_shift,
    onset_time = onset$onset_time
  )
}


# =============================================================================
# 10. Run one gene over all RNA-seq noise regimes
# =============================================================================

run_gene <- function(
  gene_i,
  truth,
  N_boot
) {

  latent <- simulate_latent_gene(
    gene_i = gene_i,
    truth = truth
  )

  out <- vector(
    "list",
    length(
      RANGE_NOISE
    )
  )

  for (
    kk in seq_along(
      RANGE_NOISE
    )
  ) {

    noise_i <- RANGE_NOISE[kk]

    set.seed(
      make_seed(
        gene_i = gene_i,
        truth = truth,
        noise = noise_i,
        stream = "measurement"
      )
    )

    noisy <- add_rnaseq_noise(
      latent$data,
      noise_i
    )

    test_seed <- make_seed(
      gene_i = gene_i,
      truth = truth,
      noise = noise_i,
      stream = "bootstrap"
    )

    t0 <- Sys.time()

    fit <- tryCatch(
      test_sigma_nested(
        tsampled_data = noisy,
        scaling_A = SCALING_A,
        t_star = T_STAR_FIT,
        B_n = N_boot,
        seed = test_seed,
        lambda_time = LAMBDA_TIME,
        lambda_diag = LAMBDA_DIAG,
        rel_floor = REL_FLOOR,
        truncate_nonnegative_boot =
          TRUNCATE_NONNEGATIVE_BOOT,
        max_failure_rate =
          MAX_FAILURE_RATE,
        return_boot = FALSE,
        verbose = FALSE
      ),
      error = function(e) {
        list(
          status = "test_error",
          error.message =
            conditionMessage(e),
          p.value = NA_real_
        )
      }
    )

    coef_full <- fit$coef_full
    coef_null <- fit$coef_null

    out[[kk]] <- data.table(
      Benchmark_version = BENCHMARK_VERSION,

      Gene = gene_i,

      truth = truth,

      Positive =
        as.integer(
          truth == "alternative"
        ),

      sigma_true =
        latent$params$sigma_c,

      Platform =
        "RNA-seq",

      Perturbation =
        "SHUTOFF",

      Exprs_noise =
        noise_i,

      N_tsamples =
        length(
          RELATIVE_TIMES
        ),

      N_replicates =
        N_REPLICATES,

      sample_times =
        paste(
          RELATIVE_TIMES,
          collapse = ","
        ),

      onset_shift =
        latent$onset_shift,

      onset_time =
        latent$onset_time,

      T_star =
        T_STAR_FIT,

      Post_R_fraction =
        POST_R_FRACTION,

      p.value =
        safe_value(
          fit,
          "p.value"
        ),

      status =
        safe_text(
          fit,
          "status"
        ),

      error_message =
        safe_text(
          fit,
          "error.message"
        ),

      T.obs =
        safe_value(
          fit,
          "T.obs"
        ),

      Sigma =
        safe_value(
          fit,
          "Sigma"
        ),

      IR =
        safe_value(
          fit,
          "IR"
        ),

      atom_zero =
        safe_value(
          fit,
          "atom.zero"
        ),

      condition_number =
        safe_value(
          fit,
          "condition.number"
        ),

      min_singular_value =
        safe_value(
          fit,
          "min.singular.value"
        ),

      rank_full =
        safe_value(
          fit,
          "rank.full"
        ),

      rank_null =
        safe_value(
          fit,
          "rank.null"
        ),

      bootstrap_failure_rate =
        safe_value(
          fit,
          "bootstrap.failure.rate"
        ),

      boot_rank_deficient_fraction =
        safe_value(
          fit,
          "bootstrap.rank.deficient.fraction"
        ),

      sigma_c_hat =
        safe_coef(
          coef_full,
          "sigma_c"
        ),

      sigma_n_hat =
        safe_coef(
          coef_full,
          "sigma_n"
        ),

      tau_hat =
        safe_coef(
          coef_full,
          "tau"
        ),

      tau_s_hat =
        safe_coef(
          coef_full,
          "tau_s"
        ),

      alpha_hat =
        safe_coef(
          coef_full,
          "alpha"
        ),

      alpha_s_hat =
        safe_coef(
          coef_full,
          "alpha_s"
        ),

      sigma_n_hat_null =
        safe_coef(
          coef_null,
          "sigma_n"
        ),

      tau_hat_null =
        safe_coef(
          coef_null,
          "tau"
        ),

      tau_s_hat_null =
        safe_coef(
          coef_null,
          "tau_s"
        ),

      alpha_hat_null =
        safe_coef(
          coef_null,
          "alpha"
        ),

      alpha_s_hat_null =
        safe_coef(
          coef_null,
          "alpha_s"
        ),

      Base_R =
        latent$params$R,

      Base_tau =
        latent$params$tau,

      Base_tau_s =
        latent$params$tau_s,

      Base_sigma_c =
        latent$params$sigma_c,

      Base_sigma_n =
        latent$params$sigma_n,

      Base_alpha =
        latent$params$alpha,

      Base_alpha_s =
        latent$params$alpha_s,

      runtime_seconds =
        as.numeric(
          difftime(
            Sys.time(),
            t0,
            units = "secs"
          )
        )
    )
  }

  rbindlist(
    out,
    use.names = TRUE,
    fill = TRUE
  )
}


# =============================================================================
# 11. Task/checkpoint helpers
# =============================================================================

tasks <- rbindlist(
  list(
    data.table(
      truth = "null",
      Gene = seq_len(
        N_NULL
      )
    ),
    data.table(
      truth = "alternative",
      Gene = seq_len(
        N_ALT
      )
    )
  )
)

tasks[
  ,
  task_id :=
    .I
]


checkpoint_file <- function(
  truth,
  gene_i
) {

  file.path(
    CHECKPOINT_DIR,
    sprintf(
      "%s_gene_%06d.rds",
      truth,
      gene_i
    )
  )
}


valid_checkpoint <- function(
  x,
  truth,
  gene_i
) {

  if (
    is.null(x) ||
      !is.data.table(x)
  ) {
    return(FALSE)
  }

  required <- c(
    "Benchmark_version",
    "truth",
    "Gene",
    "Exprs_noise",
    "onset_time"
  )

  if (
    any(
      !required %in%
        names(x)
    )
  ) {
    return(FALSE)
  }

  if (
    nrow(x) !=
      EXPECTED_TESTS_PER_GENE
  ) {
    return(FALSE)
  }

  if (
    uniqueN(
      x$Benchmark_version
    ) != 1L ||
      x$Benchmark_version[1] !=
        BENCHMARK_VERSION
  ) {
    return(FALSE)
  }

  if (
    uniqueN(
      x$truth
    ) != 1L ||
      x$truth[1] !=
        truth
  ) {
    return(FALSE)
  }

  if (
    uniqueN(
      x$Gene
    ) != 1L ||
      x$Gene[1] !=
        gene_i
  ) {
    return(FALSE)
  }

  if (
    uniqueN(
      x$Exprs_noise
    ) !=
      EXPECTED_TESTS_PER_GENE
  ) {
    return(FALSE)
  }

  if (
    uniqueN(
      x$onset_time
    ) != 1L
  ) {
    return(FALSE)
  }

  TRUE
}


run_task_checkpoint <- function(
  task_row,
  N_boot
) {

  truth_i <- as.character(
    task_row$truth
  )

  gene_i <- as.integer(
    task_row$Gene
  )

  ff <- checkpoint_file(
    truth_i,
    gene_i
  )

  if (
    file.exists(ff)
  ) {

    old <- tryCatch(
      readRDS(ff),
      error = function(e) NULL
    )

    if (
      valid_checkpoint(
        old,
        truth_i,
        gene_i
      )
    ) {

      return(
        data.table(
          task_id =
            task_row$task_id,
          truth =
            truth_i,
          Gene =
            gene_i,
          status =
            "already_checkpointed",
          error_message =
            NA_character_
        )
      )
    }
  }

  ans <- tryCatch(
    run_gene(
      gene_i = gene_i,
      truth = truth_i,
      N_boot = N_boot
    ),
    error = function(e) {
      structure(
        list(
          error_message =
            conditionMessage(e)
        ),
        class =
          "gene_error"
      )
    }
  )

  if (
    inherits(
      ans,
      "gene_error"
    )
  ) {

    return(
      data.table(
        task_id =
          task_row$task_id,
        truth =
          truth_i,
        Gene =
          gene_i,
        status =
          "gene_error",
        error_message =
          ans$error_message
      )
    )
  }

  if (
    !valid_checkpoint(
      ans,
      truth_i,
      gene_i
    )
  ) {

    return(
      data.table(
        task_id =
          task_row$task_id,
        truth =
          truth_i,
        Gene =
          gene_i,
        status =
          "invalid_result",
        error_message =
          "Invalid gene-level result."
      )
    )
  }

  tmp <- paste0(
    ff,
    ".tmp.",
    Sys.getpid()
  )

  saveRDS(
    ans,
    tmp,
    compress = "gzip"
  )

  if (
    !file.rename(
      tmp,
      ff
    )
  ) {

    unlink(
      tmp
    )

    return(
      data.table(
        task_id =
          task_row$task_id,
        truth =
          truth_i,
        Gene =
          gene_i,
        status =
          "checkpoint_write_failed",
        error_message =
          "Atomic checkpoint rename failed."
      )
    )
  }

  data.table(
    task_id =
      task_row$task_id,
    truth =
      truth_i,
    Gene =
      gene_i,
    status =
      "completed",
    error_message =
      NA_character_
  )
}


# =============================================================================
# 12. Serial sanity tests
# =============================================================================

for (
  truth_i in c(
    "null",
    "alternative"
  )
) {

  cat(
    "\n============================================================\n",
    "SERIAL SANITY TEST: ",
    truth_i,
    "\n",
    "============================================================\n",
    sep = ""
  )

  test <- run_gene(
    gene_i = 1L,
    truth = truth_i,
    N_boot = 19L
  )

  print(
    test[
      ,
      .(
        truth,
        Exprs_noise,
        sigma_true,
        onset_time,
        status,
        p.value,
        Sigma,
        IR
      )
    ]
  )

  if (
    nrow(test) !=
      EXPECTED_TESTS_PER_GENE
  ) {
    stop(
      paste(
        "Serial sanity test failed for",
        truth_i
      )
    )
  }
}


# =============================================================================
# 13. Determine tasks remaining
# =============================================================================

done <- vapply(
  seq_len(
    nrow(tasks)
  ),
  function(ii) {

    ff <- checkpoint_file(
      tasks$truth[ii],
      tasks$Gene[ii]
    )

    if (!file.exists(ff)) {
      return(FALSE)
    }

    x <- tryCatch(
      readRDS(ff),
      error = function(e) NULL
    )

    valid_checkpoint(
      x,
      tasks$truth[ii],
      tasks$Gene[ii]
    )
  },
  logical(1)
)

tasks_remaining <- tasks[
  !done
]

cat(
  "\nTasks complete: ",
  sum(done),
  " / ",
  nrow(tasks),
  "\n",
  sep = ""
)


# =============================================================================
# 14. Parallel run
# =============================================================================

overall_start <- Sys.time()

if (
  nrow(
    tasks_remaining
  ) > 0L
) {

  workers_use <- min(
    N_WORKERS,
    nrow(
      tasks_remaining
    )
  )

  cl <- makePSOCKcluster(
    workers_use,
    outfile = "",
    timeout =
      PSOCK_TIMEOUT_SECONDS,
    setup_timeout =
      120
  )

  tryCatch(
    {

      clusterEvalQ(
        cl,
        {

          Sys.setenv(
            OMP_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1",
            MKL_NUM_THREADS = "1",
            VECLIB_MAXIMUM_THREADS = "1",
            NUMEXPR_NUM_THREADS = "1"
          )

          library(data.table)
          library(nnls)
          library(MASS)
          library(deSolve)

          data.table::setDTthreads(1)

          source("../ode_model/ode.r")
          source("../commons/nested_test2.r")
          source("../commons/platforms.r")

          NULL
        }
      )

      export_names <- c(
        "N_NULL",
        "N_ALT",
        "N_BOOT",
        "BENCHMARK_VERSION",
        "RELATIVE_TIMES",
        "N_REPLICATES",
        "T_STAR_ABSOLUTE",
        "MAX_ONSET_SHIFT",
        "POST_R_FRACTION",
        "RANGE_NOISE",
        "SCALING_A",
        "LAMBDA_TIME",
        "LAMBDA_DIAG",
        "REL_FLOOR",
        "TRUNCATE_NONNEGATIVE_BOOT",
        "MAX_FAILURE_RATE",
        "T_STAR_FIT",
        "CHECKPOINT_DIR",
        "EXPECTED_TESTS_PER_GENE",
        "stable_seed_from_string",
        "make_seed",
        "add_rnaseq_noise",
        "safe_value",
        "safe_text",
        "safe_coef",
        "make_base_parameters",
        "make_onset_time",
        "simulate_latent_gene",
        "run_gene",
        "checkpoint_file",
        "valid_checkpoint",
        "run_task_checkpoint"
      )

      clusterExport(
        cl,
        export_names,
        envir = .GlobalEnv
      )

      # PSOCK smoke test before long run.
      smoke_rows <- rbindlist(
        list(
          tasks_remaining[
            truth == "null"
          ][1],
          tasks_remaining[
            truth == "alternative"
          ][1]
        ),
        fill = TRUE
      )

      smoke_rows <- smoke_rows[
        !is.na(
          Gene
        )
      ]

      if (
        nrow(
          smoke_rows
        ) > 0L
      ) {

        smoke_list <- split(
          smoke_rows,
          seq_len(
            nrow(
              smoke_rows
            )
          )
        )

        smoke <- parLapplyLB(
          cl,
          smoke_list,
          run_task_checkpoint,
          N_boot = N_BOOT
        )

        smoke <- rbindlist(
          smoke,
          fill = TRUE
        )

        print(
          smoke
        )

        if (
          any(
            !smoke$status %in%
              c(
                "completed",
                "already_checkpointed"
              )
          )
        ) {
          stop(
            "PSOCK smoke test failed."
          )
        }
      }

      # Refresh remaining tasks.
      done2 <- vapply(
        seq_len(
          nrow(tasks)
        ),
        function(ii) {

          ff <- checkpoint_file(
            tasks$truth[ii],
            tasks$Gene[ii]
          )

          if (!file.exists(ff)) {
            return(FALSE)
          }

          x <- tryCatch(
            readRDS(ff),
            error =
              function(e) NULL
          )

          valid_checkpoint(
            x,
            tasks$truth[ii],
            tasks$Gene[ii]
          )
        },
        logical(1)
      )

      rem2 <- tasks[
        !done2
      ]

      if (
        nrow(rem2) > 0L
      ) {

        batch_ids <- split(
          seq_len(
            nrow(rem2)
          ),
          ceiling(
            seq_len(
              nrow(rem2)
            ) /
              TASK_BATCH_SIZE
          )
        )

        for (
          bb in seq_along(
            batch_ids
          )
        ) {

          rows <- rem2[
            batch_ids[[bb]]
          ]

          task_list <- split(
            rows,
            seq_len(
              nrow(rows)
            )
          )

          cat(
            "\nBatch ",
            bb,
            " / ",
            length(batch_ids),
            " | tasks = ",
            length(task_list),
            "\n",
            sep = ""
          )

          tbb <- Sys.time()

          stat <- parLapplyLB(
            cl,
            task_list,
            run_task_checkpoint,
            N_boot = N_BOOT
          )

          stat <- rbindlist(
            stat,
            fill = TRUE
          )

          stat[
            ,
            `:=`(
              task_batch = bb,
              batch_elapsed_min =
                as.numeric(
                  difftime(
                    Sys.time(),
                    tbb,
                    units = "mins"
                  )
                ),
              timestamp =
                as.character(
                  Sys.time()
                )
            )
          ]

          fwrite(
            stat,
            PROGRESS_FILE,
            sep = "\t",
            append =
              file.exists(
                PROGRESS_FILE
              ),
            col.names =
              !file.exists(
                PROGRESS_FILE
              )
          )

          print(
            stat[
              ,
              .N,
              by = .(
                truth,
                status
              )
            ]
          )

          errs <- stat[
            !is.na(
              error_message
            )
          ]

          if (
            nrow(errs) > 0L
          ) {
            print(
              unique(
                errs[
                  ,
                  .(
                    status,
                    error_message
                  )
                ]
              )
            )
          }
        }
      }

    },

    finally = {

      try(
        stopCluster(
          cl
        ),
        silent = TRUE
      )
    }
  )
}


# =============================================================================
# 15. Collect checkpoints
# =============================================================================

res <- list()

for (
  ii in seq_len(
    nrow(tasks)
  )
) {

  ff <- checkpoint_file(
    tasks$truth[ii],
    tasks$Gene[ii]
  )

  if (!file.exists(ff)) {
    next
  }

  x <- tryCatch(
    readRDS(ff),
    error = function(e) NULL
  )

  if (
    valid_checkpoint(
      x,
      tasks$truth[ii],
      tasks$Gene[ii]
    )
  ) {

    res[[
      length(res) +
        1L
    ]] <- x
  }
}

if (
  length(res) == 0L
) {
  stop(
    "No valid checkpoints recovered."
  )
}

results <- rbindlist(
  res,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  results,
  truth,
  Gene,
  Exprs_noise
)

results[
  ,
  valid :=
    is.finite(
      p.value
    ) &
    p.value >= 0 &
    p.value <= 1 &
    !status %in%
      c(
        "test_error",
        "measurement_error"
      )
]


# =============================================================================
# 16. Wilson interval
# =============================================================================

wilson_interval <- function(
  successes,
  n,
  conf = 0.95
) {

  if (
    n <= 0L
  ) {
    return(
      c(
        low = NA_real_,
        high = NA_real_
      )
    )
  }

  z <- qnorm(
    1 -
      (1 - conf) /
        2
  )

  p <- successes / n

  den <- 1 + z^2 / n

  ctr <- (
    p +
      z^2 /
        (2 * n)
  ) /
    den

  half <- (
    z *
      sqrt(
        p *
          (1 - p) /
          n +
          z^2 /
            (4 * n^2)
      )
  ) /
    den

  c(
    low =
      max(
        0,
        ctr - half
      ),
    high =
      min(
        1,
        ctr + half
      )
  )
}


# =============================================================================
# 17. Matched summary
# =============================================================================

summary_by_noise <- results[
  ,
  {

    null <- .SD[
      valid &
        truth == "null"
    ]

    alt <- .SD[
      valid &
        truth == "alternative"
    ]

    n0 <- nrow(null)
    n1 <- nrow(alt)

    typeI <- if (
      n0 > 0L
    ) {
      mean(
        null$p.value <=
          0.05
      )
    } else {
      NA_real_
    }

    ci <- if (
      n0 > 0L
    ) {
      wilson_interval(
        sum(
          null$p.value <=
            0.05
        ),
        n0
      )
    } else {
      c(
        low = NA_real_,
        high = NA_real_
      )
    }

    power <- if (
      n1 > 0L
    ) {
      mean(
        alt$p.value <=
          0.05
      )
    } else {
      NA_real_
    }

    sigma_rmse <- if (
      n1 > 0L
    ) {
      sqrt(
        mean(
          (
            alt$Sigma -
              alt$sigma_true
          )^2,
          na.rm = TRUE
        )
      )
    } else {
      NA_real_
    }

    .(
      N_null_valid =
        n0,

      N_alt_valid =
        n1,

      TypeI_005 =
        typeI,

      TypeI_005_Wilson_low =
        unname(
          ci["low"]
        ),

      TypeI_005_Wilson_high =
        unname(
          ci["high"]
        ),

      TypeI_CI_contains_005 =
        is.finite(
          ci["low"]
        ) &&
        ci["low"] <= 0.05 &&
        ci["high"] >= 0.05,

      Power_005 =
        power,

      sigma_c_RMSE_alt =
        sigma_rmse,

      Null_fraction_p1 =
        mean(
          null$p.value == 1,
          na.rm = TRUE
        ),

      Null_fraction_sigma0 =
        mean(
          null$Sigma <= 1e-12,
          na.rm = TRUE
        ),

      Median_condition =
        median(
          condition_number,
          na.rm = TRUE
        ),

      Finite_condition_fraction =
        mean(
          is.finite(
            condition_number
          )
        ),

      Mean_bootstrap_failure_rate =
        mean(
          bootstrap_failure_rate,
          na.rm = TRUE
        ),

      Mean_boot_rank_deficient =
        mean(
          boot_rank_deficient_fraction,
          na.rm = TRUE
        )
    )
  },
  by = Exprs_noise
]

summary_by_noise[
  ,
  noise_order :=
    match(
      Exprs_noise,
      RANGE_NOISE
    )
]

setorder(
  summary_by_noise,
  noise_order
)

summary_by_noise[
  ,
  noise_order :=
    NULL
]


# =============================================================================
# 18. Save
# =============================================================================

fwrite(
  results,
  RAW_FILE,
  sep = "\t"
)

fwrite(
  summary_by_noise,
  SUMMARY_FILE,
  sep = "\t"
)

matched_benchmark <- list(
  version =
    BENCHMARK_VERSION,

  design =
    list(
      relative_times =
        RELATIVE_TIMES,
      N_replicates =
        N_REPLICATES,
      T_star =
        T_STAR_FIT,
      max_onset_shift =
        MAX_ONSET_SHIFT,
      post_R_fraction =
        POST_R_FRACTION,
      platform =
        "RNA-seq"
    ),

  results =
    results,

  summary =
    summary_by_noise
)

save(
  matched_benchmark,
  file =
    RDATA_FILE
)

sink(
  SESSION_FILE
)

print(
  sessionInfo()
)

sink()


# =============================================================================
# 19. Final report
# =============================================================================

cat(
  "\n============================================================\n",
  "GSE256335 MATCHED BENCHMARK COMPLETE\n",
  "============================================================\n",
  sep = ""
)

print(
  summary_by_noise
)

cat(
  "\nSaved:\n",
  "  ",
  RAW_FILE,
  "\n",
  "  ",
  SUMMARY_FILE,
  "\n",
  "  ",
  RDATA_FILE,
  "\n",
  "  ",
  PROGRESS_FILE,
  "\n",
  "  ",
  SESSION_FILE,
  "\n",
  sep = ""
)
