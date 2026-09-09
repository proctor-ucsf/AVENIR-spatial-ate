#--------------------------------------------
#
# AVENIR-spatial-ate/analysis/07x-functions.R
#
# Helper functions for the targeting-rule addendum: does prioritising grappes
# by predicted mortality recover the benefit of prioritising by predicted
# treatment effect?
#
# Sourced by 07-targeting-analysis.qmd only, and sourced AFTER 0-config.R and
# 0-functions.R: the geospatial parameters (global_crs, proj_crs), the package
# attachments and the shared raster helpers all come from there rather than
# being redeclared here. 07 is deliberately not part of run-all.R.
#
# Design note: every function below except build_prediction_table() takes a
# grappe-level *prediction table*, never a model object, and carries a
# `replicate` column through untouched. That is what lets the bootstrap be
# layered on by supplying B prediction tables and mapping these functions over
# them, with no function body needing to change. The same property is what
# makes the three-model version of the document cheap: the model is just
# another dimension the caller maps over.
#
#--------------------------------------------


#--------------------------------------
# build_prediction_table()
#
# Assembles the grappe-level prediction table that every other function in
# this file consumes. For each grappe it predicts the child mortality rate
# under both arms from a fitted spaMM geostatistical model, differences them
# to get a predicted risk difference, attaches the doses actually delivered,
# and converts to predicted deaths averted.
#
# Arguments:
# - spamm_model_fit: spaMM HLfit object; the fitted arm x space Poisson GLMM.
#                    Its formula must use `lon` and `lat` as the coordinates
#                    and `isTreated_child` as the arm indicator.
# - grappe_dat:      data.frame, one row per grappe. Must contain site_code,
#                    grappe, csi, district, lat, lon (decimal degrees), and
#                    isTreated_child (factor, for its level ordering).
# - phase_dat:       data.frame, one row per grappe per phase. Must contain
#                    site_code, phase, started_phase_1to59. Used for the
#                    eligible-child denominator only.
# - dose_dat:        data.frame, one row per grappe per distribution round.
#                    Must contain site_code, phase, n_doses. May cover more
#                    grappes than grappe_dat (the AVENIR dose file is
#                    trial-wide); it is subset to grappe_dat internally.
# - phase_years:     numeric; protection conferred by one round, in years
#                    (default 0.5, i.e. six-monthly rounds).
# - global_crs:      integer EPSG code of the coordinates in grappe_dat
#                    (default 4326, WGS84 lat/lon).
# - proj_crs:        integer EPSG code the model was fitted in
#                    (default 32631, UTM zone 31N). Predictions must be made
#                    in the same projection the Matern field was fitted in.
# - replicate:       integer; replicate identifier carried through every
#                    downstream function (default 1 = the observed data).
#
# Details:
# - person_time is the grappe's observed child-years at risk (pt_year_1to59
#   from grappe_dat, one value per grappe). It is the single exposure unit the
#   downstream targeting curves run on: benefit is RD * person_time / 1000, the
#   curve accumulates person_time on the x-axis, and every weighted mean and
#   concordance downstream weights by it. Because cost and benefit then share
#   person_time, ordering by RD is the exact optimal (benefit-per-person-time)
#   greedy order and the `rd` curve is a true envelope.
# - N_i is n_doses: study-drug doses actually delivered to children 1-59
#   months, summed over a grappe's rounds. Doses are counted in every arm, not
#   just child azithromycin - a placebo grappe's dose count measures how many
#   children a distribution round in that grappe reaches. It is retained as a
#   delivery / coverage descriptor only; since the person-time reframing it is
#   no longer a weight, a ranking signal, or the curve's budget.
# - n_children is the sum of started_phase_1to59 across a grappe's phases.
#   started_phase_1to59 counts a child once per phase, so this is a child-
#   *round* total on the same footing as n_doses, not a count of distinct
#   children: n_doses / n_children is coverage per round, not the share of
#   children ever reached. It is reported for context and is not used as a
#   weight anywhere.
# - A grappe-round with no row in dose_dat contributes zero doses. That is the
#   intended reading (no record = no distribution) and it correctly pushes
#   hard-to-reach grappes down the targeting order, but it means a data-entry
#   gap would be indistinguishable from a genuine non-distribution.
# - Predictions use a fixed offset of pt_year_1to59 = 1000, so mu0 and mu1 are
#   on a deaths-per-1000-child-years scale (the same convention as
#   get_grid_preds() in 0-functions.R). Note this prediction offset is a
#   constant 1000 and is unrelated to the `person_time` column returned below,
#   which carries each grappe's *observed* child-years and drives the curves.
# - The prediction newdata carries district and site_code even though the
#   interaction model does not use them, and a numeric isTreated_child_num
#   alongside the factor arm indicator. Between them these cover the
#   interaction, district-fixed-effect and two-Matern specifications from 02A.
#   mordorstate is not constructed, so the MORDOR model is not supported here.
# - Predicting at the grappes' own locations means the Matern BLUP is
#   informed by those grappes' own outcomes. That is deliberate (see the Setup
#   notes in 07-targeting-analysis.qmd) but it is the main reason these results
#   are an upper bound.
#
# Requirements:
# - Packages: dplyr, sf, spaMM, tibble
#
# Returns:
# - A tibble with one row per grappe and columns: replicate, site_code,
#   grappe, csi, district, lat, lon, n_phases, n_children, n_rounds, n_doses,
#   person_time, mu0, mu1, rd, deaths_averted.
#--------------------------------------

build_prediction_table <- function(spamm_model_fit,
                                   grappe_dat,
                                   phase_dat,
                                   dose_dat,
                                   phase_years = 0.5,
                                   global_crs = 4326,
                                   proj_crs = 32631,
                                   replicate = 1) {

  required_grappe <- c("site_code", "grappe", "csi", "district",
                       "lat", "lon", "isTreated_child", "pt_year_1to59")
  required_phase <- c("site_code", "phase", "started_phase_1to59")
  required_dose <- c("site_code", "phase", "n_doses")

  if (!all(required_grappe %in% names(grappe_dat))) {
    stop("grappe_dat is missing: ",
         paste(setdiff(required_grappe, names(grappe_dat)), collapse = ", "))
  }
  if (!all(required_phase %in% names(phase_dat))) {
    stop("phase_dat is missing: ",
         paste(setdiff(required_phase, names(phase_dat)), collapse = ", "))
  }
  if (!all(required_dose %in% names(dose_dat))) {
    stop("dose_dat is missing: ",
         paste(setdiff(required_dose, names(dose_dat)), collapse = ", "))
  }

  replicate_id <- replicate

  # Eligible child population per grappe: summed over phases -----
  pop <- phase_dat %>%
    group_by(site_code) %>%
    summarise(
      n_phases   = dplyr::n_distinct(phase),
      n_children = sum(started_phase_1to59, na.rm = TRUE),
      .groups    = "drop"
    )

  # Doses delivered per grappe: summed over rounds -----
  # dose_dat is trial-wide, so it is subset to the analytic grappes first.
  doses <- dose_dat %>%
    dplyr::filter(site_code %in% grappe_dat$site_code) %>%
    group_by(site_code) %>%
    summarise(
      n_rounds = dplyr::n_distinct(phase),
      n_doses  = sum(n_doses, na.rm = TRUE),
      .groups  = "drop"
    )

  # person_time is the grappe's observed child-years at risk. Renamed from
  # pt_year_1to59 on the way in so it cannot be confused with the constant
  # pt_year_1to59 = 1000 prediction offset set on `newd` below.
  base_dat <- grappe_dat %>%
    ungroup() %>%
    dplyr::select(site_code, grappe, csi, district, lat, lon,
                  person_time = pt_year_1to59) %>%
    # one-to-one: a row-multiplying join would sail past the is.na() checks
    # below and silently double-count children or doses
    dplyr::left_join(pop, by = "site_code", relationship = "one-to-one") %>%
    dplyr::left_join(doses, by = "site_code", relationship = "one-to-one")

  if (any(is.na(base_dat$n_children))) {
    stop("No phase-level population for ",
         sum(is.na(base_dat$n_children)), " grappe(s); check the site_code join.")
  }
  if (any(is.na(base_dat$n_doses))) {
    stop("No dose records at all for ",
         sum(is.na(base_dat$n_doses)), " grappe(s); check the site_code join. ",
         "(A grappe missing only *some* rounds is fine - those count as zero ",
         "doses - but a grappe absent from dose_dat entirely is a join error.)")
  }

  # Project to the CRS the Matern field was fitted in -----
  # The model stores lon/lat in metres, so predicting in decimal degrees would
  # silently place every grappe outside the fitted field.
  coords_proj <- base_dat %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = global_crs) %>%
    sf::st_transform(crs = proj_crs)

  arm_levels <- levels(grappe_dat$isTreated_child)
  if (!all(c("Control", "Treatment") %in% arm_levels)) {
    stop("isTreated_child must have levels 'Control' and 'Treatment'.")
  }

  newd <- base_dat %>%
    mutate(
      lon           = sf::st_coordinates(coords_proj)[, 1],
      lat           = sf::st_coordinates(coords_proj)[, 2],
      district      = as.factor(district),
      site_code     = as.character(site_code),
      pt_year_1to59 = 1000 # offset: predictions are per 1000 child-years
    )

  predict_arm <- function(arm) {
    newd_arm <- newd %>%
      mutate(
        isTreated_child = factor(arm, levels = arm_levels),
        # Numeric duplicate of the arm indicator. The interaction model does not
        # reference it, but the two-Matern specification carries a random slope
        # on it - Matern(0 + isTreated_child_num | lon + lat) - and predict()
        # errors on a missing column rather than ignoring it. Same construction
        # as get_grid_preds() in 0-functions.R.
        isTreated_child_num = as.integer(arm == "Treatment")
      )
    as.numeric(predict(object  = spamm_model_fit,
                       type    = "response",
                       newdata = newd_arm))
  }

  # Assemble -----
  # Built from base_dat, not newd, so the returned lat/lon stay in decimal
  # degrees for mapping; only the prediction newdata uses projected metres.
  # rd is signed so that positive = benefit, i.e. the negative of the ATE
  # (predTreatment - predControl) reported elsewhere in this repo.
  #
  # deaths_averted uses observed person_time as the exposure denominator, not
  # n_doses * phase_years. Because cost (person_time on the curve's x-axis) and
  # benefit then share person_time, deaths averted per child-year is exactly
  # rd / 1000, so ordering by rd is the optimal greedy order and the `rd` curve
  # is a true envelope. phase_years is retained as an argument for the coverage
  # descriptives but no longer enters the benefit.
  out <- base_dat %>%
    mutate(
      replicate      = replicate_id,
      mu0            = predict_arm("Control"),
      mu1            = predict_arm("Treatment"),
      rd             = mu0 - mu1,
      deaths_averted = rd * person_time / 1000
    ) %>%
    dplyr::select(replicate, site_code, grappe, csi, district, lat, lon,
                  n_phases, n_children, n_rounds, n_doses, person_time,
                  mu0, mu1, rd, deaths_averted)

  tibble::as_tibble(out)
}


#--------------------------------------
# aggregate_predictions()
#
# Dose-weighted aggregation of a grappe-level prediction table up to the
# CSI or department scale. Coarser scales are aggregations of the same grappe
# predictions rather than separate model refits, so all three scales rest on
# one prediction table and are mutually consistent.
#
# Arguments:
# - dat:   tibble from build_prediction_table(). Must contain replicate,
#          site_code, csi, district, n_children, n_doses, person_time, mu0,
#          mu1, rd, deaths_averted.
# - level: character; one of "grappe", "csi", "department".
#
# Details:
# - mu0, mu1 and rd are averaged with weights person_time; n_children, n_doses,
#   person_time and deaths_averted are summed. Because rd is person_time-
#   weighted, the summed deaths_averted equals rd * person_time / 1000
#   recomputed at the coarse scale exactly, so the three scales cannot drift
#   apart. Person-time weighting (rather than doses) makes the aggregated rate
#   the unit's actual person-time-weighted rate, consistent with the crude
#   arm-specific rates it is compared against.
# - level = "grappe" is a pass-through that returns the same column set as the
#   other two levels, so callers can map over all three uniformly. site_code
#   and grappe survive only as unit_id / unit_label, and lat, lon, n_phases and
#   n_rounds are dropped. rank_clusters() and compute_targeting_curve() both
#   reference site_code by name, so they must be given the raw
#   build_prediction_table() output, not this.
# - Aggregation is within replicate, so a stacked bootstrap table aggregates
#   correctly without any regrouping by the caller.
#
# Requirements:
# - Packages: dplyr, tibble
#
# Returns:
# - A tibble with one row per unit per replicate: replicate, level, unit_id,
#   unit_label, district (NA-free at every level), csi (NA at department
#   level), n_units, n_children, n_doses, person_time, mu0, mu1, rd,
#   deaths_averted.
#--------------------------------------

aggregate_predictions <- function(dat, level = c("grappe", "csi", "department")) {

  level <- match.arg(level)

  if (level == "grappe") {
    return(
      dat %>%
        mutate(level      = "grappe",
               unit_id    = site_code,
               unit_label = grappe,
               n_units    = 1L) %>%
        dplyr::select(replicate, level, unit_id, unit_label, district, csi,
                      n_units, n_children, n_doses, person_time, mu0, mu1, rd,
                      deaths_averted) %>%
        tibble::as_tibble()
    )
  }

  group_col <- if (level == "csi") "csi" else "district"

  # `level_name` is a plain local, deliberately not called `level`: once
  # mutate() creates a `level` column, any later reference to `level` inside
  # the same mutate() resolves to that length-n column, and `if (level ==
  # "csi")` would then be a condition of length n.
  level_name <- level

  # NOTE: summarise() evaluates its arguments in order and later expressions
  # see earlier *results*, so the three weighted means must be computed before
  # person_time is collapsed to a sum - otherwise weighted.mean() would be
  # handed a single number as its weight vector.
  out <- dat %>%
    mutate(unit_id = as.character(.data[[group_col]])) %>%
    group_by(replicate, unit_id) %>%
    summarise(
      district       = dplyr::first(district),
      n_units        = dplyr::n(),
      mu0            = stats::weighted.mean(mu0, w = person_time, na.rm = TRUE),
      mu1            = stats::weighted.mean(mu1, w = person_time, na.rm = TRUE),
      rd             = stats::weighted.mean(rd, w = person_time, na.rm = TRUE),
      deaths_averted = sum(deaths_averted, na.rm = TRUE),
      n_children     = sum(n_children, na.rm = TRUE),
      n_doses        = sum(n_doses, na.rm = TRUE),
      person_time    = sum(person_time, na.rm = TRUE),
      .groups        = "drop"
    ) %>%
    mutate(
      level      = level_name,
      unit_label = unit_id,
      csi        = if (level_name == "csi") unit_id else NA_character_
    ) %>%
    dplyr::select(replicate, level, unit_id, unit_label, district, csi,
                  n_units, n_children, n_doses, person_time, mu0, mu1, rd,
                  deaths_averted)

  tibble::as_tibble(out)
}


#--------------------------------------
# rank_clusters()
#
# Orders the grappes in a prediction table according to one targeting rule and
# records where the interpretable knots of the resulting curve fall.
#
# Five rules:
#   rd           descending predicted risk difference. The upper envelope:
#                because deaths averted per child-year is rd / 1000, ordering
#                by rd is ordering by benefit per unit person-time, which is
#                the optimal greedy order for a budget expressed in person-
#                time. Not achievable prospectively, since rd is not knowable
#                in advance.
#   mort_grappe  descending predicted control-arm mortality. The operational
#                rule of interest.
#   mort_csi     CSIs ordered by person-time-weighted mu0; all grappes in a CSI
#                treated as one block.
#   mort_dept    departments as blocks. Coarsest targeting.
#   random       population-proportional allocation; the no-information
#                reference, which is the 45-degree line by construction.
#
# Arguments:
# - dat:  grappe-level tibble from build_prediction_table(). Must contain
#         replicate, site_code, csi, district, n_doses, person_time, mu0, rd.
# - rule: character; one of "rd", "mort_grappe", "mort_csi", "mort_dept",
#         "random".
#
# Details:
# - Block rules order grappes arbitrarily *within* a block, so only the block
#   boundaries are interpretable. block_boundary marks the last grappe of each
#   block in rank order; compute_targeting_curve() keeps only those knots and
#   interpolates linearly between them. Within-block order is set by site_code
#   so the arbitrary part is at least reproducible.
# - For "random", the whole region is one block, which reduces the curve to
#   the straight line from (0, 0) to (1, 1) - exactly population-proportional
#   allocation, with no simulation and therefore no seed dependence.
# - Grappe-level rules break ties on site_code for reproducibility.
# - Ranking is within replicate.
#
# Requirements:
# - Packages: dplyr
#
# Returns:
# - dat with four columns added: rule, rank_order (integer, 1 = treated
#   first), block_id (character), block_boundary (logical). A block_score
#   column is created internally and dropped before returning.
#--------------------------------------

rank_clusters <- function(dat,
                          rule = c("rd", "mort_grappe", "mort_csi",
                                   "mort_dept", "random")) {

  rule <- match.arg(rule)

  # Score each grappe by whatever the rule orders on. For block rules the
  # score is constant within a block, so every grappe in a block sorts
  # together and the block moves as a unit.
  scored <- switch(
    rule,
    rd = dat %>%
      mutate(block_id = as.character(site_code), block_score = rd),
    mort_grappe = dat %>%
      mutate(block_id = as.character(site_code), block_score = mu0),
    mort_csi = dat %>%
      group_by(replicate, csi) %>%
      mutate(block_id    = as.character(csi),
             block_score = stats::weighted.mean(mu0, w = person_time,
                                                na.rm = TRUE)) %>%
      ungroup(),
    mort_dept = dat %>%
      group_by(replicate, district) %>%
      mutate(block_id    = as.character(district),
             block_score = stats::weighted.mean(mu0, w = person_time,
                                                na.rm = TRUE)) %>%
      ungroup(),
    random = dat %>%
      mutate(block_id = "all", block_score = 0)
  )

  scored %>%
    arrange(replicate, dplyr::desc(block_score), block_id, site_code) %>%
    group_by(replicate) %>%
    mutate(
      # .env$ so this keeps working if a caller ever passes a table that
      # already has a `rule` column, as a stacked bootstrap table might
      rule           = .env$rule,
      rank_order     = dplyr::row_number(),
      # last grappe of each block in rank order; TRUE for every grappe under
      # the two grappe-level rules, once per block otherwise
      block_boundary = block_id != dplyr::lead(block_id,
                                               default = "__end_of_ranking__")
    ) %>%
    ungroup() %>%
    dplyr::select(-block_score)
}


#--------------------------------------
# compute_targeting_curve()
#
# Cumulative predicted deaths averted against cumulative person-time covered
# for one ranking rule, evaluated on a caller-supplied shared x-grid.
#
# The y-axis is a percentage of the universal-treatment total, so any common
# scaling of the benefit cancels and the curve is a pure share. Note that
# the curve is not necessarily monotone: wherever the model predicts rd < 0
# (harm), treating those grappes subtracts from the cumulative total, so the
# rd curve typically rises above 100% and comes back down to 100% at full
# coverage. That is a feature of the predicted surface, not a coding error.
#
# Arguments:
# - dat:    grappe-level tibble from build_prediction_table().
# - rule:   character; passed to rank_clusters().
# - x_grid: numeric vector of cumulative coverage values in [0, 1] expressed
#           as a share of all person-time, the shared grid all rules are
#           evaluated on (default seq(0, 1, by = 0.005)).
#
# Details:
# - Only block boundaries are used as knots, with (0, 0) prepended, and
#   approx() interpolates linearly between them. For the grappe-level rules
#   every grappe is a boundary, so this is the exact step-free curve; for
#   block rules it is the linear interpolation across whole blocks that the
#   published convention for comparing targeting scales uses.
# - Evaluating on a fixed shared grid rather than at observed knots is what
#   makes the deferred bootstrap summarisable pointwise into quantile ribbons
#   with no interpolation mismatch.
# - approx(rule = 2) clamps rather than returning NA at the ends, so x = 0 and
#   x = 1 are well defined.
#
# Requirements:
# - Depends on: rank_clusters()
# - Packages: dplyr, purrr, tibble
#
# Returns:
# - A tibble with one row per replicate per grid point: replicate, rule, x
#   (proportion of doses delivered), doses_delivered (cumulative doses),
#   prop_averted (proportion of universal-treatment deaths averted),
#   deaths_averted (count).
#--------------------------------------

compute_targeting_curve <- function(dat, rule, x_grid = seq(0, 1, by = 0.005)) {

  if (any(x_grid < 0 | x_grid > 1)) {
    stop("x_grid must lie in [0, 1]; it is a proportion of person-time.")
  }

  ranked <- rank_clusters(dat, rule = rule)

  curve_one <- function(rep_dat) {

    total_pt      <- sum(rep_dat$person_time, na.rm = TRUE)
    total_averted <- sum(rep_dat$deaths_averted, na.rm = TRUE)

    if (total_pt <= 0 || total_averted == 0) {
      warning("Degenerate totals for rule '", rule, "'; returning NA curve.")
      return(tibble::tibble(x = x_grid, prop_averted = NA_real_,
                            person_time = NA_real_,
                            deaths_averted = NA_real_))
    }
    if (total_averted < 0) {
      warning("Total predicted deaths averted is negative; the percentage ",
              "y-axis is inverted and the curve should not be read as usual.")
    }

    knots <- rep_dat %>%
      arrange(rank_order) %>%
      mutate(cum_x = cumsum(person_time) / total_pt,
             cum_y = cumsum(deaths_averted) / total_averted) %>%
      dplyr::filter(block_boundary)

    tibble::tibble(
      x            = x_grid,
      # ties = "ordered": cum_x is nondecreasing by construction, and this
      # stops approx() from warning (and averaging) if a zero-population unit
      # ever produces a duplicated knot.
      prop_averted = stats::approx(x    = c(0, knots$cum_x),
                                   y    = c(0, knots$cum_y),
                                   xout = x_grid,
                                   rule = 2,
                                   ties = "ordered")$y
    ) %>%
      mutate(person_time    = x * total_pt,
             deaths_averted  = prop_averted * total_averted)
  }

  ranked %>%
    dplyr::group_split(replicate) %>%
    purrr::map_dfr(function(rep_dat) {
      curve_one(rep_dat) %>%
        mutate(replicate = rep_dat$replicate[1], rule = rule)
    }) %>%
    dplyr::select(replicate, rule, x, person_time, prop_averted,
                  deaths_averted)
}


#--------------------------------------
# rule_concordance()
#
# The AUROC of a targeting rule against the predicted-effect ranking: for a
# randomly drawn pair of villages, weighted by person-time, the probability
# that the rule ranks the village with the larger predicted RD ahead of the
# other one.
#
# This is the per-rule counterpart of the efficiency ratio, and it lives on the
# same scale with the same two anchors:
#
#   rule = "rd"      -> C = 1   exactly (the rule is the comparator)
#   rule = "random"  -> C = 0.5 exactly (one block, every pair tied)
#
# Arguments:
# - dat:  grappe-level tibble from build_prediction_table().
# - rule: character; passed to rank_clusters().
#
# Details:
# - The score a rule is judged on is its *block* position, not the grappe's
#   position within the block, so every pair of villages inside one block is a
#   tie. That is the point: a block rule genuinely cannot discriminate within a
#   block, and scoring those pairs at 0.5 (the AUROC convention for ties in the
#   predictor) is what charges the rule for it. Taking rank_order directly would
#   silently credit the arbitrary within-block ordering that rank_clusters()
#   imposes on site_code, inflating C for the two block rules.
# - Block scores come from rank_clusters(), so this statistic and the targeting
#   curves are guaranteed to be ranking on the same thing.
# - Not the same as a concordance computed on aggregated CSI or department
#   units. That would answer "do CSIs rank consistently against each other",
#   which is a different question from "what does the CSI block rule cost",
#   because it never sees the within-block pairs at all.
# - Implemented with outer products: O(n^2) memory. At the AVENIR grappe count
#   (n ~ 1400) that is ~60 MB; the function stops above max_n rather than
#   silently exhausting memory.
#
# Requirements:
# - Depends on: rank_clusters()
# - Packages: dplyr
#
# Returns:
# - A single numeric in [0, 1], or NA_real_ if the statistic is undefined.
#--------------------------------------

rule_concordance <- function(dat, rule, max_n = 5000) {

  ranked <- rank_clusters(dat, rule = rule) %>%
    group_by(replicate, block_id) %>%
    mutate(block_rank = min(rank_order)) %>%
    ungroup()

  # higher score = treated earlier, so negate the rank
  score <- -ranked$block_rank
  rd    <- ranked$rd
  w     <- ranked$person_time

  keep  <- is.finite(score) & is.finite(rd) & is.finite(w) & w > 0
  score <- score[keep]; rd <- rd[keep]; w <- w[keep]

  n <- length(score)
  if (n < 2) {
    warning("Fewer than 2 usable units; returning NA.")
    return(NA_real_)
  }
  if (n > max_n) {
    stop("n = ", n, " exceeds max_n = ", max_n,
         "; the O(n^2) pairwise sum would be too large.")
  }

  weight_mat <- outer(w, w)
  diag(weight_mat) <- 0

  ahead <- outer(score, score, ">")            # rule puts i before j
  tied  <- outer(score, score, "==")           # same block: no discrimination
  diag(tied) <- FALSE

  rd_agree <- outer(rd, rd, ">") + 0.5 * outer(rd, rd, "==")
  diag(rd_agree) <- 0

  # tied pairs are counted twice by the symmetric outer product, hence the
  # halving; each contributes 0.5 to the numerator
  denominator <- sum(weight_mat * ahead) + 0.5 * sum(weight_mat * tied)

  if (!is.finite(denominator) || denominator <= 0) {
    warning("No usable pairs; returning NA.")
    return(NA_real_)
  }

  (sum(weight_mat * ahead * rd_agree) + 0.25 * sum(weight_mat * tied)) / denominator
}


#--------------------------------------
# weighted_correlation()
#
# Weighted Spearman and Pearson correlations between two vectors. The weights
# are caller-supplied; every call in this project passes person-time.
#
# Arguments:
# - x, y: numeric vectors of equal length.
# - w:    numeric vector of the same length; weights, here person-time.
#
# Details:
# - The weighted Spearman coefficient is the weighted Pearson correlation of
#   the mid-ranks (rank(ties.method = "average")); ranks themselves are
#   unweighted, which is the standard construction.
# - Covariances use the weighted (biased) form, sum w (a - abar)(b - bbar) /
#   sum w. The normalisation cancels in the correlation, so this choice does
#   not affect the reported values.
# - Returns NA rather than erroring - on fewer than three usable observations,
#   on zero weighted variance in either vector, or on any non-finite result -
#   so it is safe to call inside a future resampling loop where an occasional
#   degenerate replicate should not abort the run.
#
# Requirements:
# - Packages: tibble
#
# Returns:
# - A one-row tibble: n (units used), rho (weighted Spearman), r (weighted
#   Pearson).
#--------------------------------------

weighted_correlation <- function(x, y, w) {

  if (length(x) != length(y) || length(x) != length(w)) {
    stop("x, y and w must be the same length.")
  }

  na_out <- tibble::tibble(n = NA_integer_, rho = NA_real_, r = NA_real_)

  keep <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  x <- x[keep]
  y <- y[keep]
  w <- w[keep]

  n <- length(x)

  if (n < 3) {
    warning("Fewer than 3 usable observations; returning NA.")
    return(na_out)
  }

  weighted_pearson <- function(a, b) {
    a_bar <- sum(w * a) / sum(w)
    b_bar <- sum(w * b) / sum(w)
    cov_ab <- sum(w * (a - a_bar) * (b - b_bar)) / sum(w)
    var_a  <- sum(w * (a - a_bar)^2) / sum(w)
    var_b  <- sum(w * (b - b_bar)^2) / sum(w)
    if (var_a <= 0 || var_b <= 0) return(NA_real_)
    cov_ab / sqrt(var_a * var_b)
  }

  rho <- weighted_pearson(rank(x, ties.method = "average"),
                          rank(y, ties.method = "average"))
  r   <- weighted_pearson(x, y)

  if (!is.finite(rho) || !is.finite(r)) {
    warning("Degenerate weighted variance; returning NA.")
    return(na_out)
  }

  tibble::tibble(n = n, rho = rho, r = r)
}


# ==========================================================================
# Uncertainty: predictive-draw bootstrap
# ==========================================================================

#--------------------------------------
# predictive_draw_tables()
#
# B replicate prediction tables drawn from the fitted model's predictive
# distribution, without refitting.
#
# WHY NOT A REFIT BOOTSTRAP: the two-Matern model this document runs on took
# 2.77 hours to fit, so B refits is not a feasible design. Instead the linear
# predictor is drawn from its joint predictive distribution with the covariance
# parameters (rho, nu, lambda) held at their point estimates. This propagates
# fixed-effect and Matern-field uncertainty but NOT uncertainty in the
# covariance parameters themselves, so the ribbons are mildly optimistic. It is
# the same machinery 02B-posterior-prob-sim.Rmd uses for posterior probabilities.
#
# WHY THE DRAW MUST BE JOINT ACROSS ARMS: rd_i = mu0_i - mu1_i, and the shared
# Matern field w(s_i) is common to both arms, so it cancels in that contrast.
# Drawing the arms independently would destroy the cancellation and inflate the
# variance of rd enormously - the measured cross-arm correlation for a single
# village under the two-Matern model is about 0.79. The covariance is therefore
# requested over all villages x both arms at once and the draw is multivariate
# normal on the linear predictor scale.
#
# Poisson count dispersion is deliberately NOT included: the targeting curves
# are functionals of the predicted rate surface, not of realised counts.
#
# Arguments:
# - spamm_model_fit: spaMM HLfit object.
# - pred_table:      the point-estimate table from build_prediction_table();
#                    supplies unit identifiers, n_doses, and the village order.
# - grappe_dat:      as passed to build_prediction_table(); supplies the
#                    isTreated_child factor levels and the unprojected coords.
# - B:               number of replicates.
# - phase_years, global_crs, proj_crs: as in build_prediction_table().
# - seed:            RNG seed for the draws.
#
# Details:
# - The point predictions implied by the constructed newdata are checked
#   against pred_table$mu0/mu1 before drawing. If the two disagree the village
#   ordering has drifted and every downstream ribbon would be silently wrong,
#   so this errors rather than warns.
# - Cholesky is tried first and falls back to a clamped eigendecomposition,
#   because a predictive covariance can come back positive semi-definite.
#
# Requirements:
# - Packages: spaMM, dplyr, sf, tibble, purrr
#
# Returns:
# - A stacked tibble with the same columns as build_prediction_table(), with
#   replicate running 1..B.
#--------------------------------------

predictive_draw_tables <- function(spamm_model_fit,
                                   pred_table,
                                   grappe_dat,
                                   B = 200,
                                   phase_years = 0.5,
                                   global_crs = 4326,
                                   proj_crs = 32631,
                                   seed = 20260812) {

  arm_levels <- levels(grappe_dat$isTreated_child)

  # newdata in EXACTLY pred_table's village order, both arms stacked:
  # rows 1..n are Control, rows n+1..2n are Treatment.
  base <- pred_table %>%
    dplyr::select(site_code, grappe, csi, district, lat, lon, n_doses)

  coords_proj <- base %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = global_crs) %>%
    sf::st_transform(crs = proj_crs)

  newd_one <- base %>%
    mutate(lon           = sf::st_coordinates(coords_proj)[, 1],
           lat           = sf::st_coordinates(coords_proj)[, 2],
           district      = as.factor(district),
           site_code     = as.character(site_code),
           pt_year_1to59 = 1000)

  newd <- purrr::map_dfr(c("Control", "Treatment"), function(arm) {
    newd_one %>%
      mutate(isTreated_child     = factor(arm, levels = arm_levels),
             isTreated_child_num = as.integer(arm == "Treatment"))
  })

  n <- nrow(pred_table)

  message("Computing joint predictive covariance over ", nrow(newd),
          " rows (villages x arms); this is the slow step.")

  mu_hat <- predict(spamm_model_fit, newdata = newd, type = "response",
                    variances = list(linPred = TRUE, cov = TRUE))
  V      <- attr(mu_hat, "predVar")
  mu_hat <- as.numeric(mu_hat)

  # Ordering guard: these must be the same predictions build_prediction_table()
  # produced, or the draws are attached to the wrong villages.
  stopifnot(
    max(abs(mu_hat[seq_len(n)]              - pred_table$mu0)) < 1e-6,
    max(abs(mu_hat[n + seq_len(n)]          - pred_table$mu1)) < 1e-6
  )

  eta_hat <- log(mu_hat)

  # Matrix square root of the predictive covariance -----
  V <- (V + t(V)) / 2  # enforce exact symmetry before factorising
  L <- tryCatch(t(chol(V)), error = function(e) {
    message("Cholesky failed (covariance is only PSD); ",
            "falling back to a clamped eigendecomposition.")
    ev <- eigen(V, symmetric = TRUE)
    ev$vectors %*% diag(sqrt(pmax(ev$values, 0)))
  })

  set.seed(seed)
  z   <- matrix(stats::rnorm(nrow(V) * B), nrow = nrow(V), ncol = B)
  eta <- as.vector(eta_hat) + (L %*% z)   # nrow(V) x B, recycled by column
  mu  <- exp(eta)

  purrr::map_dfr(seq_len(B), function(b) {
    mu0_b <- mu[seq_len(n), b]
    mu1_b <- mu[n + seq_len(n), b]
    pred_table %>%
      mutate(replicate      = b,
             mu0            = mu0_b,
             mu1            = mu1_b,
             rd             = mu0_b - mu1_b,
             # observed person-time is the exposure denominator, matching
             # build_prediction_table(); phase_years no longer enters.
             deaths_averted = (mu0_b - mu1_b) * person_time / 1000)
  })
}


#--------------------------------------
# degradation_tables()
#
# Recombine bootstrap replicates so that villages are RANKED by the replicate's
# noisy draw but SCORED against a fixed reference surface.
#
# WHY THIS IS NEEDED: computing a targeting curve entirely within one replicate
# makes the `rd` rule select on the same quantity it accumulates, which is a
# winner's curse. Villages whose rd was perturbed upward move to the front of
# the queue and are then credited with that inflated rd, so the normalised
# curve gets steeper. Injecting pure noise into rd (leaving the surface
# untouched) lifts the rd curve at 30% coverage from 0.574 to 0.694 at sd = 2,
# while leaving mort_grappe flat at ~0.487 - the mortality rules rank on mu0
# and score on rd, so selection acts on a different variable and they are
# unaffected. The consequence is that the within-replicate band for `rd` floats
# entirely above the point-estimate curve, and because that inflated envelope
# is the denominator of the efficiency ratio, the efficiency band is biased
# downward too.
#
# Ranking on the draw and scoring on the reference surface removes the
# selection: a noisier draw can only produce a WORSE ordering, so replicate
# curves sit at or below the point-estimate curve and the band means "what
# estimation error costs you". Every band in 07-targeting-analysis.qmd is built
# this way.
#
# Arguments:
# - boot_tables: stacked replicate tables from predictive_draw_tables().
# - pred_table:  the point-estimate table; supplies the reference surface.
# - score:       "curve" keeps rd from the draw (so the `rd` RULE still ranks
#                on the draw) and takes deaths_averted from the reference, for
#                targeting curves. "concordance" instead takes rd from the
#                reference, so the rule's ordering is judged against the
#                reference effect ordering.
#
# Requirements:
# - Packages: dplyr
#
# Returns:
# - A table of the same shape as boot_tables.
#--------------------------------------

degradation_tables <- function(boot_tables, pred_table,
                               score = c("curve", "concordance")) {

  score <- match.arg(score)

  # Join on site_code, never positionally.
  ref <- pred_table %>%
    dplyr::select(site_code, rd_ref = rd, deaths_averted_ref = deaths_averted)

  out <- boot_tables %>%
    dplyr::left_join(ref, by = "site_code")

  if (any(is.na(out$rd_ref))) {
    stop("Some replicate villages have no match in the reference table; ",
         "the site_code join is wrong.")
  }

  if (score == "curve") {
    # mu0 and rd stay as drawn (they are the ranking signals); only the benefit
    # accrued is taken from the reference surface. The normalising total is
    # then identical across replicates, so the curves are directly comparable.
    out <- out %>% mutate(deaths_averted = deaths_averted_ref)
  } else {
    # rd becomes the yardstick the rule's ordering is judged against, so it
    # comes from the reference; mu0 stays as drawn.
    out <- out %>% mutate(rd = rd_ref, deaths_averted = deaths_averted_ref)
  }

  out %>% dplyr::select(-rd_ref, -deaths_averted_ref)
}
