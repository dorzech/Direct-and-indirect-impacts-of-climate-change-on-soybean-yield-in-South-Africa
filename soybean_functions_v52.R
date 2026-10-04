palSeq <- function(n) hcl.colors(n, "Viridis")

palDiv <- function(n) hcl.colors(n, "Vik")

symBreaks <- function(x, n = 10) {
  lim <- max(abs(x), na.rm = TRUE)
  seq(-lim, lim, length.out = n + 1)
}

commonBreaks <- function(..., n = 10) {
  rng <- range(c(...), na.rm = TRUE)
  seq(rng[1], rng[2], length.out = n + 1)
}

plotCategorical <- function(geom, classes, cols, title) {
  cl <- factor(classes, levels = names(cols))
  op <- par(mar = c(1, 1, 3, 1))
  on.exit(par(op))
  plot(sf::st_geometry(geom), col = cols[as.character(cl)], border = NA, main = title)
  legend("bottomleft", legend = names(cols), fill = cols, border = NA,
         bty = "n", cex = 0.68)
}

vifOf <- function(df) {
  df <- as.data.frame(df)
  v <- vapply(seq_along(df), function(j) {
    1 / (1 - summary(lm(df[[j]] ~ ., data = df[-j]))$r.squared)
  }, numeric(1))
  names(v) <- names(df)
  v
}


# climate

readBioclimDir <- function(dir, predictors, crop_ext = sdm_ext) {
  files <- list.files(dir, pattern = "\\.(tif|bil)$",
                      full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) stop("No rasters in ", dir)

  if (length(files) == 1) {
    r <- terra::rast(files)
  } else {
    # future files end in bi70<band>, current files in <band>
    fn <- basename(files)
    fut <- regmatches(fn, regexec("bi([0-9]{2})([0-9]{1,2})\\.(tif|bil)$", fn))
    bnum <- suppressWarnings(vapply(seq_along(fn), function(k) {
      m <- fut[[k]]
      if (length(m) == 4) as.integer(m[3])
      else as.integer(sub(".*[^0-9]([0-9]+)\\.(tif|bil)$", "\\1", fn[k]))
    }, integer(1)))
    if (anyNA(bnum) || !setequal(bnum, 1:19)) stop("Cannot read band numbers in ", dir)
    r <- terra::rast(files[order(bnum)])
  }

  if (terra::nlyr(r) != 19) stop("Expected 19 layers in ", dir)
  names(r) <- paste0("bio", 1:19)
  terra::crop(r[[predictors]], crop_ext)
}

loadFutureCmip5 <- function(code, template, root = CLIM_FUTURE_ROOT,
                            predictors = predictor_vars) {
  r <- readBioclimDir(file.path(root, code), predictors)
  if (!isTRUE(terra::compareGeom(r, template, stopOnError = FALSE))) {
    r <- terra::resample(r, template, method = "bilinear")
  }
  r
}


# GBIF

firstOf <- function(df, v) if (v %in% names(df)) as.character(df[[v]][1]) else NA_character_

normaliseBackbone <- function(bb, input_names) {
  # column names differ between rgbif versions
  if (!"verbatim_name" %in% names(bb)) {
    alt <- intersect(c("verbatim_scientificName", "verbatim_scientificname"), names(bb))
    bb$verbatim_name <- if (length(alt)) bb[[alt[1]]] else input_names[bb$verbatim_index]
  }
  for (v in c("usageKey", "scientificName", "rank", "status", "matchType", "confidence",
              "kingdom", "phylum", "class", "order", "family")) {
    if (!v %in% names(bb)) bb[[v]] <- NA
  }
  tibble(verbatim_name = as.character(bb$verbatim_name),
         usageKey = suppressWarnings(as.numeric(as.character(bb$usageKey))),
         scientificName = as.character(bb$scientificName),
         rank = toupper(as.character(bb$rank)),
         status = as.character(bb$status),
         matchType = toupper(as.character(bb$matchType)),
         confidence = suppressWarnings(as.numeric(bb$confidence)),
         kingdom = as.character(bb$kingdom), phylum = as.character(bb$phylum),
         class = as.character(bb$class), order = as.character(bb$order),
         family = as.character(bb$family))
}

lookupGenus <- function(g, fam) {
  r1 <- tryCatch(
    normaliseBackbone(name_backbone_checklist(
      data.frame(name = g, rank = "GENUS", kingdom = "Animalia",
                 order = "Hymenoptera", family = fam)), g),
    error = function(e) NULL)
  if (!is.null(r1) && isTRUE(r1$rank[1] == "GENUS") &&
      isTRUE(r1$family[1] %in% BEE_FAMILIES) && !is.na(r1$usageKey[1])) {
    return(r1 %>% mutate(source = "match with family"))
  }

  r2 <- tryCatch(name_lookup(query = g, rank = "GENUS", datasetKey = GBIF_BACKBONE_KEY,
                             limit = 50)$data,
                 error = function(e) NULL)
  if (!is.null(r2) && nrow(r2)) {
    if (!"canonicalName" %in% names(r2)) r2$canonicalName <- sub(" .*", "", r2$scientificName)
    if (!"family" %in% names(r2)) r2$family <- NA_character_
    if (!"taxonomicStatus" %in% names(r2)) r2$taxonomicStatus <- NA_character_
    r2 <- r2 %>% filter(tolower(canonicalName) == tolower(g), family %in% BEE_FAMILIES)
    if (nrow(r2)) {
      acc <- r2 %>% arrange(desc(taxonomicStatus %in% "ACCEPTED")) %>% slice(1)
      if (!is.na(firstOf(acc, "acceptedKey")) && !(acc$taxonomicStatus %in% "ACCEPTED")) {
        key <- firstOf(acc, "acceptedKey")
      } else if (!is.na(firstOf(acc, "nubKey"))) {
        key <- firstOf(acc, "nubKey")
      } else {
        key <- firstOf(acc, "key")
      }
      return(tibble(verbatim_name = g, usageKey = as.numeric(key),
                    scientificName = firstOf(acc, "scientificName"), rank = "GENUS",
                    status = firstOf(acc, "taxonomicStatus"), matchType = "LOOKUP",
                    confidence = NA_real_, kingdom = firstOf(acc, "kingdom"),
                    phylum = firstOf(acc, "phylum"), class = firstOf(acc, "class"),
                    order = firstOf(acc, "order"), family = firstOf(acc, "family"),
                    source = "name_lookup"))
    }
  }

  tibble(verbatim_name = g, usageKey = NA_real_, scientificName = NA_character_,
         rank = NA_character_, status = NA_character_, matchType = "NONE",
         confidence = NA_real_, kingdom = NA_character_, phylum = NA_character_,
         class = NA_character_, order = NA_character_, family = NA_character_,
         source = "unresolved")
}

downloadKeys <- function(keys, out_path, meta_path) {
  if (!all(nzchar(Sys.getenv(c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL"))))) {
    stop("GBIF credentials not set")
  }
  dl_key <- occ_download(
    pred_in("taxonKey", keys),
    pred_in("basisOfRecord", ALLOWED_BASIS),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred_gte("year", YEAR_MIN),
    pred_within(SDM_WKT),
    format = "SIMPLE_CSV"
  )
  occ_download_wait(dl_key)
  dl_meta <- occ_download_meta(dl_key)
  cat("DOI:", dl_meta$doi, "\n")
  cat("Records:", dl_meta$totalRecords, "\n")

  raw_all <- occ_download_import(occ_download_get(dl_key, overwrite = TRUE))
  recs <- raw_all %>%
    filter(family %in% BEE_FAMILIES, genus != "Apis") %>%
    transmute(key = gbifID, datasetKey, basisOfRecord, scientificName,
              family, genus, species, decimalLatitude, decimalLongitude,
              coordinateUncertaintyInMeters, countryCode, year)
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(recs, out_path)
  saveRDS(list(download_key = dl_key, doi = dl_meta$doi, created = dl_meta$created,
               total_raw = dl_meta$totalRecords, n_genus_keys = length(keys), keys = keys,
               wkt = SDM_WKT, basis = ALLOWED_BASIS, year_min = YEAR_MIN), meta_path)
  invisible(recs)
}

cleanRecords <- function(raw, clim_template) {
  occ <- raw %>%
    filter(basisOfRecord %in% ALLOWED_BASIS,
           !is.na(decimalLongitude), !is.na(decimalLatitude),
           !is.na(species), species != "")
  n_basis <- nrow(occ)

  occ <- occ %>%
    filter(!is.na(year), year >= YEAR_MIN,
           is.na(coordinateUncertaintyInMeters) |
             coordinateUncertaintyInMeters <= MAX_COORD_UNC_M)
  n_year <- nrow(occ)

  # seas test uses the coastline shipped with CoordinateCleaner
  data("buffland", package = "CoordinateCleaner", envir = environment())
  ref_sea <- buffland
  if (inherits(ref_sea, "PackedSpatVector")) ref_sea <- terra::unwrap(ref_sea)
  if (!inherits(ref_sea, "SpatVector")) ref_sea <- terra::vect(ref_sea)

  cc <- clean_coordinates(
    x = as.data.frame(occ),
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    species = "species",
    tests = c("capitals", "centroids", "equal", "gbif",
              "institutions", "seas", "zeros", "duplicates"),
    seas_ref = ref_sea,
    value = "spatialvalid",
    verbose = FALSE
  )
  keep <- as.logical(cc$.summary)
  keep[is.na(keep)] <- FALSE
  occ <- occ[keep, , drop = FALSE]
  n_cc <- nrow(occ)

  occ <- occ %>%
    mutate(.cell = terra::cellFromXY(clim_template[[1]],
                                     cbind(decimalLongitude, decimalLatitude))) %>%
    filter(!is.na(.cell)) %>%
    distinct(species, .cell, .keep_all = TRUE) %>%
    select(-.cell)

  attr(occ, "counts") <- c(after_basis = n_basis, after_year = n_year,
                           after_cc = n_cc, after_thin = nrow(occ))
  occ
}


# SDM

aucSafe <- function(obs, pred) {
  if (length(unique(obs)) < 2 || all(is.na(pred))) return(NA_real_)
  as.numeric(auc(roc(obs, pred, quiet = TRUE)))
}

tssFromScores <- function(obs, pred) {
  if (all(is.na(pred)) || length(unique(obs)) < 2) {
    return(list(thr = 0.5, tss = NA_real_))
  }
  thr_grid <- seq(0.01, 0.99, by = 0.01)
  tss_vec <- vapply(thr_grid, function(t) {
    pp <- as.integer(pred >= t)
    tp <- sum(pp == 1 & obs == 1, na.rm = TRUE)
    fn <- sum(pp == 0 & obs == 1, na.rm = TRUE)
    tn <- sum(pp == 0 & obs == 0, na.rm = TRUE)
    fp <- sum(pp == 1 & obs == 0, na.rm = TRUE)
    sens <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else 0
    sens + spec - 1
  }, numeric(1))
  i <- which.max(tss_vec)
  list(thr = thr_grid[i], tss = tss_vec[i])
}

sensSpec <- function(obs, pred, thr) {
  if (all(is.na(pred))) return(c(sensitivity = NA_real_, specificity = NA_real_))
  pp <- as.integer(pred >= thr)
  c(sensitivity = sum(pp == 1 & obs == 1, na.rm = TRUE) / sum(obs == 1),
    specificity = sum(pp == 0 & obs == 0, na.rm = TRUE) / sum(obs == 0))
}

predictMaxnetRaster <- function(stk, model) {
  terra::predict(stk, model, na.rm = TRUE, fun = function(model, x, ...) {
    df <- as.data.frame(x)
    names(df) <- predictor_vars
    out <- tryCatch(predict(model, newdata = df, type = "cloglog"),
                    error = function(e) rep(NA_real_, nrow(df)))
    as.numeric(out)
  })
}

sampleBackground <- function(presence_df, target_df, clim_template, n_bg, calib, seed) {
  set.seed(seed)
  tmpl <- terra::mask(clim_template[[1]], calib)
  pres_cells <- terra::cellFromXY(tmpl, cbind(presence_df$decimalLongitude,
                                              presence_df$decimalLatitude))
  bg_type <- "target_group"
  xy <- NULL
  tg_cells <- terra::cellFromXY(tmpl, cbind(target_df$decimalLongitude,
                                            target_df$decimalLatitude))
  tg_cells <- setdiff(unique(na.omit(tg_cells)), pres_cells)
  if (length(tg_cells)) tg_cells <- tg_cells[!is.na(tmpl[tg_cells][, 1])]
  if (length(tg_cells) >= MIN_BG_CELLS) {
    keep <- if (length(tg_cells) > n_bg) sample(tg_cells, n_bg) else tg_cells
    xy <- terra::xyFromCell(tmpl, keep)
  } else {
    bg_type <- "random_within_calibration"
  }
  if (is.null(xy)) {
    tmpl[pres_cells] <- NA
    xy <- as.matrix(terra::spatSample(tmpl, size = n_bg, method = "random",
                                      xy = TRUE, na.rm = TRUE)[, c("x", "y")])
  }
  if (nrow(xy) < MIN_RECORDS) {
    stop(sprintf("only %d background cells (%s)", nrow(xy), bg_type))
  }
  tibble(decimalLongitude = as.numeric(xy[, 1]),
         decimalLatitude = as.numeric(xy[, 2]),
         n_bg_realised = as.integer(nrow(xy)),
         bg_type = bg_type)
}

makeFolds <- function(pa, k, seed) {
  set.seed(seed)
  f <- integer(length(pa))
  f[pa == 1] <- sample(rep(seq_len(k), length.out = sum(pa == 1)))
  f[pa == 0] <- sample(rep(seq_len(k), length.out = sum(pa == 0)))
  f
}

fitAlgorithms <- function(train) {
  n1 <- min(sum(train$pa == 1), sum(train$pa == 0))
  y <- factor(train$pa, levels = c(0, 1))
  gam_formula <- as.formula(
    paste("pa ~", paste0("s(", predictor_vars, ", k = ", GAM_K, ")", collapse = " + "))
  )
  # equal total weight for presences and background
  gam_w <- ifelse(train$pa == 1, 1, sum(train$pa == 1) / sum(train$pa == 0))
  list(
    rf = randomForest(x = train[, predictor_vars], y = y, ntree = 500,
                      sampsize = c(n1, n1), strata = y),
    gam = suppressWarnings(gam(gam_formula, data = train, family = binomial(),
                               method = "REML", weights = gam_w, select = TRUE)),
    mx = tryCatch(
      maxnet(p = train$pa, data = train[, predictor_vars],
             f = maxnet.formula(p = train$pa, data = train[, predictor_vars],
                                classes = "default")),
      error = function(e) NULL
    )
  )
}

predictAlgorithms <- function(mods, newdata) {
  list(
    rf = as.numeric(predict(mods$rf, newdata[, predictor_vars], type = "prob")[, "1"]),
    gam = as.numeric(predict(mods$gam, newdata, type = "response")),
    mx = if (is.null(mods$mx)) {
      rep(NA_real_, nrow(newdata))
    } else {
      as.numeric(predict(mods$mx, newdata[, predictor_vars], type = "cloglog"))
    }
  )
}

makeSpatialFolds <- function(dat, k, seed) {
  pts <- sf::st_as_sf(as.data.frame(dat[, c("decimalLongitude", "decimalLatitude", "pa")]),
                      coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  cv <- blockCV::cv_spatial(x = pts, column = "pa", k = k, size = BLOCK_SIZE_M,
                            selection = "random", iteration = 200, seed = seed,
                            progress = FALSE, report = FALSE, plot = FALSE)
  f <- cv$folds_ids
  attr(f, "block_m") <- BLOCK_SIZE_M
  f
}

foldsValid <- function(pa, folds) {
  !is.null(folds) && length(unique(folds)) >= 2 &&
    all(tapply(pa, folds, function(v) sum(v == 1)) >= MIN_TEST_PRES) &&
    all(tapply(pa, folds, function(v) sum(v == 0)) >= 1)
}

chooseFolds <- function(dat, cv_type, seed) {
  if (cv_type == "spatial") {
    for (k in SPATIAL_K_TRY) {
      f <- tryCatch(makeSpatialFolds(dat, k = k, seed = seed), error = function(e) NULL)
      if (foldsValid(dat$pa, f)) {
        return(list(folds = f, cv_used = sprintf("spatial_k%d", k),
                    block_m = attr(f, "block_m")))
      }
    }
    message("  spatial folds failed, using random folds")
    return(list(folds = makeFolds(dat$pa, k = K_FOLDS, seed = seed),
                cv_used = "random_fallback", block_m = NA_real_))
  }
  list(folds = makeFolds(dat$pa, k = K_FOLDS, seed = seed),
       cv_used = "random", block_m = NA_real_)
}

evaluateKfold <- function(dat, seed, cv_type) {
  # thresholds (max TSS) come from the pooled out-of-fold predictions
  fc <- chooseFolds(dat, cv_type, seed)
  folds <- fc$folds
  algs <- c("rf", "gam", "mx")
  oof <- matrix(NA_real_, nrow(dat), length(algs), dimnames = list(NULL, algs))
  for (f in sort(unique(folds))) {
    test_idx <- which(folds == f)
    mods <- fitAlgorithms(dat[-test_idx, ])
    preds <- predictAlgorithms(mods, dat[test_idx, ])
    for (a in algs) oof[test_idx, a] <- preds[[a]]
  }
  mx_ok <- !is.na(oof[, "mx"])
  ens <- (oof[, "rf"] + oof[, "gam"] + ifelse(mx_ok, oof[, "mx"], 0)) / (2 + mx_ok)
  oof_all <- cbind(oof, ensemble = ens)
  thr <- vapply(colnames(oof_all), function(a) tssFromScores(dat$pa, oof_all[, a])$thr,
                numeric(1))
  rows <- list()
  for (f in sort(unique(folds))) {
    test_idx <- which(folds == f)
    for (a in colnames(oof_all)) {
      pred <- oof_all[test_idx, a]
      ss <- sensSpec(dat$pa[test_idx], pred, thr[[a]])
      rows[[length(rows) + 1]] <- tibble(
        fold = f, algorithm = a,
        n_test_pres = sum(dat$pa[test_idx] == 1), n_test_bg = sum(dat$pa[test_idx] == 0),
        auc = aucSafe(dat$pa[test_idx], pred), threshold = thr[[a]],
        tss = unname(ss[["sensitivity"]] + ss[["specificity"]] - 1),
        sensitivity = ss[["sensitivity"]], specificity = ss[["specificity"]])
    }
  }
  list(eval = bind_rows(rows), oof = oof_all, thresholds = thr, folds = folds,
       k_used = length(unique(folds)), cv_used = fc$cv_used, block_m = fc$block_m,
       mx_folds_failed = sum(tapply(!mx_ok, folds, all)))
}

boyceContinuous <- function(fit, obs, res = 100, window_w = NULL, method = "spearman") {
  # same algorithm as ecospat::ecospat.boyce()
  fit <- fit[!is.na(fit)]
  obs <- obs[!is.na(obs)]
  if (length(obs) < 5 || length(fit) < 20) return(NA_real_)
  mini <- min(fit)
  maxi <- max(fit)
  if (is.null(window_w)) window_w <- (maxi - mini) / 10
  vec_mov <- seq(mini, maxi - window_w, by = (maxi - mini - window_w) / res)
  vec_mov[res + 1] <- vec_mov[res + 1] + 1
  lo <- vec_mov
  hi <- vec_mov + window_w
  pe <- vapply(seq_along(lo), function(i) {
    p <- mean(obs >= lo[i] & obs <= hi[i])
    e <- mean(fit >= lo[i] & fit <= hi[i])
    if (e == 0) NaN else p / e
  }, numeric(1))
  keep <- which(!is.nan(pe))
  pe <- pe[keep]
  pos <- vec_mov[keep]
  if (length(pe) < 3) return(NA_real_)
  r <- seq_along(pe)[pe != c(pe[-1], TRUE)]
  suppressWarnings(stats::cor(pe[r], pos[r], method = method))
}

clampStack <- function(stk, rng) {
  out <- terra::rast(lapply(predictor_vars, function(v)
    terra::clamp(stk[[v]], lower = rng[1, v], upper = rng[2, v], values = TRUE)))
  names(out) <- predictor_vars
  out
}

noveltyByPredictor <- function(dat, clim_fut_list) {
  rng <- sapply(predictor_vars, function(v) range(dat[[v]], na.rm = TRUE))
  out <- lapply(predictor_vars, function(v) {
    terra::app(terra::rast(lapply(clim_fut_list, function(stk)
      stk[[v]] < rng[1, v] | stk[[v]] > rng[2, v])), max, na.rm = TRUE)
  })
  out <- terra::rast(out)
  names(out) <- predictor_vars
  out
}

responseCurves <- function(models, dat, weights, n_pts = 50) {
  med <- vapply(predictor_vars, function(v) median(dat[[v]], na.rm = TRUE), numeric(1))
  map_dfr(predictor_vars, function(v) {
    grid <- seq(min(dat[[v]], na.rm = TRUE), max(dat[[v]], na.rm = TRUE), length.out = n_pts)
    nd <- as.data.frame(matrix(rep(med, each = n_pts), nrow = n_pts,
                               dimnames = list(NULL, predictor_vars)))
    nd[[v]] <- grid
    pr <- predictAlgorithms(models, nd)
    ens <- weights[["rf"]] * pr$rf + weights[["gam"]] * pr$gam +
      weights[["mx"]] * replace_na(pr$mx, 0)
    tibble(predictor = v, value = grid, rf = pr$rf, gam = pr$gam, mx = pr$mx, ensemble = ens)
  })
}

noveltyLayers <- function(dat, clim_fut_list) {
  rng <- sapply(predictor_vars, function(v) range(dat[[v]], na.rm = TRUE))
  nov <- lapply(clim_fut_list, function(stk) {
    Reduce(`|`, lapply(predictor_vars, function(v) stk[[v]] < rng[1, v] | stk[[v]] > rng[2, v]))
  })
  nov_stack <- terra::rast(nov)
  names(nov_stack) <- names(clim_fut_list)
  nov_stack
}

ensureOutputDirs <- function(gcm_names) {
  dirs <- c(file.path("points", "background"),
            file.path("predictions", c("rf", "gam", "mx")),
            file.path("predictions", "ensemble"),
            "evaluation",
            file.path("future", gcm_names, "predictions", "ensemble"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

fitSdm <- function(presences, background, clim_hist, clim_fut_list, label, seed, calib) {
  pres_env <- terra::extract(clim_hist,
                             cbind(presences$decimalLongitude,
                                   presences$decimalLatitude))
  bg_env <- terra::extract(clim_hist,
                           cbind(background$decimalLongitude,
                                 background$decimalLatitude))
  pres <- cbind(presences, pres_env) %>% mutate(pa = 1L)
  bg <- cbind(background[, c("decimalLongitude", "decimalLatitude")], bg_env) %>%
    mutate(pa = 0L)
  dat <- bind_rows(pres, bg) %>% drop_na(all_of(predictor_vars))

  n_pres <- sum(dat$pa == 1)
  if (n_pres < MIN_RECORDS) {
    return(list(skip = TRUE, label = label,
                reason = sprintf("only %d clean presences", n_pres)))
  }

  ev <- evaluateKfold(dat, seed = seed, cv_type = CV_TYPE)
  eval_tbl <- ev$eval %>% mutate(species = label, cv = ev$cv_used, .before = 1)
  eval_random <- NULL
  if (CV_TYPE == "spatial" && RUN_RANDOM_COMPARISON) {
    ev_r <- evaluateKfold(dat, seed = seed, cv_type = "random")
    eval_random <- ev_r$eval %>% mutate(species = label, cv = "random", .before = 1)
  }
  eval_mean <- eval_tbl %>%
    group_by(algorithm) %>%
    summarise(across(c(auc, threshold, tss, sensitivity, specificity),
                     ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  weights <- c(rf = 1, gam = 1, mx = as.numeric(!all(is.na(ev$oof[, "mx"]))))
  weights <- weights / sum(weights)
  ens_oof <- ev$oof[, "ensemble"]
  threshold <- unname(ev$thresholds[["ensemble"]])
  boyce <- boyceContinuous(fit = ens_oof, obs = ens_oof[dat$pa == 1])

  set.seed(seed)
  final <- fitAlgorithms(dat)
  if (is.null(final$mx)) weights[["mx"]] <- 0
  weights <- weights / sum(weights)

  novelty <- noveltyLayers(dat, clim_fut_list)
  novelty_pred <- noveltyByPredictor(dat, clim_fut_list)
  cal_mat <- as.matrix(dat[, predictor_vars])
  cal_mu <- colMeans(cal_mat)
  cal_S <- stats::cov(cal_mat)
  d2_max <- max(stats::mahalanobis(cal_mat, cal_mu, cal_S))
  nt2 <- terra::rast(lapply(clim_fut_list, function(stk)
    terra::app(stk[[predictor_vars]], function(x) stats::mahalanobis(x, cal_mu, cal_S)) / d2_max))
  names(nt2) <- names(clim_fut_list)
  rf_importance <- randomForest::importance(final$rf)[, 1]

  predictStack <- function(stk) {
    p_rf <- terra::predict(stk, final$rf, na.rm = TRUE,
                           fun = function(model, x, ...)
                             predict(model, x, type = "prob")[, "1"])
    p_gam <- terra::predict(stk, final$gam, type = "response", na.rm = TRUE)
    out <- weights[["rf"]] * p_rf + weights[["gam"]] * p_gam
    p_mx <- NULL
    if (!is.null(final$mx)) {
      p_mx <- predictMaxnetRaster(stk, final$mx)
      out <- out + weights[["mx"]] * terra::ifel(is.na(p_mx), 0, p_mx)
    }
    list(rf = p_rf, gam = p_gam, mx = p_mx, ensemble = out)
  }

  hist_pred <- predictStack(clim_hist)
  fut_preds <- lapply(clim_fut_list, predictStack)
  fut_ens <- terra::rast(lapply(fut_preds, `[[`, "ensemble"))
  names(fut_ens) <- names(clim_fut_list)
  cal_rng <- sapply(predictor_vars, function(v) range(dat[[v]], na.rm = TRUE))
  fut_ens_clamped <- terra::rast(lapply(clim_fut_list, function(stk)
    predictStack(clampStack(stk, cal_rng))$ensemble))
  names(fut_ens_clamped) <- names(clim_fut_list)

  rcl <- matrix(c(-Inf, threshold, 0, threshold, Inf, 1), ncol = 3, byrow = TRUE)
  fut_bin <- terra::classify(fut_ens, rcl = rcl)
  names(fut_bin) <- names(clim_fut_list)
  hist_bin <- terra::classify(hist_pred$ensemble, rcl = rcl)
  fut_bin_primary <- terra::classify(terra::app(fut_ens, mean, na.rm = TRUE), rcl = rcl)
  fut_bin_noexp_primary <- fut_bin_primary * hist_bin

  sp_file <- gsub(" ", "_", label)
  for (a in c("rf", "gam", "mx")) {
    if (!is.null(hist_pred[[a]])) {
      terra::writeRaster(hist_pred[[a]],
                         file.path("predictions", a, paste0(sp_file, ".tif")),
                         overwrite = TRUE)
    }
    if (!is.null(final[[a]])) {
      saveRDS(final[[a]], file.path("evaluation", paste0(sp_file, "_", a, "_model.rds")))
    }
  }
  terra::writeRaster(hist_pred$ensemble,
                     file.path("predictions", "ensemble", paste0(sp_file, ".tif")),
                     overwrite = TRUE)
  for (g in names(clim_fut_list)) {
    terra::writeRaster(fut_preds[[g]]$ensemble,
                       file.path("future", g, "predictions", "ensemble", paste0(sp_file, ".tif")),
                       overwrite = TRUE)
  }
  saveRDS(list(species = label, eval_by_fold = eval_tbl, eval_mean = eval_mean,
               weights = weights, threshold = threshold, thresholds = ev$thresholds,
               n_pres = n_pres, cv_used = ev$cv_used, block_m = ev$block_m,
               n_bg_realised = background$n_bg_realised[1],
               bg_type = background$bg_type[1]),
          file.path("evaluation", paste0(sp_file, "_eval.rds")))

  list(skip = FALSE, label = label,
       n_bg_realised = as.integer(background$n_bg_realised[1]),
       bg_type = background$bg_type[1],
       weights = weights,
       threshold = threshold,
       eval_by_fold = eval_tbl,
       ens_hist = hist_pred$ensemble,
       ens_fut = terra::app(fut_ens, mean, na.rm = TRUE),
       ens_fut_clamped = terra::app(fut_ens_clamped, mean, na.rm = TRUE),
       fut_bin_gcm = fut_bin,
       hist_bin = hist_bin,
       fut_bin_noexp_primary = fut_bin_noexp_primary,
       calib = terra::wrap(calib),
       boyce = boyce,
       novelty_by_predictor = novelty_pred,
       nt2_gcm = nt2,
       eval_random = eval_random,
       k_used = ev$k_used,
       mx_folds_failed = ev$mx_folds_failed,
       rf_importance = rf_importance,
       models = final,
       dat = dat,
       novelty = novelty,
       cv_used = ev$cv_used,
       block_m = ev$block_m,
       maxnet_ok = !is.null(final$mx))
}

runSdmSet <- function(taxa_tbl, cleaned_records, clim_hist, clim_fut_list,
                      base_seed = 123L, n_bg = 10000L) {
  ensureOutputDirs(names(clim_fut_list))
  set.seed(base_seed)
  results <- list()
  skipped <- tibble(species = character(), reason = character())

  for (i in seq_len(nrow(taxa_tbl))) {
    sp <- taxa_tbl$species[i]
    pres <- cleaned_records %>% filter(species == sp)

    t0 <- Sys.time()
    message(sprintf("[%d/%d] %s (n=%d) %s", i, nrow(taxa_tbl), sp, nrow(pres),
                    format(t0, "%H:%M:%S")))

    res <- tryCatch({
      pres_v <- terra::vect(as.data.frame(pres[, c("decimalLongitude", "decimalLatitude")]),
                            geom = c("decimalLongitude", "decimalLatitude"),
                            crs = "EPSG:4326")
      calib_sp <- terra::aggregate(terra::buffer(pres_v, width = CALIB_BUFFER_M))
      bg <- sampleBackground(pres,
                             target_df = cleaned_records %>% filter(species != sp),
                             clim_template = clim_hist, n_bg = n_bg,
                             calib = calib_sp, seed = base_seed + i)
      write_csv(bg, file.path("points", "background", paste0(gsub(" ", "_", sp), ".csv")))
      if (EXCLUDE_NO_TARGET_GROUP && bg$bg_type[1] != "target_group") {
        list(skip = TRUE, label = sp,
             reason = sprintf("fewer than %d target-group background cells", MIN_BG_CELLS))
      } else {
        fitSdm(pres, bg, clim_hist, clim_fut_list, label = sp,
               seed = base_seed + 1000L + i, calib = calib_sp)
      }
    }, error = function(e) list(skip = TRUE, label = sp,
                                reason = paste("error:", conditionMessage(e))))
    message(sprintf("  %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

    if (isTRUE(res$skip)) {
      skipped <- bind_rows(skipped, tibble(species = sp, reason = res$reason))
    } else {
      results[[sp]] <- res
    }
  }

  list(results = results, skipped = skipped)
}


# richness and areas

thresholdEach <- function(res, period = c("hist", "fut")) {
  period <- match.arg(period)
  r <- if (period == "hist") res$ens_hist else res$ens_fut
  terra::classify(r, rcl = matrix(c(-Inf, res$threshold, 0, res$threshold, Inf, 1),
                                  ncol = 3, byrow = TRUE))
}

sumMasked <- function(layers) {
  terra::mask(terra::app(terra::rast(layers), sum, na.rm = TRUE), clim_hist[[1]])
}

stackRichness <- function(res_list) {
  list(n = length(res_list),
       hist = sumMasked(lapply(res_list, thresholdEach, period = "hist")),
       fut = sumMasked(lapply(res_list, thresholdEach, period = "fut")))
}

stackRichnessField <- function(res_list, field) {
  list(n = length(res_list),
       hist = sumMasked(lapply(res_list, thresholdEach, period = "hist")),
       fut = sumMasked(lapply(res_list, `[[`, field)))
}

stackRichnessByGcm <- function(res_list, gcm_names) {
  hist_r <- sumMasked(lapply(res_list, thresholdEach, period = "hist"))
  lapply(setNames(gcm_names, gcm_names), function(g) {
    list(n = length(res_list), hist = hist_r,
         fut = sumMasked(lapply(res_list, function(r) r$fut_bin_gcm[[g]])))
  })
}

areaKm2 <- function(r, csz) as.numeric(terra::global(r * csz, "sum", na.rm = TRUE))

validAreaKm2 <- function(r, csz) areaKm2(terra::ifel(is.na(r), NA, 1), csz)

sharePct <- function(part, whole) if (whole > 0) 100 * part / whole else 0

areaWithin <- function(bin, mask_v = NULL) {
  b <- if (is.null(mask_v)) bin else terra::mask(bin, mask_v)
  areaKm2(b, terra::cellSize(bin, unit = "km"))
}

speciesAreaChange <- function(res, mask_v = NULL) {
  a_h <- areaWithin(thresholdEach(res, "hist"), mask_v)
  a_f <- areaWithin(thresholdEach(res, "fut"), mask_v)
  tibble(species = res$label, area_hist_km2 = a_h, area_fut_km2 = a_f,
         change_km2 = a_f - a_h, pct_change = 100 * (a_f - a_h) / a_h)
}


# catchments

catchmentRichness <- function(rich, k_eq = 1) {
  # exact = TRUE weights cells by the share of the polygon they cover
  h <- terra::extract(rich$hist, catch_v, fun = mean, na.rm = TRUE, exact = TRUE, ID = FALSE)[, 1]
  f <- terra::extract(rich$fut, catch_v, fun = mean, na.rm = TRUE, exact = TRUE, ID = FALSE)[, 1]
  loss <- round(h - f, 9)  # removes floating-point error from the weighted means
  tibble(rich_hist = h, rich_fut = f, rich_change = f - h,
         pct_retained = ifelse(h > 0, 100 * f / h, NA_real_),
         bee_decline = if (k_eq > 0) loss >= k_eq else loss > 0)
}

overlapShare <- function(poll, yield_gain) {
  ok <- !is.na(yield_gain) & !is.na(poll$bee_decline)
  with_decline <- sum(map_aea$area_km2[which(ok & yield_gain & poll$bee_decline)], na.rm = TRUE)
  without_decline <- sum(map_aea$area_km2[which(ok & yield_gain & !poll$bee_decline)], na.rm = TRUE)
  tibble(area_gain_km2 = with_decline + without_decline, area_gain_with_decline_km2 = with_decline,
         overlap_share_pct = 100 * with_decline / (with_decline + without_decline))
}

classifyRisk <- function(pct_change, bee_decline, decline_pct = YIELD_DECLINE_PCT) {
  yield_drop <- (pct_change < -decline_pct) & !is.na(pct_change)
  case_when(
    is.na(bee_decline) ~ "Not assessed",
    is.na(pct_change) ~ FLOOR_LABEL,
    yield_drop & bee_decline ~ "Both decline",
    yield_drop & !bee_decline ~ "Yield decline only",
    !yield_drop & bee_decline ~ "Richness decline only",
    TRUE ~ "No projected decline in either"
  )
}

classifyGain <- function(pct_change, yield_gain, bee_decline) {
  case_when(
    is.na(bee_decline) ~ "Not assessed",
    is.na(pct_change) ~ FLOOR_LABEL,
    !yield_gain ~ "No projected gain",
    yield_gain & bee_decline ~ "Gain with richness decline",
    TRUE ~ "Gain without richness decline"
  )
}

areaOf <- function(tbl, cls) {
  v <- tbl$area_km2[tbl$gain_status == cls]
  if (length(v) == 0) 0 else v
}
