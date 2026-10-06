# =============================================================================
# Datierungsteam: Shiny-App fuer das Seminar "Datierungsmethoden"
# 14C-Kalibration + einfache Bayes'sche Verfeinerung (Stratigrafie, TPQ/TAQ,
# Dendro-Fenster) mit Annahmen-Schaltern, Szenarien, Altholz-Check,
# Messwert-Umrechnung (F14C), Wiggle-Matching, SPD, OxCal-Export
#
# Start:   install.packages(c("shiny", "rcarbon"))
#          shiny::runApp("datierung_app")
#
# Die Beispieldaten sind FIKTIV und nur zur Demonstration gedacht.
# Das Modell ist ein didaktisches Monte-Carlo-Modell, kein Ersatz fuer OxCal.
# =============================================================================

library(shiny)

# ---- IntCal20 laden ---------------------------------------------------------
# Reihenfolge: (1) Datei intcal20.14c neben app.R, (2) Paket rcarbon,
# (3) Download von intcal.org (wird neben app.R gespeichert).
parse_curve <- function(df) {
  df <- as.data.frame(df)[, 1:3]
  names(df) <- c("calBP", "c14", "err")
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
  df <- df[stats::complete.cases(df), ]
  if (nrow(df) < 1000) stop("Kurve unvollstaendig")
  df
}

load_curve <- function() {
  ic <- NULL
  f <- file.path(getwd(), "intcal20.14c")

  # (1) lokale Datei
  if (file.exists(f)) {
    ic <- tryCatch(parse_curve(read.csv(f, comment.char = "#", header = FALSE)), error = function(e) NULL)
  }
  # (2) rcarbon (nur lokal; im Browser-Betrieb mit Shinylive/webR nicht verfuegbar)
  is_wasm <- grepl("emscripten", R.version$platform, fixed = TRUE)
  RC <- "rcarbon"   # Name bewusst als Variable, damit Shinylive das Paket nicht einbinden will
  if (is.null(ic) && !is_wasm && requireNamespace(RC, quietly = TRUE)) {
    for (nm in c("intcal20", "IntCal20")) {
      x <- tryCatch(getExportedValue(RC, nm), error = function(e) NULL)
      if (!is.null(x)) {
        ic <- tryCatch(parse_curve(x), error = function(e) NULL)
        if (!is.null(ic)) break
      }
    }
  }
  # (3) Download
  if (is.null(ic) && !is_wasm) {
    ok <- tryCatch({
      download.file("https://intcal.org/curves/intcal20.14c", f, quiet = TRUE, mode = "wb")
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (ok && file.exists(f)) {
      ic <- tryCatch(parse_curve(read.csv(f, comment.char = "#", header = FALSE)), error = function(e) NULL)
    }
  }
  if (is.null(ic)) {
    stop(paste0(
      "IntCal20 konnte nicht geladen werden. Optionen: (a) install.packages('rcarbon') und App neu starten, ",
      "(b) die Datei https://intcal.org/curves/intcal20.14c herunterladen und als 'intcal20.14c' in den App-Ordner legen ",
      "(bei der Web-Version muss die Datei im Ordner app/ des Repositories liegen). ",
      "Arbeitsverzeichnis: ", getwd()))
  }
  ic <- ic[ic$calBP <= 15000, ]
  ic$year <- 1950 - ic$calBP          # Kalenderjahr, negativ = v. Chr. (1 Jahr Unschaerfe: kein Jahr 0)
  ic[order(ic$year), ]
}
CURVE <- load_curve()

# ---- Jahre: Eingabe (astronomisch ODER "480 v. Chr." / "101 n. Chr.") und Anzeige ----
# Interne Zaehlung wie in OxCal: astronomisch, Jahr 0 = 1 v. Chr., -479 = 480 v. Chr.
parse_year1 <- function(x) {
  if (length(x) != 1 || is.na(x)) return(NA_real_)
  if (is.numeric(x)) return(as.numeric(x))
  s <- trimws(as.character(x))
  if (!nzchar(s)) return(NA_real_)
  m <- regmatches(s, regexec("^(-?[0-9]+)\\s*(v\\.?\\s*Chr\\.?|BCE?|n\\.?\\s*Chr\\.?|CE|AD)$", s, ignore.case = TRUE))[[1]]
  if (length(m) == 3) {
    n <- abs(as.numeric(m[2]))
    era <- tolower(gsub("[ .]", "", m[3]))
    if (era %in% c("vchr", "bc", "bce")) return(1 - n)
    return(n)
  }
  suppressWarnings(as.numeric(s))
}
parse_year <- function(x) vapply(x, parse_year1, numeric(1), USE.NAMES = FALSE)

fmt_year <- function(y, mode = "astro") {
  if (identical(mode, "bcad")) ifelse(y <= 0, paste0(1 - y, " v. Chr."), paste0(y, " n. Chr."))
  else as.character(y)
}
xlab_years <- function(mode) if (identical(mode, "bcad")) "Kalenderjahr" else "Kalenderjahr (astronomisch: -479 = 480 v. Chr.)"
axis_years <- function(rng, mode = "astro") {
  if (identical(mode, "bcad")) {
    # Ticks auf runden historischen Zahlen (z. B. 500 v. Chr.), intern astronomisch gesetzt
    h_rng <- ifelse(rng <= 0, rng - 1, rng)
    h <- pretty(h_rng, n = 6)
    h <- h[h != 0]
    at <- ifelse(h < 0, h + 1, h)
    keep <- at >= rng[1] & at <= rng[2]
    h <- h[keep]; at <- at[keep]
    labels <- ifelse(h < 0, paste0(-h, " v. Chr."), paste0(h, " n. Chr."))
  } else {
    at <- pretty(rng, n = 6)
    at <- at[at >= rng[1] & at <= rng[2]]
    labels <- as.character(at)
  }
  axis(1, at = at, labels = labels, cex.axis = 0.9)
}

# ---- Hilfsfunktionen --------------------------------------------------------
read_txt <- function(txt) {
  txt <- trimws(txt)
  if (!nzchar(txt)) return(NULL)
  first <- strsplit(txt, "\n")[[1]][1]
  sep <- if (grepl(";", first)) ";" else if (grepl("\t", first)) "\t" else ","
  read.table(text = txt, header = TRUE, sep = sep, strip.white = TRUE,
             stringsAsFactors = FALSE, comment.char = "",
             na.strings = c("", "NA"), fill = TRUE)
}

curve_at <- function(yrs) {
  list(mu = approx(CURVE$year, CURVE$c14, xout = yrs, rule = 2)$y,
       sc = approx(CURVE$year, CURVE$err, xout = yrs, rule = 2)$y)
}

# Kalibration: Wahrscheinlichkeit pro Kalenderjahr (Summe = 1)
dens_c14 <- function(age, err, yrs) {
  cv <- curve_at(yrs)
  s2 <- err^2 + cv$sc^2
  d <- exp(-0.5 * (age - cv$mu)^2 / s2) / sqrt(s2)
  d / sum(d)
}

smooth5 <- function(d, k = 5) {
  s <- as.numeric(stats::filter(d, rep(1 / k, k), sides = 2))
  s[is.na(s)] <- 0
  if (sum(s) > 0) s / sum(s) else s
}

# Glaettung mit Randauffuellung (fuer Abstandsverteilungen, die bei 0 beginnen)
smooth_edge <- function(d, k = 5) {
  n <- length(d)
  if (n < 6) return(d)
  d2 <- c(rev(d[1:2]), d, rev(d[(n - 1):n]))
  s <- as.numeric(stats::filter(d2, rep(1 / k, k), sides = 2))[3:(n + 2)]
  s[is.na(s)] <- 0
  if (sum(s) > 0) s / sum(s) else s
}

# 95,4 %-HPD-Bereiche als Text
hpd <- function(d, yrs, p = 0.954, gap = 3, mode = "astro") {
  if (sum(d) <= 0) return(NA_character_)
  d <- d / sum(d)
  o <- order(d, decreasing = TRUE)
  keep <- o[seq_len(which(cumsum(d[o]) >= p)[1])]
  y <- sort(yrs[keep])
  br <- c(0, which(diff(y) > gap), length(y))
  segs <- vapply(seq_len(length(br) - 1), function(j) {
    a <- y[br[j] + 1]; b <- y[br[j + 1]]
    if (a == b) fmt_year(a, mode) else paste0(fmt_year(a, mode), " bis ", fmt_year(b, mode))
  }, character(1))
  paste(segs, collapse = "; ")
}

is_c14 <- function(typ) toupper(trimws(typ)) %in% c("14C", "C14")

# Zeitraster zu einem einzelnen 14C-Alter (fuer das Umrechnungs-Tab)
age_grid <- function(age, err) {
  s <- sqrt(err^2 + CURVE$err^2)
  idx <- which(abs(age - CURVE$c14) < 4 * s)
  if (length(idx) == 0) stop("Alter ausserhalb der Kalibrationskurve")
  seq(floor(min(CURVE$year[idx])) - 25, ceiling(max(CURVE$year[idx])) + 25)
}

# Zeitraster aus den Daten ableiten
make_grid <- function(smp, lay, upper) {
  pts <- upper
  for (i in seq_len(nrow(smp))) {
    if (is_c14(smp$typ[i])) {
      s <- sqrt(smp$wert2[i]^2 + CURVE$err^2)
      idx <- which(abs(smp$wert1[i] - CURVE$c14) < 4 * s)
      if (length(idx) == 0) stop(paste0("Probe ", smp$id[i], ": 14C-Alter liegt ausserhalb der Kalibrationskurve."))
      pts <- c(pts, range(CURVE$year[idx]))
    } else {
      pts <- c(pts, smp$wert1[i], smp$wert2[i])
    }
  }
  pts <- c(pts, lay$tpq[!is.na(lay$tpq)], lay$taq[!is.na(lay$taq)])
  lo <- floor(min(pts)) - 25
  hi <- ceiling(max(pts)) + 25
  if (hi - lo > 5000) stop("Der betrachtete Zeitraum ist zu gross (> 5000 Jahre). Obergrenze oder Daten pruefen.")
  seq(lo, hi)
}

# ---- Tabellen pruefen und umwandeln ------------------------------------------
# Jahre duerfen astronomisch ("-479") oder als Text ("480 v. Chr.") eingegeben werden.
# Bei Dendro-Proben sind wert1/wert2 Jahre, bei 14C-Proben Alter BP und Fehler.
coerce_tables <- function(smp, lay) {
  need_s <- c("id", "schicht", "typ", "wert1", "wert2")
  need_l <- c("schicht", "rang", "tpq", "taq")
  if (is.null(smp) || !all(need_s %in% names(smp))) stop("Proben-Tabelle: Spalten id, schicht, typ, wert1, wert2 erforderlich.")
  if (is.null(lay) || !all(need_l %in% names(lay))) stop("Schicht-Tabelle: Spalten schicht, rang, tpq, taq erforderlich.")
  smp <- smp[!is.na(smp$id), , drop = FALSE]
  lay <- lay[!is.na(lay$schicht), , drop = FALSE]
  smp$schicht <- as.character(smp$schicht); lay$schicht <- as.character(lay$schicht)
  smp$typ <- as.character(smp$typ)
  dend <- !is_c14(smp$typ)
  w1 <- as.character(smp$wert1); w2 <- as.character(smp$wert2)
  smp$wert1 <- ifelse(dend, parse_year(w1), suppressWarnings(as.numeric(w1)))
  smp$wert2 <- ifelse(dend, parse_year(w2), suppressWarnings(as.numeric(w2)))
  lay$rang <- suppressWarnings(as.numeric(as.character(lay$rang)))
  lay$tpq <- parse_year(lay$tpq)
  lay$taq <- parse_year(lay$taq)
  if (nrow(smp) == 0) stop("Proben-Tabelle ist leer.")
  if (anyNA(smp$wert1) || anyNA(smp$wert2)) stop("Proben-Tabelle: wert1/wert2 muessen Zahlen sein (bei Dendro auch '480 v. Chr.' moeglich).")
  if (anyNA(lay$rang)) stop("Schicht-Tabelle: 'rang' muss fuer jede Schicht eine Zahl sein (1 = aelteste/tiefste).")
  if (!all(smp$schicht %in% lay$schicht)) stop("Es gibt Proben mit einer Schicht, die nicht in der Schicht-Tabelle steht.")
  list(smp = smp, lay = lay)
}

# Annahmen ein- und ausschalten (fuer den Vorher/nachher-Vergleich)
apply_assumptions <- function(smp, lay, opt) {
  if (!isTRUE(opt$tpq)) lay$tpq <- NA_real_
  if (!isTRUE(opt$taq)) lay$taq <- NA_real_
  if (!isTRUE(opt$dendro)) smp <- smp[is_c14(smp$typ), , drop = FALSE]
  if (!isTRUE(opt$strat)) lay$rang <- 1
  list(smp = smp, lay = lay)
}
OPT_ALL <- list(strat = TRUE, tpq = TRUE, taq = TRUE, dendro = TRUE)
OPT_NONE <- list(strat = FALSE, tpq = FALSE, taq = FALSE, dendro = FALSE)
OFFMAX <- 500L

# ---- Das Modell -------------------------------------------------------------
# Annahmen (didaktisch):
#  * Jede Schicht hat ein Ablagerungsdatum D.
#  * Jede Probe (14C oder Dendro) datiert ein Ereignis S, das NICHT nach D liegt (S <= D).
#  * Muenz-TPQ: D >= TPQ.  Schriftquelle/Versiegelung TAQ: D <= TAQ.
#  * Stratigrafie: Ablagerung einer tieferen (kleinerer Rang) Schicht <= der hoeheren.
#  * Gleichverteilte Prioren; Berechnung per Monte-Carlo-Ziehen und Verwerfen.
run_model <- function(smp, lay, upper, N = 40000, opt = OPT_ALL, yrs = NULL) {
  ct <- coerce_tables(smp, lay)
  smp <- ct$smp; lay <- ct$lay
  if (is.null(yrs)) yrs <- make_grid(smp, lay, upper)   # Raster aus ALLEN Daten, damit Schalter die Achse nicht verschieben
  ap <- apply_assumptions(smp, lay, opt)
  smp <- ap$smp; lay <- ap$lay
  if (nrow(smp) == 0) stop("Keine Proben uebrig (alle Dendro-Proben ausgeschaltet?).")
  ng <- length(yrs); ns <- nrow(smp); nl <- nrow(lay)

  L <- lapply(seq_len(ns), function(i) {
    if (is_c14(smp$typ[i])) {
      dens_c14(smp$wert1[i], smp$wert2[i], yrs)
    } else {
      d <- as.numeric(yrs >= min(smp$wert1[i], smp$wert2[i]) & yrs <= max(smp$wert1[i], smp$wert2[i]))
      if (sum(d) == 0) stop(paste0("Probe ", smp$id[i], ": Dendro-Fenster liegt nicht im Raster."))
      d / sum(d)
    }
  })
  cdf <- lapply(L, cumsum)
  members <- lapply(seq_len(nl), function(k) which(smp$schicht == lay$schicht[k]))

  # Randverteilung der Ablagerung D je Schicht: P(D=d) ~ Produkt der Verteilungsfunktionen
  Dm <- matrix(0L, N, nl)
  for (k in seq_len(nl)) {
    tpq <- if (is.na(lay$tpq[k])) -Inf else lay$tpq[k]
    taq <- if (is.na(lay$taq[k])) upper else lay$taq[k]
    if (is.finite(tpq) && tpq > taq) {
      stop(paste0("Konflikt in Schicht ", lay$schicht[k], ": Das TPQ (", tpq, ") liegt nach dem TAQ (", taq,
                  "). Ist die Muenze intrusiv, das TAQ falsch oder die Schicht nicht einheitlich? (Jahre astronomisch gezaehlt)"))
    }
    w <- as.numeric(yrs >= tpq & yrs <= taq)
    for (i in members[[k]]) w <- w * cdf[[i]]
    if (sum(w) <= 0) {
      stop(paste0("Konflikt in Schicht ", lay$schicht[k],
                  ": TPQ/TAQ und die Proben schliessen sich aus. Mindestens eine Probe ist juenger als das TAQ oder aelter als erlaubt."))
    }
    Dm[, k] <- sample.int(ng, N, replace = TRUE, prob = w / sum(w))
  }

  # Stratigrafie: Verwerfen aller Ziehungen, die die Reihenfolge verletzen
  ok <- rep(TRUE, N)
  prevmax <- rep(-Inf, N)
  for (r in sort(unique(lay$rang))) {
    ks <- which(lay$rang == r)
    for (k in ks) ok <- ok & (Dm[, k] >= prevmax)
    prevmax <- pmax(prevmax, apply(Dm[, ks, drop = FALSE], 1, max))
  }
  Dacc <- Dm[ok, , drop = FALSE]
  na <- nrow(Dacc)
  if (na < 50) {
    stop("Die stratigrafische Reihenfolge widerspricht den Daten (fast alle Ziehungen verworfen). Rang, TPQ/TAQ und Proben pruefen.")
  }

  postD <- matrix(sapply(seq_len(nl), function(k) smooth5(tabulate(Dacc[, k], ng))), nrow = ng)
  postS <- matrix(0, ng, ns)
  off_med <- rep(NA_real_, ns); off_q95 <- rep(NA_real_, ns); off_dens <- vector("list", ns)   # Abstand zur Ablagerung
  rel_med <- rep(NA_real_, ns); rel_q95 <- rep(NA_real_, ns); rel_dens <- vector("list", ns)   # Abstand zur juengsten Probe der Schicht
  for (k in seq_len(nl)) {
    mem <- members[[k]]
    if (!length(mem)) next
    Di <- Dacc[, k]
    IDX <- matrix(0L, na, length(mem))
    for (j in seq_along(mem)) {
      i <- mem[j]
      u <- runif(na) * cdf[[i]][Di]
      idx <- pmin(findInterval(u, cdf[[i]]) + 1L, Di)
      IDX[, j] <- idx
      postS[, i] <- smooth5(tabulate(idx, ng))
      off <- Di - idx                                   # Jahre zwischen Probenereignis und Ablagerung
      qs <- as.numeric(quantile(off, c(0.5, 0.95)))
      off_med[i] <- qs[1]; off_q95[i] <- qs[2]
      off_dens[[i]] <- smooth_edge(tabulate(pmin(off, OFFMAX) + 1L, OFFMAX + 1L))
    }
    if (length(mem) > 1) {
      ref <- apply(IDX, 1, max)                         # juengstes Probenereignis der Schicht in jeder Ziehung
      for (j in seq_along(mem)) {
        i <- mem[j]
        rel <- ref - IDX[, j]
        qs <- as.numeric(quantile(rel, c(0.5, 0.95)))
        rel_med[i] <- qs[1]; rel_q95[i] <- qs[2]
        rel_dens[[i]] <- smooth_edge(tabulate(pmin(rel, OFFMAX) + 1L, OFFMAX + 1L))
      }
    }
  }
  list(yrs = yrs, smp = smp, lay = lay, L = L, postS = postS, postD = postD,
       members = members, acc = na / N, upper = upper, opt = opt,
       off_med = off_med, off_q95 = off_q95, off_dens = off_dens,
       rel_med = rel_med, rel_q95 = rel_q95, rel_dens = rel_dens)
}

prep_tables <- function(smp, lay) coerce_tables(smp, lay)

# ---- OxCal-Code erzeugen ----------------------------------------------------
# Uebersetzung des App-Modells in eine OxCal-Sequenz:
#  * Rang -> Sequence mit Boundary zwischen den Schichten (gleicher Rang: Phase-Gruppe)
#  * 14C -> R_Date, Dendro-Fenster -> Date mit U(von,bis)
#  * Muenz-TPQ -> Date mit U(tpq, Obergrenze) innerhalb der Phase der Schicht
#  * TAQ -> Date mit U(Untergrenze, taq) direkt NACH der Phase der Schicht
# Hinweis: OxCal modelliert Phasen mit Boundaries, die App die Ablagerung D nach
# den Proben. Die Ergebnisse sind daher aehnlich, aber nicht identisch.
oxcal_code <- function(smp, lay, upper, lo, outlier = FALSE) {
  out <- character()
  add <- function(...) out <<- c(out, paste0(...))
  q <- function(x) gsub("[\"\\\\]", "", as.character(x))
  num <- function(x) format(x, scientific = FALSE, trim = TRUE)

  layer_block <- function(k, ind) {
    sid <- which(smp$schicht == lay$schicht[k])
    nm <- q(lay$schicht[k])
    has_taq <- !is.na(lay$taq[k])
    if (has_taq) {
      add(ind, 'Sequence("', nm, ' mit TAQ")')
      add(ind, "{")
      ind <- paste0(ind, " ")
    }
    add(ind, 'Phase("', nm, '")')
    add(ind, "{")
    for (i in sid) {
      if (is_c14(smp$typ[i])) {
        if (outlier) {
          add(ind, ' R_Date("', q(smp$id[i]), '",', num(smp$wert1[i]), ",", num(smp$wert2[i]), ") { Outlier(0.05); };")
        } else {
          add(ind, ' R_Date("', q(smp$id[i]), '",', num(smp$wert1[i]), ",", num(smp$wert2[i]), ");")
        }
      } else {
        add(ind, ' Date("', q(smp$id[i]), '",U(', num(min(smp$wert1[i], smp$wert2[i])), ",",
            num(max(smp$wert1[i], smp$wert2[i])), "));")
      }
    }
    if (!is.na(lay$tpq[k])) {
      add(ind, ' Date("TPQ Muenze ', nm, '",U(', num(lay$tpq[k]), ",", num(max(upper, lay$tpq[k] + 1)), "));")
    }
    add(ind, "};")
    if (has_taq) {
      add(ind, 'Date("TAQ ', nm, '",U(', num(lo), ",", num(lay$taq[k]), "));")
      ind <- substring(ind, 2)
      add(ind, "};")
    }
  }

  add("// OxCal-Code, erzeugt aus der Shiny-App 'Datierungsteam'")
  add("// Einfuegen unter https://c14.arch.ox.ac.uk/oxcal/OxCal.html (Reiter Input), dann Run.")
  add("// Jahreszahlen astronomisch wie in OxCal (Jahr 0 = 1 v. Chr., -78 = 79 v. Chr.). Kurve: Standard von OxCal (IntCal20).")
  add("Options()")
  add("{")
  add(" Resolution=1;")
  add("};")
  if (outlier) add('Outlier_Model("General",T(5),U(0,4),"t");')
  add("Plot()")
  add("{")
  add(' Sequence("Stratigrafie")')
  add(" {")
  add('  Boundary("Beginn");')
  ranks <- sort(unique(lay$rang))
  for (ri in seq_along(ranks)) {
    ks <- which(lay$rang == ranks[ri])
    if (length(ks) > 1) {
      add('  Phase("Rang ', num(ranks[ri]), '")')
      add("  {")
      for (k in ks) layer_block(k, "   ")
      add("  };")
    } else {
      layer_block(ks, "  ")
    }
    if (ri < length(ranks)) {
      add('  Boundary("', num(ranks[ri]), "/", num(ranks[ri + 1]), '");')
    }
  }
  add('  Boundary("Ende");')
  add(" };")
  add("};")
  out
}

# ---- Wiggle-Matching --------------------------------------------------------
wig_T <- function(start, d) {
  cv <- curve_at(start + d$abstand)
  sum((d$c14 - cv$mu)^2 / (d$fehler^2 + cv$sc^2))
}

# ---- Dendro-Crossdating (simulierte Referenz, Schema) ------------------------
# Die Referenzchronologie ist SIMULIERT (gleiche Referenz fuer alle). Sie zeigt das Prinzip,
# ersetzt aber keine echte Standardchronologie.
with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old, envir = .GlobalEnv)
  })
  set.seed(seed)
  expr
}
dendro_make_ref <- function(n = 300, seed = 2026) {
  with_seed(seed, {
    clim <- as.numeric(arima.sim(list(ar = 0.35), n = n))
    clim <- (clim - mean(clim)) / sd(clim)
    yy <- seq_len(n)
    low <- 0.25 * sin(2 * pi * yy / 170) + 0.15 * sin(2 * pi * yy / 60)
    w <- exp(0.30 * clim + low)
    w / mean(w)
  })
}
dendro_make_sample <- function(ref, n, end, noise, seed) {
  with_seed(seed, {
    seg <- ref[(end - n + 1):end]
    j <- seq_len(n)
    trend <- 0.7 + 1.3 * exp(-j / 35)            # Alterstrend: junge Ringe sind breit
    w <- seg * trend * exp(rnorm(n, 0, noise))  # Standortrauschen
    w / mean(w)
  })
}
# Detrending: Verhaeltnis zum gleitenden Mittel (9 Jahre); an den Raendern mit verkuerztem Fenster,
# damit die Probe ihre volle Laenge behaelt (das letzte Ring ist das Faelljahr!)
dendro_index <- function(x, k = 9) {
  n <- length(x); h <- (k - 1) %/% 2
  ma <- vapply(seq_len(n), function(i) mean(x[max(1, i - h):min(n, i + h)]), numeric(1))
  x / ma
}
# Statistik fuer eine Position e (= Endjahr der Probe in der Referenz)
dendro_stats <- function(refidx, smpidx, e) {
  na4 <- c(r = NA_real_, t = NA_real_, glk = NA_real_, n = NA_real_)
  n <- length(smpidx)
  if (e < n || e > length(refidx)) return(na4)
  a <- refidx[(e - n + 1):e]; b <- smpidx
  ok <- !is.na(a) & !is.na(b)
  m <- sum(ok)
  if (m < 10) return(na4)
  r <- suppressWarnings(cor(a[ok], b[ok]))
  tt <- if (is.finite(r) && abs(r) < 1) r * sqrt(m - 2) / sqrt(1 - r^2) else NA_real_
  da <- diff(a); db <- diff(b)
  ok2 <- !is.na(da) & !is.na(db)
  glk <- if (any(ok2)) 100 * mean(sign(da[ok2]) == sign(db[ok2])) else NA_real_
  c(r = r, t = tt, glk = glk, n = m)
}
DENDRO_REF <- dendro_make_ref()
DENDRO_REF_IDX <- dendro_index(DENDRO_REF)
# Splintholzringe bei Eiche, statistische Spannen (nach Tegel et al. 2022)
DENDRO_SAPWOOD <- list(
  "Großbritannien: 10–55 Splintringe" = c(10, 55),
  "Westdeutschland: 9–33 Splintringe" = c(9, 33),
  "Norddeutschland: 10–30 Splintringe" = c(10, 30),
  "Polen: 9–23 Splintringe" = c(9, 23))
# Faelljahr-Fenster aus dem letzten gemessenen Ring
dendro_window <- function(end, mode, lo, hi, s, margin) {
  if (identical(mode, "wk")) return(c(end - margin, end + margin))
  if (identical(mode, "splint")) return(c(end + max(lo - s, 0), end + max(hi - s, 0)))
  c(end + lo, NA_real_)   # nur Kernholz: nur TPQ
}
nz <- function(x, d) if (is.null(x) || length(x) != 1 || is.na(x)) d else x

# ---- Szenarien (alle fiktiv) -------------------------------------------------
SCEN <- list(
  "Kastell Musterberg (mit Konflikt in Schicht C)" = list(
    smp = "id,schicht,typ,wert1,wert2
A1,A,Dendro,101,103
A2,A,14C,1930,25
A3,A,14C,1900,25
B1,B,14C,1850,25
B2,B,14C,1820,25
B3,B,14C,1790,30
C1,C,14C,1760,25
C2,C,14C,1740,30
C3,C,14C,1900,25",
    lay = "schicht,rang,tpq,taq
A,1,,
B,2,138 n. Chr.,
C,3,271 n. Chr.,260 n. Chr.",
    upper = 400,
    desc = "Holz-Erde-Kastell mit drei Schichten. In Schicht C liegt eine Münze, die frühestens 271 geprägt wurde, die Schriftquelle sagt aber: spätestens 260 aufgegeben. Aufgabe: Konflikt finden, Hypothesen prüfen (Münze intrusiv? TAQ falsch?). Die 14C-Alter sind illustrativ gewählt."),
  "Kastell Musterberg (Konflikt gelöst: Münze intrusiv)" = list(
    smp = "id,schicht,typ,wert1,wert2
A1,A,Dendro,101,103
A2,A,14C,1930,25
A3,A,14C,1900,25
B1,B,14C,1850,25
B2,B,14C,1820,25
B3,B,14C,1790,30
C1,C,14C,1760,25
C2,C,14C,1740,30
C3,C,14C,1900,25",
    lay = "schicht,rang,tpq,taq
A,1,,
B,2,138 n. Chr.,
C,3,,260 n. Chr.",
    upper = 400,
    desc = "Wie oben, aber die Münze in Schicht C wird als intrusiv gewertet (TPQ gestrichen). Aufgabe: Schalten Sie die Annahmen einzeln aus und beobachten Sie, was sich für B und C ändert. Ist C3 Altholz?"),
  "Plateau-Beispiel (Hallstatt-Zeit)" = list(
    smp = "id,schicht,typ,wert1,wert2
A1,A,14C,2520,30
A2,A,14C,2480,30
B1,B,14C,2450,30
B2,B,14C,2480,25
B3,B,Dendro,545 v. Chr.,540 v. Chr.
C1,C,14C,2420,30
C2,C,14C,2600,30",
    lay = "schicht,rang,tpq,taq
A,1,,
B,2,,
C,3,480 v. Chr.,",
    upper = -300,
    desc = "Alle 14C-Alter liegen auf dem Plateau (ca. 800–400 v. Chr.): breite, mehrgipfelige Verteilungen. Aufgabe: Wie stark verengen Stratigrafie, Münze (480 v. Chr.) und das Dendro-Fenster die Daten? Probieren Sie die Schalter einzeln aus."),
  "Altholz-Übung (Brandgrube)" = list(
    smp = "id,schicht,typ,wert1,wert2
G1_Samen,Grube,14C,1850,25
G2_Samen,Grube,14C,1840,25
G3_Balken,Grube,14C,2050,25",
    lay = "schicht,rang,tpq,taq
Grube,1,,300 n. Chr.",
    upper = 400,
    desc = "Eine Brandgrube mit zwei Samen und einem Balkenstück. Der Balken ist deutlich älter. Aufgabe: Im Tab „5 · Altholz-Check“ ansehen, wie viel älter jede Probe ist als die jüngste Probe der Grube. Welche Probe datiert die Grube, welche nur das Holz?")
)
SCEN[["Kollektivgrab (Übung nach M. Hinz 2012)"]] <- list(
  smp = "id,schicht,typ,wert1,wert2
V1,Vornutzung,14C,4680,50
V2,Vornutzung,14C,4485,40
N1,Phase1,14C,4415,29
N2,Phase1,14C,4395,34
N3,Phase2,14C,4355,40
N4,Phase2,14C,4375,31
O1,Opfer,14C,4325,40
O2,Opfer,14C,4340,40
O3,Opfer,14C,4335,40",
  lay = "schicht,rang,tpq,taq
Vornutzung,1,,
Phase1,2,,
Phase2,3,,
Opfer,4,,",
  upper = -2500,
  desc = "Kollektivgrab mit Vornutzung (unter dem Grab), zwei Nutzungsphasen und einem Opfer (drei Getreidekörner in einem Topf) nach der Nutzung. Aufgabe: Gesamtmodell rechnen, Dauer der Nutzung abschätzen. Tab 3, „Kombination“: die drei Körner stammen vom selben Ereignis. Übungsdaten nach M. Hinz (Seminar 2012).")
SCEN_FIRST <- SCEN[[1]]

# ---- UI ---------------------------------------------------------------------
PAL <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")

ui <- fluidPage(
  tags$head(tags$style(HTML("
    textarea { font-family: monospace; font-size: 12px; }
    .hint { color: #555; font-size: 13px; }
    .box { background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 6px; padding: 8px 12px; margin-bottom: 10px; }
  "))),
  titlePanel("Datierungsteam: 14C, Dendro, Münzen und Stratigrafie"),
  fluidRow(column(12, div(style = "text-align:right",
    radioButtons("yrfmt", "Jahre anzeigen als", choices = c("v. Chr. / n. Chr." = "bcad", "astronomisch (−479)" = "astro"), selected = "bcad", inline = TRUE)))),
  tabsetPanel(
    tabPanel("1 · Kalibration & Modell",
      sidebarLayout(
        sidebarPanel(width = 4,
          h4("Szenario"),
          selectInput("scen", NULL, choices = c(names(SCEN), "Eigene Daten (Tabellen unverändert lassen)")),
          actionButton("load_scen", "Szenario laden"),
          uiOutput("scen_desc"),
          h4("Proben"),
          p(class = "hint", "typ = 14C: wert1 = 14C-Alter (BP), wert2 = Fehler. typ = Dendro: wert1/wert2 = Fenster des Fälljahrs."),
          textAreaInput("txt_smp", NULL, SCEN_FIRST$smp, rows = 9, width = "100%"),
          fileInput("f_smp", "oder CSV hochladen", accept = c(".csv", ".txt")),
          h4("Schichten"),
          p(class = "hint", "rang: 1 = älteste/tiefste Schicht. tpq: jüngste Münze. taq: Schriftquelle/Versiegelung. Jahre als „480 v. Chr.“, „101 n. Chr.“ oder astronomisch (−479)."),
          textAreaInput("txt_lay", NULL, SCEN_FIRST$lay, rows = 5, width = "100%"),
          fileInput("f_lay", "oder CSV hochladen", accept = c(".csv", ".txt")),
          div(class = "box",
            strong("Jahreswandler"),
            fluidRow(
              column(5, numericInput("yw_n", NULL, 480, min = 0, step = 1)),
              column(7, radioButtons("yw_era", NULL, choices = c("v. Chr." = "bc", "n. Chr." = "ad"), inline = TRUE))),
            textOutput("yw_out")),
          numericInput("upper", "Späteste denkbare Ablagerung (wenn kein TAQ), astronomisch", SCEN_FIRST$upper, step = 10),
          div(class = "box",
            strong("Annahmen"),
            checkboxInput("a_strat", "Stratigrafische Reihenfolge anwenden", TRUE),
            checkboxInput("a_tpq", "Münz-TPQ anwenden", TRUE),
            checkboxInput("a_taq", "Schriftquelle / TAQ anwenden", TRUE),
            checkboxInput("a_dendro", "Dendro-Fenster verwenden", TRUE),
            checkboxInput("cmp", "Vergleich: Ablagerung nur aus den 14C-Daten (grau)", TRUE)),
          actionButton("run", "Modell rechnen", class = "btn-primary"),
          br(), br(),
          p(class = "hint", "Gestrichelt: reine Kalibration. Farbig gefüllt: Probe im Modell. Schwarz: modelliertes Ablagerungsdatum. Grau: Ablagerung nur aus den 14C-Daten. Rot: TPQ, blau: TAQ.")
        ),
        mainPanel(width = 8,
          textOutput("acc"),
          plotOutput("plot_model", height = "auto"),
          h4("Ergebnis (95,4 %-Bereiche)"),
          tableOutput("tab"),
          downloadButton("dl", "Tabelle als CSV")
        )
      )
    ),
    tabPanel("2 · Wiggle-Matching",
      sidebarLayout(
        sidebarPanel(width = 4,
          h4("Daten simulieren"),
          numericInput("w_true", "Wahres Startjahr (nur für die Simulation), astronomisch", -640, min = -1100, max = -150),
          numericInput("w_n", "Anzahl Proben", 7, min = 3, max = 15),
          numericInput("w_step", "Jahre zwischen den Proben (Ringe)", 10, min = 1, max = 50),
          numericInput("w_err", "Messfehler (14C-Jahre)", 25, min = 10, max = 60),
          actionButton("sim", "Neue Simulation"),
          h4("Daten (oder eigene einfügen)"),
          p(class = "hint", "abstand = Jahre seit der ersten Probe (Ringzählung)."),
          textAreaInput("w_txt", NULL, "", rows = 9, width = "100%"),
          hr(),
          sliderInput("w_start", "Kalenderjahr der ersten Probe (Schieber, astronomisch)", min = -1200, max = -100, value = -800, step = 1, width = "100%"),
          checkboxInput("w_show_T", "Passungskurve (T) für alle Startjahre zeigen", FALSE),
          checkboxInput("w_show_true", "Lösung einblenden", FALSE)
        ),
        mainPanel(width = 8,
          p("Aufgabe: Die Proben stammen aus einer Ringsequenz mit bekanntem Abstand. Schiebt die Sequenz entlang der Kalibrationskurve, bis die Punkte der Kurve folgen."),
          plotOutput("plot_wig", height = "380px"),
          textOutput("wig_text"),
          conditionalPanel("input.w_show_T", plotOutput("plot_T", height = "260px"))
        )
      )
    ),
    tabPanel("3 · SPD & Kombination",
      sidebarLayout(
        sidebarPanel(width = 4,
          radioButtons("sum_mode", NULL, choices = c("Summe: verschiedene Ereignisse" = "sum", "Kombination: gleiches Ereignis" = "cmb")),
          conditionalPanel("input.sum_mode == 'sum'",
            checkboxInput("spd_mod", "Modellierte Verteilungen zusätzlich zeigen", FALSE),
            p(class = "hint", "Es werden die 14C-Proben aus Tab 1 summiert. Eine echte SPD braucht Dutzende bis Hunderte Daten und einen Permutationstest (siehe Code rechts)."),
            p(class = "hint", "Summenkurven zeigen nicht automatisch Bevölkerung oder Aktivität: Probenauswahl, Plateaus und Taphonomie prägen die Form.")),
          conditionalPanel("input.sum_mode == 'cmb'",
            uiOutput("cmb_layer_ui"),
            p(class = "hint", "Kombination gilt nur, wenn alle Messungen dasselbe Ereignis datieren (z. B. mehrere Körner aus einem Topf, mehrere Messungen an einem Skelett). Dann werden die Alter gewichtet gemittelt und erst danach kalibriert. Der χ²-Test prüft, ob die Messungen zusammenpassen."))
        ),
        mainPanel(width = 8,
          conditionalPanel("input.sum_mode == 'sum'",
            plotOutput("plot_spd", height = "340px"),
            h4("Für größere Datensätze in R (rcarbon)"),
            tags$pre("library(rcarbon)
x <- calibrate(x = c14$alter, errors = c14$fehler, calCurves = 'intcal20')
s <- spd(x, timeRange = c(2800, 2200))       # cal BP
plot(s)
# Permutationstest gegen ein Nullmodell:
# modelTest(x, errors = c14$fehler, bins = bins, nsim = 500,
#           timeRange = c(2800, 2200), model = 'exponential', runm = 100)")),
          conditionalPanel("input.sum_mode == 'cmb'",
            plotOutput("plot_cmb", height = "340px"),
            h4(textOutput("cmb_head")),
            textOutput("cmb_chi"),
            textOutput("cmb_hpd"))
        )
      )
    ),
    tabPanel("4 · OxCal-Code",
      br(),
      fluidRow(column(12,
        p("Der Code wird aus den Tabellen in Tab 1 erzeugt (live) und berücksichtigt die dort gewählten Annahmen. In OxCal online im Reiter „Input“ einfügen und „Run“ drücken. Danach die Ergebnisse mit denen der App vergleichen."),
        checkboxInput("ox_out", "Ausreißermodell ergänzen (General, 5 % a priori)", FALSE),
        tags$button("In Zwischenablage kopieren", class = "btn btn-default",
                    onclick = "navigator.clipboard.writeText(document.getElementById('oxcal_code').innerText)"),
        downloadButton("dl_ox", "Als .oxcal speichern"),
        tags$a("OxCal online öffnen", href = "https://c14.arch.ox.ac.uk/oxcal/OxCal.html", target = "_blank", class = "btn btn-link"),
        br(), br(),
        verbatimTextOutput("oxcal_code"),
        p(class = "hint", "Übersetzung: Rang → Sequence mit Boundaries; ¹⁴C → R_Date; Dendro-Fenster → Date mit U(von,bis); Münze → Date U(TPQ,Obergrenze) in der Phase; TAQ → Date U(Untergrenze,TAQ) direkt nach der Phase. OxCal modelliert Phasen mit Boundaries, die App die Ablagerung nach den Proben: Unterschiede in den Ergebnissen sind erwartbar und eine gute Diskussionsfrage.")
      ))
    ),
    tabPanel("5 · Altholz-Check",
      sidebarLayout(
        sidebarPanel(width = 4,
          p("Für jede Probe: Wie viele Jahre liegt das datierte Ereignis (Absterben, Fälljahr) von einem Bezugspunkt entfernt? Große Abstände deuten auf Altholz oder umgelagertes Material hin."),
          radioButtons("off_ref", "Bezugspunkt", choices = c("jüngste Probe derselben Schicht" = "rel", "Ablagerung der Schicht" = "dep"), selected = "rel"),
          sliderInput("off_max", "Anzeige bis (Jahre Abstand)", min = 50, max = 500, value = 300, step = 25),
          sliderInput("off_thr", "Auffällig ab Medianabstand (Jahre)", min = 10, max = 300, value = 60, step = 10),
          p(class = "hint", "„Jüngste Probe der Schicht“ ist robust, braucht aber mindestens zwei Proben in der Schicht. „Ablagerung“ hängt stark vom TAQ bzw. der Obergrenze ab: Ohne feste Obergrenze wird der Abstand für alle Proben groß."),
          p(class = "hint", "Auf Plateaus wird auch bei kurzlebigen Proben der Abstand breit, weil die Kalibration unscharf ist. Ein hoher Median ist ein Hinweis, kein Beweis."),
          p(class = "hint", "Das Modell aus Tab 1 wird verwendet (dort „Modell rechnen“ drücken).")
        ),
        mainPanel(width = 8,
          plotOutput("plot_off", height = "auto"),
          h4("Abstand der Proben"),
          tableOutput("tab_off")
        )
      )
    ),
    tabPanel("6 · Messwert → Alter",
      sidebarLayout(
        sidebarPanel(width = 4,
          radioButtons("f_mode", "Eingabe", choices = c("F14C → 14C-Alter" = "f2a", "14C-Alter → F14C" = "a2f")),
          conditionalPanel("input.f_mode == 'f2a'",
            numericInput("f14", "F14C", 0.737, min = 0.0001, step = 0.001),
            numericInput("f14s", "Messfehler von F14C (1σ)", 0.002, min = 0.00001, step = 0.0005)),
          conditionalPanel("input.f_mode == 'a2f'",
            numericInput("a_in", "14C-Alter (BP)", 2450, step = 10),
            numericInput("a_err", "Fehler (1σ)", 25, min = 1, step = 5)),
          hr(),
          h4("In Tab 1 übernehmen"),
          textInput("f_id", "Proben-ID", "NEU1"),
          textInput("f_layer", "Schicht", "A"),
          actionButton("f_add", "Zur Proben-Tabelle hinzufügen"),
          br(), br(),
          p(class = "hint", "Das 14C-Alter wird mit dem Libby-Wert berechnet: Alter = −8033 · ln(F14C). Die Kalibration gleicht den Unterschied zur wahren Halbwertszeit mit aus.")
        ),
        mainPanel(width = 8,
          h4(textOutput("f14_text")),
          plotOutput("plot_f14", height = "320px"),
          textOutput("f14_hpd"),
          h4("Zum Üben"),
          tableOutput("tab_f14")
        )
      )
    ),
    tabPanel("7 · Dendro-Crossdating",
      sidebarLayout(
        sidebarPanel(width = 4,
          p("Eine Holzprobe wird gegen eine (simulierte) Referenzchronologie verschoben. Wo passt das Ringmuster? Danach: Wie genau ist das Fälljahr?"),
          h4("Probe simulieren"),
          numericInput("d_n", "Anzahl Ringe", 60, min = 20, max = 150, step = 5),
          sliderInput("d_noise", "Rauschen der Probe (Standort)", min = 0.05, max = 0.8, value = 0.25, step = 0.05),
          numericInput("d_end", "Wahres Endjahr (die Lösung), n. Chr.", 102, min = 20, max = 300),
          numericInput("d_seed", "Zufallszahl (gleiche Zahl = gleiche Probe)", 1, min = 1, step = 1),
          actionButton("d_sim", "Probe simulieren"),
          hr(),
          sliderInput("d_pos", "Endjahr der Probe in der Referenz (Schieber)", min = 20, max = 300, value = 180, step = 1, width = "100%"),
          sliderInput("d_thr", "t-Wert-Schwelle (Faustwert, laborabhängig)", min = 2, max = 6, value = 3.5, step = 0.1),
          checkboxInput("d_show_scan", "Gütekurve über alle Positionen zeigen", FALSE),
          checkboxInput("d_show_true", "Lösung einblenden", FALSE),
          hr(),
          h4("Fälljahr bestimmen"),
          radioButtons("d_mode", NULL, choices = c("Waldkante erhalten" = "wk", "Splint, aber keine Waldkante" = "splint", "nur Kernholz" = "kern")),
          selectInput("d_reg", "Region (Splintholz-Spanne)", choices = names(DENDRO_SAPWOOD)),
          numericInput("d_sap", "Gemessene Splintringe", 0, min = 0, max = 80),
          numericInput("d_margin", "Spielraum ± Jahre (bei Waldkante)", 1, min = 0, max = 5),
          hr(),
          h4("In Tab 1 übernehmen"),
          textInput("d_id", "Proben-ID", "A1"),
          textInput("d_layer", "Schicht", "A"),
          actionButton("d_add", "Als Dendro-Probe hinzufügen")
        ),
        mainPanel(width = 8,
          plotOutput("d_plot", height = "330px"),
          textOutput("d_stat"),
          conditionalPanel("input.d_show_scan", plotOutput("d_scan", height = "250px")),
          h4("Beste Positionen"),
          tableOutput("d_tab"),
          h4("Fälljahr"),
          textOutput("d_fell_txt"),
          p(class = "hint", "Die Referenz ist simuliert (Schema). t-Wert nach Baillie-Pilcher, Gleichläufigkeit (Glk) in %. Ringbreiten werden vor dem Vergleich um das gleitende Mittel (9 Jahre) bereinigt. Kurze Proben erzeugen schnell Zufallstreffer: probieren Sie 25 gegen 100 Ringe. Splintholz-Spannen für Eiche nach Tegel et al. 2022.")
        )
      )
    ),
    tabPanel("Hilfe & Grenzen",
      br(),
      h4("Was das Modell macht"),
      tags$ul(
        tags$li("Jede Schicht hat ein Ablagerungsdatum D. Jede Probe datiert ein Ereignis S (Tod des Baums, Tod des Tieres, Fälljahr), das nicht nach D liegt."),
        tags$li("Münzen setzen eine Untergrenze (TPQ) für D, Schriftquellen oder versiegelnde Befunde eine Obergrenze (TAQ)."),
        tags$li("Die Stratigrafie erzwingt D(tiefere Schicht) ≤ D(höhere Schicht)."),
        tags$li("Berechnet wird per Monte-Carlo (Ziehen und Verwerfen), nicht per MCMC. Für den Seminarzweck genügt das, die Ergebnisse sind aber nicht identisch mit OxCal.")
      ),
      h4("Annahmen-Schalter"),
      p("Mit den Schaltern in Tab 1 lässt sich jede Annahme einzeln aus- und einschalten. Die graue Linie zeigt die Ablagerung nur aus den ¹⁴C-Daten, ohne Stratigrafie, TPQ, TAQ und Dendro. So wird sichtbar, was jede Zusatzinformation bringt."),
      h4("Jahreszählung"),
      tags$ul(
        tags$li("Intern astronomisch wie in OxCal (Jahr 0 = 1 v. Chr., −78 = 79 v. Chr.). In den Tabellen sind auch Eingaben wie „480 v. Chr.“ und „101 n. Chr.“ möglich, sie werden umgerechnet."),
        tags$li("Mit der Auswahl oben rechts wechselt die Anzeige zwischen v./n. Chr. und astronomischer Zählung. Felder wie „Späteste denkbare Ablagerung“ erwarten astronomische Jahre; der Jahreswandler hilft.")
      ),
      h4("Dendro-Crossdating"),
      p("Tab 7 arbeitet mit einer simulierten Referenzchronologie. Sie zeigt das Prinzip (Detrending, Korrelation, t-Wert, Gleichläufigkeit), ersetzt aber keine echte Standardchronologie."),
      h4("Nicht enthalten"),
      tags$ul(
        tags$li("Ausreißermodelle (nur im OxCal-Export), Reservoireffekte, Phasen mit Boundary-Prioren, Marine20/SHCal20."),
        tags$li("Fehlerhafte Eingaben werden als Meldung angezeigt. Bei Konflikten (z. B. Münze jünger als TAQ) rechnet das Modell bewusst nicht weiter.")
      ),
      h4("Vorschlag für den Vergleich"),
      p("Dieselben Daten in OxCal (Sequence, Date, Boundary, After/Before) modellieren und die Ergebnisse gegenüberstellen. Wo weichen sie ab, und warum?")
    )
  )
)

# ---- Server -----------------------------------------------------------------
server <- function(input, output, session) {

  get_opt <- function() list(strat = isTRUE(input$a_strat), tpq = isTRUE(input$a_tpq),
                             taq = isTRUE(input$a_taq), dendro = isTRUE(input$a_dendro))

  # ---- Szenarien
  output$scen_desc <- renderUI({
    s <- SCEN[[input$scen]]
    if (is.null(s)) p(class = "hint", "Eigene Daten: Tabellen unten bearbeiten oder CSV hochladen.")
    else p(class = "hint", s$desc)
  })
  observeEvent(input$load_scen, {
    s <- SCEN[[input$scen]]
    if (is.null(s)) return()
    updateTextAreaInput(session, "txt_smp", value = s$smp)
    updateTextAreaInput(session, "txt_lay", value = s$lay)
    updateNumericInput(session, "upper", value = s$upper)
    showNotification("Szenario geladen. Jetzt „Modell rechnen“ drücken.", type = "message", duration = 4)
  })

  # ---- Jahreswandler
  output$yw_out <- renderText({
    n <- input$yw_n
    if (is.null(n) || is.na(n)) return("")
    y <- if (identical(input$yw_era, "bc")) 1 - abs(n) else abs(n)
    paste0(n, if (identical(input$yw_era, "bc")) " v. Chr." else " n. Chr.", "  →  ", y, " (astronomisch)")
  })

  observeEvent(input$f_smp, {
    updateTextAreaInput(session, "txt_smp",
      value = paste(readLines(input$f_smp$datapath, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  })
  observeEvent(input$f_lay, {
    updateTextAreaInput(session, "txt_lay",
      value = paste(readLines(input$f_lay$datapath, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  })

  # ---- Modell
  model <- eventReactive(input$run, {
    tryCatch({
      smp <- read_txt(input$txt_smp)
      lay <- read_txt(input$txt_lay)
      m <- run_model(smp, lay, upper = input$upper, opt = get_opt())
      m$base <- NULL
      if (isTRUE(input$cmp)) {
        m$base <- tryCatch(run_model(smp, lay, upper = input$upper, N = 20000, opt = OPT_NONE, yrs = m$yrs),
                           error = function(e) NULL)
      }
      m
    }, error = function(e) list(error = conditionMessage(e)))
  }, ignoreNULL = FALSE)

  nlayers <- reactive({
    lay <- tryCatch(read_txt(input$txt_lay), error = function(e) NULL)
    if (is.null(lay)) 1 else max(1, nrow(lay))
  })

  output$acc <- renderText({
    m <- model()
    if (!is.null(m$error)) return("")
    on <- c(if (m$opt$strat) "Stratigrafie", if (m$opt$tpq) "TPQ", if (m$opt$taq) "TAQ", if (m$opt$dendro) "Dendro")
    paste0("Aktive Annahmen: ", if (length(on)) paste(on, collapse = ", ") else "keine",
           "   |   Akzeptierte Ziehungen (Stratigrafie erfüllt): ", round(100 * m$acc, 1), " %")
  })

  output$plot_model <- renderPlot({
    m <- model()
    validate(need(is.null(m$error), m$error))
    mode <- input$yrfmt
    nl <- nrow(m$lay)
    op <- par(mfrow = c(nl, 1), mar = c(3.4, 3, 2.2, 1), mgp = c(2, 0.7, 0))
    on.exit(par(op))
    for (k in seq_len(nl)) {
      mem <- m$members[[k]]
      has_base <- !is.null(m$base) && length(m$base$members[[k]]) > 0
      sc <- max(c(m$postD[, k], unlist(m$L[mem]), m$postS[, mem], if (has_base) m$base$postD[, k]))
      plot(NA, xlim = range(m$yrs), ylim = c(0, 1.08), yaxt = "n", xaxt = "n",
           xlab = xlab_years(mode), ylab = "rel. Wahrscheinlichkeit",
           main = paste0("Schicht ", m$lay$schicht[k], " (Rang ", m$lay$rang[k], ")"))
      axis_years(range(m$yrs), mode)
      if (has_base) lines(m$yrs, m$base$postD[, k] / sc, col = "grey55", lwd = 2.5)
      for (j in seq_along(mem)) {
        i <- mem[j]; col <- PAL[(i - 1) %% length(PAL) + 1]
        polygon(c(m$yrs, rev(m$yrs)), c(m$postS[, i] / sc, rep(0, length(m$yrs))),
                col = adjustcolor(col, 0.35), border = col)
        lines(m$yrs, m$L[[i]] / sc, col = col, lty = 2, lwd = 1.5)
      }
      lines(m$yrs, m$postD[, k] / sc, lwd = 3, col = "black")
      if (!is.na(m$lay$tpq[k])) abline(v = m$lay$tpq[k], col = "red", lwd = 2, lty = 3)
      if (!is.na(m$lay$taq[k])) abline(v = m$lay$taq[k], col = "blue", lwd = 2, lty = 3)
      leg <- c(as.character(m$smp$id[mem]), "Ablagerung (Modell)", if (has_base) "Ablagerung nur 14C")
      cols <- c(PAL[(mem - 1) %% length(PAL) + 1], "black", if (has_base) "grey55")
      legend("topright", bty = "n", cex = 0.85, legend = leg, col = cols,
             lwd = c(rep(2, length(mem)), 3, if (has_base) 2.5))
    }
  }, height = function() max(300, 230 * nlayers()))

  results <- reactive({
    m <- model()
    validate(need(is.null(m$error), m$error))
    mode <- input$yrfmt
    s <- data.frame(
      Objekt = m$smp$id, Schicht = m$smp$schicht,
      Typ = ifelse(is_c14(m$smp$typ), "14C", "Dendro"),
      Unmodelliert = vapply(seq_len(nrow(m$smp)), function(i) hpd(m$L[[i]], m$yrs, mode = mode), ""),
      Modelliert = vapply(seq_len(nrow(m$smp)), function(i) hpd(m$postS[, i], m$yrs, mode = mode), ""),
      stringsAsFactors = FALSE)
    s[["Nur 14C (Ablagerung)"]] <- "–"
    d <- data.frame(
      Objekt = paste0("Ablagerung ", m$lay$schicht), Schicht = m$lay$schicht, Typ = "Modell",
      Unmodelliert = "–",
      Modelliert = vapply(seq_len(nrow(m$lay)), function(k) hpd(m$postD[, k], m$yrs, mode = mode), ""),
      stringsAsFactors = FALSE)
    d[["Nur 14C (Ablagerung)"]] <- vapply(seq_len(nrow(m$lay)), function(k) {
      if (!is.null(m$base) && length(m$base$members[[k]]) > 0) hpd(m$base$postD[, k], m$yrs, mode = mode) else "–"
    }, "")
    rbind(s, d)
  })
  output$tab <- renderTable(results(), striped = TRUE, spacing = "s")
  output$dl <- downloadHandler(
    filename = function() "datierung_ergebnis.csv",
    content = function(file) write.csv(results(), file, row.names = FALSE, fileEncoding = "UTF-8"))

  # ---- Altholz-Check
  off_data <- function(m) {
    if (identical(input$off_ref, "dep")) list(dens = m$off_dens, med = m$off_med, q95 = m$off_q95,
                                              xlab = "Jahre zwischen Probenereignis und Ablagerung")
    else list(dens = m$rel_dens, med = m$rel_med, q95 = m$rel_q95,
              xlab = "Jahre älter als die jüngste Probe der Schicht")
  }
  output$plot_off <- renderPlot({
    m <- model()
    validate(need(is.null(m$error), m$error))
    od <- off_data(m)
    nl <- nrow(m$lay)
    xmax <- input$off_max
    op <- par(mfrow = c(nl, 1), mar = c(3.4, 3, 2.2, 1), mgp = c(2, 0.7, 0))
    on.exit(par(op))
    for (k in seq_len(nl)) {
      mem <- m$members[[k]]
      mem <- mem[!vapply(od$dens[mem], is.null, logical(1))]
      # Spitze bei 0 (jüngste Probe) abschneiden, damit die übrigen Kurven sichtbar bleiben
      ymax <- max(c(0.001, unlist(lapply(mem, function(i) od$dens[[i]][6:(xmax + 1)]))))
      plot(NA, xlim = c(0, xmax), ylim = c(0, ymax * 1.15), yaxt = "n",
           xlab = od$xlab, ylab = "Wahrscheinlichkeit",
           main = paste0("Schicht ", m$lay$schicht[k], if (!length(mem)) " (nur eine Probe: kein Vergleich möglich)" else ""))
      if (length(mem)) mtext("Spitze bei 0 abgeschnitten", side = 3, line = 0, adj = 1, cex = 0.7, col = "grey40")
      for (i in mem) {
        col <- PAL[(i - 1) %% length(PAL) + 1]
        lines(0:xmax, od$dens[[i]][1:(xmax + 1)], col = col, lwd = 2.5)
        abline(v = od$med[i], col = col, lty = 3)
      }
      if (length(mem)) legend("topright", bty = "n", cex = 0.85, legend = as.character(m$smp$id[mem]),
                              col = PAL[(mem - 1) %% length(PAL) + 1], lwd = 2.5)
    }
  }, height = function() max(300, 230 * nlayers()))

  output$tab_off <- renderTable({
    m <- model()
    validate(need(is.null(m$error), m$error))
    od <- off_data(m)
    ok <- !is.na(od$med)
    d <- data.frame(Probe = m$smp$id, Schicht = m$smp$schicht,
                    `Median Abstand (Jahre)` = ifelse(ok, as.character(round(od$med)), "–"),
                    `95 %-Quantil (Jahre)` = ifelse(ok, as.character(round(od$q95)), "–"),
                    Hinweis = ifelse(!ok, "kein Vergleich (einzige Probe der Schicht)",
                                     ifelse(od$med > input$off_thr, "auffällig: Altholz oder umgelagert?", "unauffällig")),
                    check.names = FALSE, stringsAsFactors = FALSE)
    d
  }, striped = TRUE, spacing = "s")

  # ---- Wiggle
  true_start <- reactiveVal(NA_real_)
  observeEvent(input$sim, {
    n <- max(3, input$w_n); step <- input$w_step; err <- input$w_err; ys <- input$w_true
    abst <- (0:(n - 1)) * step
    cv <- curve_at(ys + abst)
    age <- round(cv$mu + rnorm(n, 0, err))
    true_start(ys)
    updateTextAreaInput(session, "w_txt",
      value = paste(c("abstand,c14,fehler", paste(abst, age, err, sep = ",")), collapse = "\n"))
  })

  wdat <- reactive({
    d <- tryCatch(read_txt(input$w_txt), error = function(e) NULL)
    validate(need(!is.null(d) && all(c("abstand", "c14", "fehler") %in% names(d)),
                  "Daten mit den Spalten abstand, c14, fehler einfügen oder Simulation starten."))
    d <- d[stats::complete.cases(d[, c("abstand", "c14", "fehler")]), ]
    validate(need(nrow(d) >= 3, "Mindestens 3 Proben nötig."))
    d
  })

  output$plot_wig <- renderPlot({
    d <- wdat(); st <- input$w_start; mode <- input$yrfmt
    yrs <- st + d$abstand
    lo <- st - 60; hi <- max(yrs) + 60
    cv <- CURVE[CURVE$year >= lo & CURVE$year <= hi, ]
    ylim <- range(c(cv$c14 - cv$err, cv$c14 + cv$err, d$c14 - d$fehler, d$c14 + d$fehler))
    plot(NA, xlim = c(lo, hi), ylim = ylim, xaxt = "n", xlab = xlab_years(mode),
         ylab = "14C-Alter (BP)", main = "IntCal20 und verschobene Ringsequenz")
    axis_years(c(lo, hi), mode)
    polygon(c(cv$year, rev(cv$year)), c(cv$c14 - cv$err, rev(cv$c14 + cv$err)),
            col = adjustcolor("grey60", 0.4), border = NA)
    lines(cv$year, cv$c14, lwd = 2, col = "grey30")
    arrows(yrs, d$c14 - d$fehler, yrs, d$c14 + d$fehler, angle = 90, code = 3, length = 0.04, col = "firebrick")
    points(yrs, d$c14, pch = 19, col = "firebrick")
    if (isTRUE(input$w_show_true) && !is.na(true_start())) abline(v = true_start(), lty = 2, col = "darkgreen")
  })

  output$wig_text <- renderText({
    d <- wdat()
    T0 <- wig_T(input$w_start, d)
    crit <- qchisq(0.95, df = nrow(d) - 1)
    paste0("T = ", round(T0, 1), "   (kritischer Wert 5 %: ", round(crit, 1), ", df = ", nrow(d) - 1, ")  → ",
           if (T0 <= crit) "Passung statistisch akzeptabel." else "Passung nicht akzeptabel.")
  })

  output$plot_T <- renderPlot({
    d <- wdat(); mode <- input$yrfmt
    starts <- seq(-1200, -100, by = 1)
    Tv <- vapply(starts, wig_T, numeric(1), d = d)
    crit <- qchisq(0.95, df = nrow(d) - 1)
    plot(starts, Tv, type = "l", lwd = 2, log = "y", xaxt = "n", xlab = xlab_years(mode), ylab = "T (log)",
         main = "Passung der Sequenz in Abhängigkeit vom Startjahr")
    axis_years(range(starts), mode)
    abline(h = crit, col = "red", lty = 3)
    abline(v = input$w_start, col = "firebrick", lwd = 2)
    if (isTRUE(input$w_show_true) && !is.na(true_start())) abline(v = true_start(), lty = 2, col = "darkgreen")
  })

  # ---- OxCal-Export
  oxcode <- reactive({
    p <- tryCatch(prep_tables(read_txt(input$txt_smp), read_txt(input$txt_lay)),
                  error = function(e) list(error = conditionMessage(e)))
    validate(need(is.null(p$error), p$error))
    lo <- tryCatch(min(make_grid(p$smp, p$lay, input$upper)), error = function(e) NA_real_)
    if (is.na(lo)) lo <- min(c(p$lay$tpq, p$lay$taq, input$upper), na.rm = TRUE) - 3000
    ap <- apply_assumptions(p$smp, p$lay, get_opt())
    oxcal_code(ap$smp, ap$lay, input$upper, lo, outlier = isTRUE(input$ox_out))
  })
  output$oxcal_code <- renderText(paste(oxcode(), collapse = "\n"))
  output$dl_ox <- downloadHandler(
    filename = function() "modell.oxcal",
    content = function(file) writeLines(oxcode(), file, useBytes = TRUE))

  # ---- Kombination (gleiches Ereignis)
  output$cmb_layer_ui <- renderUI({
    m <- model()
    if (!is.null(m$error)) return(p(class = "hint", "Bitte zuerst in Tab 1 ein Modell rechnen."))
    n14 <- vapply(seq_len(nrow(m$lay)), function(k) sum(is_c14(m$smp$typ[m$members[[k]]])), numeric(1))
    ch <- m$lay$schicht[n14 >= 2]
    if (!length(ch)) return(p(class = "hint", "Keine Schicht mit mindestens zwei 14C-Proben."))
    selectInput("cmb_layer", "Schicht (alle Proben = dasselbe Ereignis)", choices = ch)
  })
  cmb <- reactive({
    m <- model()
    validate(need(is.null(m$error), m$error))
    validate(need(!is.null(input$cmb_layer), "Keine Schicht mit mindestens zwei 14C-Proben."))
    k <- which(m$lay$schicht == input$cmb_layer)
    idx <- m$members[[k]]
    idx <- idx[is_c14(m$smp$typ[idx])]
    validate(need(length(idx) >= 2, "Mindestens zwei 14C-Proben in der Schicht nötig."))
    ages <- m$smp$wert1[idx]; errs <- m$smp$wert2[idx]
    w <- 1 / errs^2
    mu <- sum(w * ages) / sum(w); se <- sqrt(1 / sum(w))
    chi <- sum(w * (ages - mu)^2); df <- length(idx) - 1
    list(m = m, idx = idx, mu = mu, se = se, chi = chi, df = df, crit = qchisq(0.95, df),
         dens = dens_c14(mu, se, m$yrs))
  })
  output$plot_cmb <- renderPlot({
    r <- cmb(); m <- r$m; mode <- input$yrfmt
    sc <- max(c(r$dens, unlist(m$L[r$idx])))
    plot(NA, xlim = range(m$yrs), ylim = c(0, 1.08), xaxt = "n", yaxt = "n", xlab = xlab_years(mode),
         ylab = "rel. Wahrscheinlichkeit", main = paste0("Kombination in Schicht ", input$cmb_layer))
    axis_years(range(m$yrs), mode)
    for (j in seq_along(r$idx)) {
      i <- r$idx[j]; col <- PAL[(i - 1) %% length(PAL) + 1]
      lines(m$yrs, m$L[[i]] / sc, col = col, lty = 2, lwd = 1.8)
    }
    polygon(c(m$yrs, rev(m$yrs)), c(r$dens / sc, rep(0, length(m$yrs))), col = adjustcolor("black", 0.25), border = "black", lwd = 2.5)
    legend("topright", bty = "n", cex = 0.85, legend = c(as.character(m$smp$id[r$idx]), "Kombination"),
           col = c(PAL[(r$idx - 1) %% length(PAL) + 1], "black"), lwd = c(rep(1.8, length(r$idx)), 2.5), lty = c(rep(2, length(r$idx)), 1))
  })
  output$cmb_head <- renderText({
    r <- cmb(); sprintf("Gewichtetes Mittel: %.0f ± %.0f BP  (aus %d Messungen)", r$mu, r$se, length(r$idx))
  })
  output$cmb_chi <- renderText({
    r <- cmb()
    sprintf("χ² = %.1f bei %d Freiheitsgraden (kritischer Wert 5 %%: %.1f)  →  %s", r$chi, r$df, r$crit,
            if (r$chi <= r$crit) "Messungen passen zusammen: Kombination zulässig." else "Messungen passen NICHT zusammen: nicht dasselbe Ereignis oder Ausreißer?")
  })
  output$cmb_hpd <- renderText({
    r <- cmb(); paste0("95,4 %-Bereich der Kombination: ", hpd(r$dens, r$m$yrs, mode = input$yrfmt))
  })

  # ---- SPD
  output$plot_spd <- renderPlot({
    m <- model()
    validate(need(is.null(m$error), m$error))
    mode <- input$yrfmt
    idx <- which(is_c14(m$smp$typ))
    validate(need(length(idx) > 0, "Keine 14C-Proben vorhanden."))
    spd_raw <- Reduce(`+`, m$L[idx])
    spd_raw <- spd_raw / max(spd_raw)
    plot(m$yrs, spd_raw, type = "l", lwd = 2, xaxt = "n", xlab = xlab_years(mode),
         ylab = "Summe (normiert)", main = paste0("Summenkurve von ", length(idx), " 14C-Daten"), ylim = c(0, 1.05))
    axis_years(range(m$yrs), mode)
    if (isTRUE(input$spd_mod)) {
      spd_mod <- rowSums(m$postS[, idx, drop = FALSE])
      lines(m$yrs, spd_mod / max(spd_mod), lwd = 2, col = "firebrick")
      legend("topright", bty = "n", lwd = 2, col = c("black", "firebrick"), legend = c("unmodelliert", "modelliert"))
    }
  })

  # ---- Messwert -> Alter
  f14_calc <- reactive({
    if (identical(input$f_mode, "a2f")) {
      validate(need(is.finite(input$a_in) && is.finite(input$a_err) && input$a_err > 0, "Alter und Fehler (> 0) eingeben."))
      age <- input$a_in; err <- input$a_err
      f <- exp(-age / 8033); s <- f * err / 8033
    } else {
      validate(need(is.finite(input$f14) && input$f14 > 0, "F14C muss größer als 0 sein."),
               need(is.finite(input$f14s) && input$f14s > 0, "Der F14C-Fehler muss größer als 0 sein."))
      f <- input$f14; s <- input$f14s
      age <- -8033 * log(f); err <- 8033 * s / f
    }
    list(f = f, s = s, age = age, err = err)
  })
  output$f14_text <- renderText({
    r <- f14_calc()
    sprintf("F14C = %.4f ± %.4f   →   14C-Alter = %.0f ± %.0f BP", r$f, r$s, r$age, r$err)
  })
  output$plot_f14 <- renderPlot({
    r <- f14_calc()
    yrs <- tryCatch(age_grid(r$age, r$err), error = function(e) NULL)
    validate(need(!is.null(yrs), "Dieses Alter liegt außerhalb der Kalibrationskurve (oder ist negativ)."))
    mode <- input$yrfmt
    d <- dens_c14(r$age, r$err, yrs)
    plot(yrs, d / max(d), type = "n", xaxt = "n", yaxt = "n", xlab = xlab_years(mode), ylab = "rel. Wahrscheinlichkeit",
         main = "Kalibriertes Datum (IntCal20)", ylim = c(0, 1.05))
    axis_years(range(yrs), mode)
    polygon(c(yrs, rev(yrs)), c(d / max(d), rep(0, length(yrs))), col = adjustcolor("firebrick", 0.35), border = "firebrick")
  })
  output$f14_hpd <- renderText({
    r <- f14_calc()
    yrs <- tryCatch(age_grid(r$age, r$err), error = function(e) NULL)
    if (is.null(yrs)) return("")
    paste0("95,4 %-Bereich: ", hpd(dens_c14(r$age, r$err, yrs), yrs, mode = input$yrfmt))
  })
  output$tab_f14 <- renderTable({
    f <- c(1, 0.9, 0.75, 0.5, 0.25, 0.1, 0.01)
    data.frame(F14C = f, `14C-Alter (BP)` = round(-8033 * log(f)), check.names = FALSE)
  }, digits = 3, striped = TRUE, spacing = "s")
  observeEvent(input$f_add, {
    r <- f14_calc()
    id <- gsub("[, ]", "_", input$f_id); lay <- gsub("[, ]", "_", input$f_layer)
    line <- paste(id, lay, "14C", round(r$age), max(1, round(r$err)), sep = ",")
    cur <- sub("\\s+$", "", input$txt_smp)
    updateTextAreaInput(session, "txt_smp", value = paste0(cur, "\n", line))
    showNotification(paste0("Zeile hinzugefügt: ", line, ". Schicht muss in der Schicht-Tabelle stehen."), type = "message", duration = 5)
  })
  # ---- Dendro-Crossdating
  dsim <- reactiveVal(NULL)
  observeEvent(input$d_sim, {
    L <- length(DENDRO_REF)
    n <- as.integer(max(20, min(150, nz(input$d_n, 60))))
    end_true <- as.integer(min(L, max(n, round(nz(input$d_end, 102)))))
    x <- dendro_make_sample(DENDRO_REF, n, end_true, nz(input$d_noise, 0.25), nz(input$d_seed, 1))
    dsim(list(x = x, n = n, end_true = end_true))
    cand <- setdiff(n:L, (end_true - 25):(end_true + 25))
    if (length(cand)) updateSliderInput(session, "d_pos", value = cand[sample.int(length(cand), 1)])
  })

  dd <- reactive({
    s <- dsim()
    validate(need(!is.null(s), "Bitte „Probe simulieren“ drücken."))
    smpidx <- dendro_index(s$x)
    n <- length(s$x); L <- length(DENDRO_REF)
    ends <- n:L
    M <- t(vapply(ends, function(e) dendro_stats(DENDRO_REF_IDX, smpidx, e), numeric(4)))
    list(ends = ends, M = M, smpidx = smpidx, n = n, end_true = s$end_true)
  })
  d_pos_eff <- reactive({ d <- dd(); min(max(nz(input$d_pos, 100), d$n), length(DENDRO_REF)) })

  output$d_plot <- renderPlot({
    d <- dd(); n <- d$n; L <- length(DENDRO_REF); e <- d_pos_eff()
    z <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
    zr <- z(DENDRO_REF_IDX); zs <- z(d$smpidx)
    xl <- c(max(1, e - n - 15), min(L, e + 15))
    yl <- range(c(zr[xl[1]:xl[2]], zs), na.rm = TRUE)
    plot(NA, xlim = xl, ylim = yl, xlab = "Referenzjahr (n. Chr.)", ylab = "Ringbreitenindex (standardisiert)",
         main = "Referenzchronologie (grau) und verschobene Probe (grün)")
    lines(seq_len(L), zr, col = "grey55", lwd = 2)
    lines((e - n + 1):e, zs, col = "#1E8E5A", lwd = 2.5)
    abline(v = c(e - n + 1, e), lty = 3, col = "#1E8E5A")
    if (isTRUE(input$d_show_true)) abline(v = d$end_true, col = "red", lty = 2, lwd = 2)
    legend("topright", bty = "n", cex = 0.9, legend = c("Referenz", "Probe", if (isTRUE(input$d_show_true)) "wahres Endjahr"),
           col = c("grey55", "#1E8E5A", if (isTRUE(input$d_show_true)) "red"), lwd = c(2, 2.5, if (isTRUE(input$d_show_true)) 2),
           lty = c(1, 1, if (isTRUE(input$d_show_true)) 2))
  })

  output$d_stat <- renderText({
    d <- dd(); e <- d_pos_eff()
    st <- dendro_stats(DENDRO_REF_IDX, dendro_index(dsim()$x), e)
    if (is.na(st["t"])) return("An dieser Position ist kein Vergleich möglich.")
    paste0("Endjahr ", e, " n. Chr.:  r = ", round(st["r"], 2), ",  t = ", round(st["t"], 1),
           ",  Glk = ", round(st["glk"]), " %,  Überlappung ", st["n"], " Ringe  →  ",
           if (st["t"] >= input$d_thr) "über der Schwelle." else "unter der Schwelle: noch nicht überzeugend.")
  })

  output$d_scan <- renderPlot({
    d <- dd(); tv <- d$M[, "t"]; e <- d_pos_eff()
    plot(d$ends, tv, type = "l", lwd = 2, xlab = "Endjahr der Probe in der Referenz", ylab = "t-Wert",
         main = "Passung an allen Positionen")
    abline(h = input$d_thr, col = "red", lty = 3)
    abline(v = e, col = "#1E8E5A", lwd = 2)
    if (any(is.finite(tv))) points(d$ends[which.max(tv)], max(tv, na.rm = TRUE), pch = 19, cex = 1.3)
    if (isTRUE(input$d_show_true)) abline(v = d$end_true, col = "red", lty = 2, lwd = 2)
  })

  output$d_tab <- renderTable({
    d <- dd()
    o <- order(d$M[, "t"], decreasing = TRUE, na.last = NA)
    o <- head(o, 5)
    out <- data.frame(Endjahr = as.character(d$ends[o]), r = sprintf("%.2f", d$M[o, "r"]), `t-Wert` = sprintf("%.1f", d$M[o, "t"]),
                      `Glk (%)` = sprintf("%.0f", d$M[o, "glk"]), Ueberlappung = sprintf("%.0f", d$M[o, "n"]),
                      check.names = FALSE, stringsAsFactors = FALSE)
    names(out)[5] <- "Überlappung"
    out
  }, striped = TRUE, spacing = "s")

  d_fell <- reactive({
    reg <- DENDRO_SAPWOOD[[nz(input$d_reg, names(DENDRO_SAPWOOD)[1])]]
    dendro_window(d_pos_eff(), nz(input$d_mode, "wk"), reg[1], reg[2], nz(input$d_sap, 0), nz(input$d_margin, 0))
  })
  output$d_fell_txt <- renderText({
    w <- d_fell(); mode <- input$yrfmt
    if (is.na(w[2])) paste0("Nur Kernholz: Fälljahr frühestens ", fmt_year(w[1], mode), " (nur TPQ). Als TPQ der Schicht in Tab 1 verwenden.")
    else if (w[1] == w[2]) paste0("Fälljahr: ", fmt_year(w[1], mode))
    else paste0("Fälljahr-Fenster: ", fmt_year(w[1], mode), " bis ", fmt_year(w[2], mode),
                "  (astronomisch ", w[1], " bis ", w[2], ")")
  })
  observeEvent(input$d_add, {
    w <- d_fell()
    if (is.na(w[2])) {
      showNotification("Nur Kernholz liefert nur ein TPQ, kein Fälljahr-Fenster. Bitte als TPQ in der Schicht-Tabelle eintragen.", type = "warning", duration = 7)
      return()
    }
    id <- gsub("[, ]", "_", nz(input$d_id, "D1")); lay <- gsub("[, ]", "_", nz(input$d_layer, "A"))
    line <- paste(id, lay, "Dendro", w[1], w[2], sep = ",")
    cur <- sub("\\s+$", "", input$txt_smp)
    updateTextAreaInput(session, "txt_smp", value = paste0(cur, "\n", line))
    showNotification(paste0("Zeile hinzugefügt: ", line), type = "message", duration = 5)
  })

}

shinyApp(ui, server)
