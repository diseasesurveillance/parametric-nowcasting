# Summarise results

## Convert delay matrix to long data.frame

matdelay_df <- function(X, dates, varname = "cases") {
    X |>
        data.frame() |>
        mutate(date = dates) |>
        tidyr::pivot_longer(
            starts_with("X"), names_to = "delay", values_to = varname,
            names_transform = function(x) as.numeric(sub("X", "", x)) - 1
        )
}

## Obtain predicted mean

predmean_nonparam <- function(samples) {
    q_draws <- samples$draws("q")
    q_means <- apply(as.array(q_draws), 3, mean)
    lambda_draws <- samples$draws("lambda")
    lambda_means <- apply(as.array(lambda_draws), 3, mean)
    matrix(lambda_means) %*% q_means
}

predmean_param <- function(samples, delays, qmodel = "exponential") {
    b_draws <- samples$draws("b")
    b_means <- apply(as.array(b_draws), 3, mean)
    phi_draws <- samples$draws("phi")
    phi_means <- apply(as.array(phi_draws), 3, mean)
    if (qmodel == "exponential") {
        q_means <- 1 - (1 - phi_means) * exp(-matrix(b_means) %*% delays)
    } else if (qmodel == "gompertz") {
        q_means <- exp(log(phi_means) * exp(-matrix(b_means) %*% delays))
    }
    lambda_draws <- samples$draws("lambda")
    lambda_means <- apply(as.array(lambda_draws), 3, mean)
    if (dim(q_means)[1] == 1) {
        return(matrix(lambda_means) %*% q_means)
    } else {
        return(lambda_means * q_means)
    }
}

predmean_df <- function(models, delays, dates, qmodel = "exponential") {
  # obtain predictive means
  modnames <- names(models)
  df <- list()
  for (i in seq_along(models)) {
    if (modnames[i] == "nonparam") {
      df[[modnames[i]]] <- predmean_nonparam(models[[i]])
    } else {
      df[[modnames[i]]] <- predmean_param(models[[i]], delays, qmodel)
    }
  }

  # convert to dataframe
  df <- lapply(df, function(x) matdelay_df(x, dates)) |>
    bind_rows(.id = "method") |>
    mutate(method = factor(method, levels = modnames))
}


## Visualize curves fitting comparison between models

plot_compare_curves <- function(df_train, df_fit,
    mylabs = c("exp" = "(a) Exponential", "gomp" = "(b) Gompertz")) {
    ggplot(df_train) +
        geom_line(aes(delay, q, group = date, alpha = "Empirical"), color = "gray80",
            linewidth = rel(0.3)) +
        geom_point(aes(delay, q, group = date, alpha = "Empirical"), color = "gray30",
            size = rel(0.5)) +
        geom_line(aes(delay, q, alpha = "Theoretical"), color = "red", df_fit,
            linewidth = rel(0.7), linetype = "dashed") +
        labs(title = NULL, y = TeX("$q(d)$"), x = TeX("$d$"), alpha = NULL) +
        facet_wrap(~ model, scales = "free", labeller = as_labeller(mylabs)) +
        theme_classic(9) +
        scale_alpha_manual(values = c(1, 1)) +
        scale_y_continuous(limits = c(0, NA)) +
        theme(
            legend.position = "inside",
            legend.position.inside = c(0.9, 0.15),
            legend.key.width = unit(0.8, "cm"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold")
        )
}

plot_predmean <- function(df_train, df_fit, maxdelay = NULL, model_labels = NULL, model_colors = NULL) {
    # default arguments
    if (is.null(maxdelay)) maxdelay <- max(df_train$delay)
    if (is.null(model_labels)) {
        model_labels <- c(
            "nonparam" = TeX("Non-parametric $q_d$"),
            "exp" = TeX("Parametric $q(d)$"),
            "exp_rw1" = TeX("$q_t(d)$ with random walks"),
            "exp_ou" = TeX("$q_t(d)$ with OU processes"),
            "gom" = TeX("Parametric $q(d)$"),
            "gom_rw1" = TeX("$q_t(d)$ with random walks"),
            "gom_ou" = TeX("$q_t(d)$ with OU processes")
        )
    }
    if (is.null(model_colors)) {
        model_colors <- c(
            "nonparam" = "tomato",
            "exp" = "yellowgreen",
            "exp_rw1" = "#8446c6",
            "exp_ou" = "#4682B4",
            "gom" = "yellowgreen",
            "gom_rw1" = "#8446c6",
            "gom_ou" = "#4682B4"
        )
    }
    model_names <- levels(df_fit$method)

    ggplot(df_train) +
        geom_line(aes(delay, cases), color = "gray30", linewidth = rel(0.3)) +
        geom_line(aes(delay, cases, color = method), data = df_fit, linewidth = rel(0.4)) +
        geom_point(aes(delay, cases), color = "black", size = rel(0.5)) +
        labs(x = TeX("Delay (weeks)"), y = TeX("Cumulative reported cases"),
             color = "Model") +
        facet_wrap(~ date, scales = "free_y", nrow = 4) +
        scale_y_continuous(limits = c(0, NA)) +
        scale_x_continuous(limits = c(0, maxdelay), breaks = seq(0, maxdelay, by = 5)) +
        scale_color_manual(values = model_colors[model_names], labels = model_labels[model_names]) +
        theme_bw(9) +
        theme(
            legend.position = "bottom",
            legend.key.height = unit(0.4, 'cm'),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            panel.grid.major = element_line()
        )
}

## Nowcast summaries

nowcast_summary <- function(samples, date_end, date_by = "1 day", alpha = 0.05) {
    draws <- samples$draws("N")
    date_size <- length(dimnames(draws)$variable)
    dplyr::tibble(
        date = rev(seq(as.Date(date_end), length.out = date_size, by = paste0("-", date_by))),
        mean = apply(draws, 3, mean, na.rm = TRUE),
        # median = apply(draws, 3, quantile, probs = 0.5, na.rm = TRUE),
        lower = apply(draws, 3, quantile, probs = alpha / 2, na.rm = TRUE),
        upper = apply(draws, 3, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
    )
}

nowcast_summaries <- function(models, date_end, date_by = "1 day", alpha = 0.05) {
  modnames <- names(models)
    lapply(models, nowcast_summary, date_end = date_end, date_by = date_by, alpha = alpha) |>
        bind_rows(.id = "model") |>
        mutate(model = factor(model, levels = modnames))
}

nowcast_periods <- function(models_by_date, date_by = "1 day", alpha = 0.05) {
    fun_aux <- function(x) {
        nowcast_summaries(models_by_date[[x]], as.Date(names(models_by_date)[x]), date_by, alpha)
    }
    lapply(seq_along(models_by_date), fun_aux) |>
        setNames(names(models_by_date)) |>
        bind_rows(.id = "date_end") |>
        mutate(date_end = as.Date(date_end))
}


plot_nowcast <- function(df_nowcast, df_summary, by_model = TRUE, model_labels = NULL,
  model_colors = NULL) {
  if (is.null(model_labels)) {
    model_labels <- c(
      "nonparam"  = "bold('Non-parametric'~q[d])",
      "exp"       = "bold('Parametric'~q(d))",
      "exp_rw1"   = "bold(q[t](d)~'with random walks')",
      "exp_ou"    = "bold(q[t](d)~'with OU processes')",
      "gom"       = "bold('Parametric'~q(d))",
      "gom_rw1"   = "bold(q[t](d)~'with random walks')",
      "gom_ou"    = "bold(q[t](d)~'with OU processes')"
    )
  }
  if (is.null(model_colors)) {
    model_colors <- c(
      "nonparam" = "#008080",
      "exp" = "yellowgreen",
      "exp_rw1" = "#8446c6",
      "exp_ou" = "#4682B4",
      "gom" = "yellowgreen",
      "gom_rw1" = "#8446c6",
      "gom_ou" = "#4682B4",
      "Current reported cases" = "gray40",
      "Eventual reported cases" = "tomato"
    )
  }

  model_names <- levels(df_nowcast$model)
  emp_names <- c("Current reported cases", "Eventual reported cases")

  gg <- ggplot(df_nowcast) +
    geom_ribbon(aes(x = date, ymin = lower, ymax = upper, fill = model), alpha = 0.4, show.legend = FALSE) +
    # geom_line(aes(date, median, color = model), linewidth = rel(0.4)) +
    geom_line(aes(date, mean, color = model), linewidth = rel(0.4)) +
    geom_line(aes(date, cases_reported, color = "Current reported cases"), df_summary, linewidth = rel(0.4)) +
    geom_point(aes(date, cases_reported, color = "Current reported cases"), df_summary, size = rel(0.5)) +
    geom_line(aes(date, cases_baseline, color = "Eventual reported cases"), df_summary, linetype = 6)

  if (by_model) {
    gg <- gg +
      facet_grid(date_end ~ model,
        labeller = labeller(model = as_labeller(model_labels, label_parsed)), scales = "free")
  } else {
    gg <- gg +
      facet_wrap(~ date_end, scales = "free_y")
  }

  gg <- gg +
    labs(fill = NULL, color = NULL, x = NULL, y = "Number of cases") +
    theme_bw(9) +
    theme(
      legend.position = c(0.88, 0.94),
      legend.key.height = unit(0.4, 'cm'),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(values = model_colors)

  if (by_model) {
    gg <- gg +
      scale_color_manual(values = model_colors, breaks = emp_names)
  } else {
    model_name <- as.character(unique(df_nowcast$model))
    gg <- gg +
      scale_color_manual(values = model_colors,
        breaks = c(emp_names, model_name),
        labels = c(emp_names, "Predicted cases"),
      )
  }
}

