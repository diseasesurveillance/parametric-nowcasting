# Exploratory Data Analysis

## Visualization of reporting curves

my_classic <- function(base_size) {
  theme_classic(base_size) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.grid.major = element_line(),
          axis.line.x = element_blank())
}

plot_reported <- function(df, maxdelay = NULL) {
    if (is.null(maxdelay)) maxdelay <- max(df$delay)
    df |>
        ggplot() +
        geom_line(aes(delay, cases, group = date), color = "gray50", linewidth = rel(0.4)) +
        geom_point(aes(delay, cases), color = "gray30", size = rel(0.5)) +
        labs(x = TeX("Delay"), y = TeX("Cumulative reported cases")) +
        facet_wrap(~ date, scales = "free_y", nrow = 4) +
        scale_y_continuous(limits = c(0, NA)) +
        scale_x_continuous(limits = c(0, maxdelay), breaks = seq(0, maxdelay, by = 5)) +
        my_classic(9) +
        annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth = rel(0.8))
}

## Cumulative reporting probabiliy curves

qexp <- function(phi, b) {
    function(delay) 1 - (1 - phi) * exp(-b * delay)
}

qgom <- function(phi, b) {
    function(delay) exp(log(phi) * exp(-b * delay))
}

## Visualize curves fitting comparison

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


