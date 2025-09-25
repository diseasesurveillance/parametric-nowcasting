# Data processing

## Convert wide dataframe of reports to long format

reporcases_longer <- function (df) {
    df |>
    tidyr::pivot_longer(matches("^delay[0-9]+"), names_to = "delay", values_to = "cases",
        names_transform = function(x) as.numeric(sub("delay", "", x)))
}

## Obtain summary of reported cases and baseline for a giving ending date

reporcases_summary <- function(df, date_end) {
    date_end <- as.Date(date_end)
    df |>
    filter(date <= date_end) |>
    mutate(available_delay = seq(n()-1, 0)) |>
    reporcases_longer() |>
    filter(delay <= available_delay) |>
    group_by(date) |>
    summarise(cases_reported = max(cases),
              cases_baseline = unique(cases_baseline)
    )
}

reporcases_periods <- function(df, date_end) {
    date_end <- setNames(date_end, date_end)
    lapply(date_end, reporcases_summary, df = df) |>
        bind_rows(.id = "date_end")
}


