# Time-varying parametric nowcasting

This repository contains the analysis for the paper:
*Bayesian nowcasting for delay adjustments using time-varying parametric functions of cumulative reporting probability*.

## Repository Structure

- `data/processed/`: Processed datasets used in the case studies.
  - `sari_41002_20090614_20091122.rds`: Weekly severe acute respiratory infections (SARI) in Paraná, Brazil (2009).
  - `huso104_20110512_20110606.rds`: Daily hospitalizations for hemolytic-uremic syndrome (HUS) in Germany (2011).
- `source/`: Source code and functions for the analysis.
  - `source/models/`: Stan implementations of the nowcasting models.
- `scripts/`: R Markdown files for data processing, exploratory analysis, simulation, modeling, and summarization. Rendered outputs for applications are available at: [parametric-nowcasting](https://erickchacon.gitlab.io/parametric-nowcasting/).

## Note

This repository is a reimplementation of the original project:
[YangX-Bit/nowcasting](https://github.com/YangX-Bit/nowcasting) and [ErickChacon/parametric-nowcasting](https://github.com/ErickChacon/parametric-nowcasting).