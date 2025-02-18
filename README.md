#   The effects of supervision upon effort during resistance training: A Bayesian analysis of prior data and an experimental study of private strength clinic members

## Abstract
Supervision during resistance training (RT) may enhance strength gains by optimizing trainee effort. We investigated supervision’s role in effort during RT in a unique setting with private strength clinics, where members train either unsupervised (“Core” membership) or supervised by a qualified exercise scientist (“Assisted” membership). Using both retrospective analysis of member training records and a prospective experimental study, we examined supervision’s impact on exercise performance, measured as time under load (TUL), rating of perceived effort (RPE), and rating of perceived discomfort (RPD). Bayesian methods were applied, using empirically informed prior distributions from retrospective data to model the experimental study. The prior sample included 1000 training sessions from each membership type, while the experimental study involved 45 Core members performing both supervised and unsupervised sessions in randomized order, using their current training loads to momentary failure. Our findings suggest that, in real-world settings (in situ), exercise performance differed little between supervised and unsupervised training. However, in our experimental study, supervision significantly improved TUL (Core = 125.12 [95%QI: 113.70, 131.90] sec; Assisted = 147.35 [95%QI: 134.29, 154.81] sec; contrast = -22.10 [95%QI: -26.60, -17.61] sec). RPE was slightly higher with supervision in both real-world (Core = 53% [95%QI: 51%, 55%]; Assisted = 59% [95%QI: 57%, 61%]; contrast = -6% [95%QI: -8%, -4%]) and experimental settings (Core = 81% [95%QI: 75%, 86%]; Assisted = 87% [95%QI: 83%, 91%]; contrast = -6% [95%QI: -10%, -4%]), suggesting trainees push closer to failure under supervision. This was further supported by higher RPD during the experimental study (Core = 6.3 [95%QI: 5.1, 7.3]; Assisted = 7.5 [95%QI: 6.5, 8.3]; contrast = -1.2 [95%QI: -1.6, -0.9]). Overall, these results reinforce prior research on the benefits of supervision in RT, indicating that unsupervised trainees—especially in real-world conditions—likely train with suboptimal effort.

## Reproducibility
This repository contains the necessary files and code to reproduce the analyses, figures, and the manuscript. 

## Usage
To reproduce the analyses, you will need to have R (https://cran.r-project.org/) and RStudio (https://www.rstudio.com/products/rstudio/download/#download) installed on your computer.

To help with reproducibility, this project uses the `renv` R package (see https://rstudio.github.io/renv/articles/renv.html). With `renv`, the state of this R project can be easily loaded as `renv` keeps track of the required R packages (including version), and (if known) the external source from which packages were retrieved (e.g., CRAN, Github). With `renv`, packages are installed to a project specific library rather than your user or system library. The `renv` package must be installed on your machine before being able to benefit from its features. The package can be installed using the following command:

``` r
install.packages("renv")
```

Once you have `renv` installed, you can get a copy of this repository on your machine by clicking the green Code button then choose Download zip. Save to your machine and extract. After extraction, double click the `supervision_tul.Rproj` file in the root directory. This will automatically open RStudio. This will ensure all paths work on your system as the working directory will be set to the location of the `.Rproj` file. Upon opening, RStudio will recognize the `renv` files and you will be informed that the project library is out of sync with the lockfile. At shown in the console pane of RStudio, running `renv::restore()` will install the packages recorded in the lockfile. This could take some time depending on your machine and internet connection.

## Targets analysis pipeline

This project also uses a function based analysis pipeline using
[`targets`](https://books.ropensci.org/targets/). Instead of script based pipelines the `targets` package makes use of functions applied to targets specified within the pipeline. The targets can be viewed in the `_targets.R` file, and any user defined functions are available in `R/functions.r`.

You can view the existing targets pipeline by downloading the `targets_pipeline.html` file and opening it in your browser.

Useful console functions:

- `tar_edit()` opens the make file
- `tar_make()` to run targets
- `tar_visnetwork()` to view pipeline

## Software and packages used

The [`grateful`](https://pakillo.github.io/grateful/index.html) package was used to create citations to all software and packages used in the analysis. The `grateful` report can be viewed by downloading the `grateful-report.pdf` file.

## License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
  [cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

