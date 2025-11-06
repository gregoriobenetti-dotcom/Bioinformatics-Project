# Bioinformatics Project — Data.table vs Data.frame Benchmark

This project explores the efficiency and readability of **R data.table** functions compared to traditional **data.frame** approaches for bioinformatic data analysis.  
It was developed starting from 12 Tasks and a final integrated workflow, in the end I propose a Benchmark analysis and a plot regarding the speed of execution of the tasks in **data.frame** and **data.table**.

---

## Project Overview

Each task performs a common bioinformatics or data-manipulation operation — such as filtering, joining, summarizing, or reshaping large datasets — implemented using **data.table** syntax.

A benchmark comparison was made between:
- Base R **data.frame** solutions, and  
- Optimized **data.table** equivalents

The performance differences (execution time and memory usage) demonstrate the speed benefits of `data.table` for large datasets.

---

## Structure

| File | Description |
|------|--------------|
| `Functionsdataframe.R` | All 12 tasks implemented using **data.frame** syntax |
| `Functionsdatatable.R` | All 12 tasks implemented using **data.table** syntax |
| `finaltask_function_df.R` | The final integrative workflow using data.frame |
| `finaltask_function_dt.R` | The final integrative workflow using data.table |
| `Test_report` | R Markdown report comparing data.table performance |
| `Test_report.html` | Rendered HTML report of the RMarkdown file |

---

## Requirements

- **R** ≥ 4.0.0  
- **Packages:**  
  - `data.table`  
  - `microbenchmark`  
  - (optional) `knitr`, `rmarkdown` for generating reports

Install dependencies in R:

```r
install.packages(c("data.table", "microbenchmark", "rmarkdown"))

---

## To reproduce or test the project result:

1) Open the project in RStudio
2) Run this command to load all the data.table functions:
 source("Functionsdatatable.R")
3) At this point you can call any function inside it
4) To generate a ful rMArkDown of the project type:
 rmarkdown::render("Test_report.Rmd")
