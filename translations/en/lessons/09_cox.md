# Lesson 09: Analyze Follow-up Time with a Cox Model

## Question for this lesson

Among a group of patients with advanced lung cancer, how are age, sex, and daily activity status (ECOG) at enrollment associated with subsequent survival time? This lesson uses a Cox proportional hazards model to address this question and teaches how to handle censoring, check the proportional hazards assumption, and interpret model results.

This lesson uses real lung cancer follow-up data to teach Cox regression. The genetic scores from 1000 Genomes in the previous lesson come from a different group of people; the two datasets must not be joined row by row.

## Data and sources

The script includes the raw data in both `data/cox/lung_ncctg228.csv` and `.tsv` formats. The data come from the `lung` dataset in the R package `survival` (version 3.8.12 in the installation environment), with 228 people and 10 fields. The dataset documentation states that patients were recruited through the North Central Cancer Treatment Group; `time` is survival time in days; `status=1` means censored and `status=2` means dead; and `sex=1` means male and `sex=2` means female. `ph.ecog` is the ECOG performance status score assessed by a physician; higher values indicate poorer status.

The original study cited in the dataset documentation is: Loprinzi CL, Laurie JA, Wieand HS, et al. *Prospective evaluation of prognostic variables from patient-completed questionnaires*. **Journal of Clinical Oncology**. 1994;12(3):601–607. The dataset credits Terry Therneau as its source. Run `?survival::lung` in R to view the package documentation. The bundled files retain the original fields and missing values; missing values are not guessed or filled in.

## Runtime environment and steps

The course environment includes `survival`, `ggplot2`, and `scales`, which are required for this lesson. Run each code cell from the [cell-based lesson](../notebooks/09_cox.ipynb). To run the complete script, use the `课程资源/course` folder:

```r
source("lessons/09_cox.R")
```

You can also run this command in a terminal from the same folder:

```text
Rscript lessons/09_cox.R
```

The script reads the TSV first, then converts `status == 2` to the event indicator `event`: 1 for an event and 0 for censoring. The primary model requires complete values for `time`, `event`, `age`, `sex`, and `ph.ecog`, so a person missing any model field is excluded. The script writes the missing-value counts for each field and the sample flow to the results directory, making the complete-case definition easy to check.

The script then constructs `Surv(time, event)` and fits this model:

```r
survival::coxph(
  Surv(time, event) ~ age + sex_f + ph.ecog,
  data = analysis,
  ties = "efron"
)
```

Here, `sex_f` explicitly labels levels 1 and 2 as male and female, so the HR for women is relative to men. Age and ECOG are modeled per one-year and one-grade increase, respectively. `ties="efron"` specifies Efron's approximation for events at the same follow-up time. The model table reports regression coefficients, HRs, 95% confidence intervals, and Wald-test P values. The script also plots Kaplan–Meier curves for ECOG 0, 1, and 2+, and runs `cox.zph()` to check the proportional hazards assumption for each term and for the model overall.

The tables and plots are written to `results/09_cox/`: `missing_by_variable.tsv`, `sample_flow.tsv`, `cox_model.tsv`, `ph_test.tsv`, as well as `km_by_ecog.pdf` and `km_by_ecog.png`. The 2+ group in the KM plot combines ECOG grades 2, 3, and 4 for display; in the primary Cox model, ECOG remains a numeric variable representing each one-grade increase. Plot labels are in English so they display reliably across base R graphics devices on different systems.

## How to read the results

The raw data used here contain 228 people, of whom 165 died. The primary model has 227 complete cases, including 164 deaths; one person was excluded because ECOG was missing. The outcome is death during follow-up; censoring must not be treated as “survived until the end of the study.”

- **Age:** HR 1.011 (95% CI 0.993–1.030, P=0.232). After adjusting for sex and ECOG, the estimated HR per additional year is close to 1; this sample does not show a clear statistical association.
- **Women relative to men:** HR 0.575 (95% CI 0.414–0.799, P=0.001). After adjusting for age and ECOG, the estimated instantaneous hazard is lower in women. This observational association alone does not show that sex caused a survival difference.
- **ECOG:** HR 1.590 (95% CI 1.273–1.986, P<0.001). After adjusting for age and sex, the estimated instantaneous hazard is multiplied by about 1.59 for each one-grade increase in ECOG. This interpretation depends on the assumption of a linear effect across grades.

The GLOBAL P value for the proportional hazards test is 0.22; the P values for age, sex, and ECOG are approximately 0.66, 0.13, and 0.15, respectively. The test found no clear evidence of a departure, but a sample-based test cannot prove that the PH assumption holds. Interpret it alongside the curves, research question, and modeling choices.

An HR describes the relative instantaneous event hazard at a given time, conditional on the model. It is not a difference in survival probability or an absolute risk difference. KM curves are unadjusted group descriptions and must not be treated as adjusted model results.

## Follow-up time scale and delayed entry

In this dataset, time starts when a patient enters observation and is measured in days; it is not an age time scale beginning at birth. Age is therefore a baseline predictor here, not the Cox time axis.

If a study uses age as its time axis, participants usually enter observation at different ages. Delayed entry/left truncation may then be needed: each person contributes to the risk set from the age at which they actually enter observation, rather than being treated as if they had been observed from an earlier age. Survival analyses should state “when time starts and when each person enters the risk set.” This lesson uses the recorded follow-up days directly. Switching to age as the time axis would require explaining how age at entry was measured and how left truncation was handled.

## Try it yourself

1. Change the formula to `Surv(time, event) ~ age + sex_f` and rerun. How do the age and sex estimates change when ECOG is no longer adjusted for? Compare the coefficients and intervals, and explain that this reflects a change in the adjustment set.
2. Change the KM groups to ECOG 0, 1, 2, and 3/4. Leave the primary Cox model unchanged for now. Explain why a KM plot can use groups while the Cox model estimates a continuous per-grade effect.
3. After `cox.zph(fit)`, inspect `$table` and report the chi-squared statistic, degrees of freedom, and P value in each row. Try `plot(ph)` and observe whether the residuals show a systematic trend over time.

## Troubleshooting common errors

- **Event count or HR direction looks wrong:** Check that `status == 2` was converted to event 1. In the raw data, 1 means censored and 2 means dead; do not pass the original `status` directly to `Surv()`.
- **The HR for women is interpreted in reverse:** Check `levels(analysis$sex_f)`. This script uses men as the reference and women as the comparison group.
- **The model reports missing values or the sample size differs from the lesson:** Check the input file, the fields included in `complete.cases()`, and the R package version; do not silently impute missing values.
- **A file cannot be found:** Set the RStudio working directory to the teaching project root and confirm that `data/cox/lung_ncctg228.tsv` exists.

## Exercise answer key

1. Without ECOG, the age and sex estimates may change because the set of adjusted covariates changes. The direction and size of any change should be determined by rerunning the model locally. Do not describe the difference as a causal effect of a single variable.
2. Grouped curves display unadjusted survival for patients at different grades; the numeric ECOG term in the primary model assumes the same change in log hazard for each one-grade increase. Grouped curves do not test this linearity assumption.
3. Each `cox.zph` variable row checks the corresponding covariate, while the GLOBAL row checks the model overall. P values test for evidence of a departure from PH; they are not the probability that the model is “correct.” Interpret them alongside residual plots and the analysis design.
