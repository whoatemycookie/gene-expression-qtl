ls---
title: "Load and explore the data"
teaching: 15
exercises: 30
questions:
- "What data are required for eqtl mapping?"
objectives:
- To provide an example and exploration of data used for eqtl mapping.
keypoints:
- ""
source: Rmd
---



Load the libraries.

~~~
library(tidyverse)
library(knitr)
library(GGally)
library(corrplot)
library(broom)
library(qtl2)
library(qtl2convert)
# the following analysis is from File S1 Attie_eQTL_paper_physiology.Rmd 
# compliments of Daniel Gatti. See Data Dryad entry for more information.
~~~
{: .language-r}

## Physiological Phenotypes

The complete data used in these analyses are available from 
[Data Dryad](https://doi.org/10.5061/dryad.pj105). 

Load in the clinical phenotypes.


~~~
# load the data
load("../data/attie_DO500_clinical.phenotypes.RData")
~~~
{: .language-r}

See the [data dictionary](../data/Attie-232_Attie_DO_Islets-dictionary.csv) to 
see a description of each of these phenotypes.  

#### Phenotype Ranges


~~~
tmp = pheno_clin %>%
  select(num_islets:weight_10wk) %>%
  summarize_all(funs(min, max), na.rm = TRUE) %>%
  gather(phenotype, value) %>%
  mutate(phenotype = str_replace(phenotype, "_min", ".min")) %>%
  mutate(phenotype = str_replace(phenotype, "_max", ".max")) %>%
  separate(phenotype, c("phenotype", "stat"), sep = "\\.") %>%
  mutate(stat = factor(stat, levels = c("min", "max"))) %>%
  spread(key = stat, value = value)
kable(tmp, caption = "Phenotype Ranges")
~~~
{: .language-r}



Table: Phenotype Ranges

|phenotype     |          min|          max|
|:-------------|------------:|------------:|
|food_ave      | 1.882500e+00| 5.536790e+00|
|Glu_10wk      | 5.394220e+01| 5.903783e+02|
|Glu_14wk      | 1.005210e+01| 5.616875e+02|
|Glu_6wk       | 5.049661e+01| 2.112279e+02|
|Glu_tAUC      | 1.585963e+04| 1.000717e+05|
|HOMA_B_0min   | 1.076550e-02| 9.856865e+00|
|HOMA_IR_0min  | 7.430000e-04| 1.244811e+00|
|Ins_10wk      | 3.236660e-02| 1.394590e+01|
|Ins_14wk      | 2.648000e-02| 3.309000e+01|
|Ins_6wk       | 6.055000e-02| 1.504120e+01|
|Ins_per_islet | 8.706667e+00| 2.902618e+02|
|Ins_tAUC      | 4.800000e+00| 1.292177e+03|
|num_islets    | 4.200000e+01| 1.096000e+03|
|TG_10wk       | 5.373719e+01| 5.330684e+02|
|TG_14wk       | 4.923035e+01| 3.508755e+02|
|TG_6wk        | 4.731155e+01| 7.811902e+02|
|weight_10wk   | 1.430000e+01| 5.140000e+01|
|weight_2wk    | 1.250000e+01| 3.350000e+01|
|weight_6wk    | 1.340000e+01| 4.450000e+01|
|WPIC          | 3.658200e+02| 2.314127e+05|

#### Univariate Boxplot


~~~
pheno_clin %>%
  select(num_islets:weight_10wk) %>%
  gather(phenotype, value) %>%
  ggplot(aes(x = phenotype, y = value)) +
    geom_boxplot() +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Distribution of Log Transformed Phenotypes")
~~~
{: .language-r}

<img src="../fig/rmd-02-pheno_boxplot-1.png" title="plot of chunk pheno_boxplot" alt="plot of chunk pheno_boxplot" width="612" style="display: block; margin: auto;" />

Subset the phenotypes to include only those analyzed in the paper.


~~~
# convert sex and DO wave (batch) to factors
pheno_clin$sex = factor(pheno_clin$sex)
pheno_clin$DOwave = factor(pheno_clin$DOwave)
~~~
{: .language-r}

### Figure 1 Boxplots


~~~
pheno_clin %>%
  select(mouse, sex, starts_with("weight")) %>%
  gather(week, value, -mouse, -sex) %>%
  separate(week, c("tmp", "week")) %>%
  mutate(week = factor(week, levels = c("2wk", "6wk", "10wk"))) %>%
  ggplot(aes(week, value, fill = sex)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(title = "Body Weight", y = "Weight")
~~~
{: .language-r}

<img src="../fig/rmd-02-bw_boxplot-1.png" title="plot of chunk bw_boxplot" alt="plot of chunk bw_boxplot" width="612" style="display: block; margin: auto;" />


~~~
pheno_clin %>%
  select(mouse, sex, starts_with("Glu")) %>%
  select(mouse, sex, ends_with("wk")) %>%
  gather(week, value, -mouse, -sex) %>%
  separate(week, c("tmp", "week")) %>%
  mutate(week = factor(week, levels = c("6wk", "10wk", "14wk"))) %>%
  ggplot(aes(week, value, fill = sex)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(title = "Glucose", y = "Glucose")
~~~
{: .language-r}

<img src="../fig/rmd-02-glucose_boxplot-1.png" title="plot of chunk glucose_boxplot" alt="plot of chunk glucose_boxplot" width="612" style="display: block; margin: auto;" />


~~~
pheno_clin %>%
  select(mouse, sex, starts_with("Ins")) %>%
  select(mouse, sex, ends_with("wk")) %>%
  gather(week, value, -mouse, -sex) %>%
  separate(week, c("tmp", "week")) %>%
  mutate(week = factor(week, levels = c("6wk", "10wk", "14wk"))) %>%
  ggplot(aes(week, value, fill = sex)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(title = "Insulin", y = "Insulin")
~~~
{: .language-r}

<img src="../fig/rmd-02-insulin_boxplot-1.png" title="plot of chunk insulin_boxplot" alt="plot of chunk insulin_boxplot" width="612" style="display: block; margin: auto;" />


~~~
pheno_clin %>%
  select(mouse, sex, starts_with("TG")) %>%
  select(mouse, sex, ends_with("wk")) %>%
  gather(week, value, -mouse, -sex) %>%
  separate(week, c("tmp", "week")) %>%
  mutate(week = factor(week, levels = c("6wk", "10wk", "14wk"))) %>%
  ggplot(aes(week, value, fill = sex)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(title = "TG", y = "TG")
~~~
{: .language-r}

<img src="../fig/rmd-02-trig_boxplot-1.png" title="plot of chunk trig_boxplot" alt="plot of chunk trig_boxplot" width="612" style="display: block; margin: auto;" />


~~~
pheno_clin %>%
  select(mouse, sex, num_islets:Ins_tAUC, food_ave) %>%
  gather(phenotype, value, -mouse, -sex) %>%
  ggplot(aes(sex, value, fill = sex)) +
    geom_boxplot() +
    scale_y_log10() +
    facet_wrap(~phenotype, scales = "free_y")
~~~
{: .language-r}

<img src="../fig/rmd-02-fig1_boxplots-1.png" title="plot of chunk fig1_boxplots" alt="plot of chunk fig1_boxplots" width="612" style="display: block; margin: auto;" />

### QA/QC

#### Proportion Missing Data


~~~
tmp = pheno_clin %>% 
  mutate_all(is.na) %>% 
  summarize_all(mean) %>%
  gather(phenotype, value)
kable(tmp, caption = "Proportion of Missing Data")
~~~
{: .language-r}



Table: Proportion of Missing Data

|phenotype          | value|
|:------------------|-----:|
|mouse              | 0.000|
|sex                | 0.000|
|sac_date           | 0.222|
|partial_inflation  | 0.994|
|coat_color         | 0.000|
|oGTT_date          | 0.010|
|FAD_NAD_paired     | 0.814|
|FAD_NAD_filter_set | 0.814|
|crumblers          | 0.962|
|birthdate          | 0.034|
|diet_days          | 0.034|
|num_islets         | 0.034|
|Ins_per_islet      | 0.036|
|WPIC               | 0.036|
|HOMA_IR_0min       | 0.010|
|HOMA_B_0min        | 0.016|
|Glu_tAUC           | 0.010|
|Ins_tAUC           | 0.010|
|Glu_6wk            | 0.000|
|Ins_6wk            | 0.000|
|TG_6wk             | 0.000|
|Glu_10wk           | 0.002|
|Ins_10wk           | 0.002|
|TG_10wk            | 0.002|
|Glu_14wk           | 0.002|
|Ins_14wk           | 0.002|
|TG_14wk            | 0.002|
|food_ave           | 0.000|
|weight_2wk         | 0.000|
|weight_6wk         | 0.002|
|weight_10wk        | 0.002|
|DOwave             | 0.000|

The phenotypes that we're mapping (on the right) are mostly free of missing 
values. The highest are `Ins_per_islet` and WPIC at 3.6%.


Log transform and standardize each phenotype. Consider setting points that are 
more than 5 std. dev. from the mean to NA. Only do this if the final 
distribution doesn't look skewed.


~~~
pheno_clin_log = pheno_clin %>%
                   mutate_if(is.numeric, log)
pheno_clin_std = pheno_clin_log %>%
                   select(mouse, num_islets:weight_10wk) %>%
                   mutate_if(is.numeric, scale)
pheno_clin_std %>%
  select(num_islets:weight_10wk) %>%
  gather(phenotype, value) %>%
  ggplot(aes(x = phenotype, y = value)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = "Distribution of Standardized Phenotypes")
~~~
{: .language-r}

<img src="../fig/rmd-02-pheno_std-1.png" title="plot of chunk pheno_std" alt="plot of chunk pheno_std" width="612" style="display: block; margin: auto;" />


~~~
outliers = pheno_clin_std %>% 
              gather(pheno, value, -mouse) %>%
              filter(abs(value) > 5)
kable(outliers, caption = "Potential Outliers")
~~~
{: .language-r}



Table: Potential Outliers

|mouse |pheno      |     value|
|:-----|:----------|---------:|
|DO093 |num_islets | -6.767430|
|DO372 |num_islets | -6.217230|
|DO093 |WPIC       | -6.161516|
|DO372 |WPIC       | -5.068911|
|DO372 |Glu_10wk   |  5.854390|
|DO213 |Glu_14wk   | -9.031752|
|DO372 |Glu_14wk   |  5.802727|

### All Pairs


~~~
knitr::include_graphics("../fig/rmd-02-pheno_all_pairs-1.png")
~~~
{: .language-r}

<img src="../fig/rmd-02-pheno_all_pairs-1.png" title="plot of chunk pheno_all_pairs_include" alt="plot of chunk pheno_all_pairs_include" width="100%" style="display: block; margin: auto;" />

~~~
ggpairs(select(pheno_clin_log, num_islets:weight_10wk)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
~~~
{: .language-r}

<img src="../fig/rmd-02-pheno_all_pairs-1.png" title="plot of chunk pheno_all_pairs" alt="plot of chunk pheno_all_pairs" width="1080" style="display: block; margin: auto;" />

The HOMA phenotypes have odd distributions.


~~~
saveRDS(pheno_clin_log, 
        file = "../data/pheno_clin_log_outliers_removed.rds")
ggplot(pheno_clin_log, aes(HOMA_IR_0min, HOMA_B_0min, color = DOwave, shape = sex)) +
  geom_point()
~~~
{: .language-r}

<img src="../fig/rmd-02-homair_vs_homab-1.png" title="plot of chunk homair_vs_homab" alt="plot of chunk homair_vs_homab" width="612" style="display: block; margin: auto;" />

There doesn't appear to be a batch effect, but there are a large number of low 
values. Is there some lower bound to the HOMA measurements?


### Tests for sex, wave and diet_days.


~~~
tmp = pheno_clin_log %>%
        select(mouse, sex, DOwave, diet_days, num_islets:weight_10wk) %>%
        gather(phenotype, value, -mouse, -sex, -DOwave, -diet_days) %>%
        group_by(phenotype) %>%
        nest()
mod_fxn = function(df) {
  lm(value ~ sex + DOwave + diet_days, data = df)
}
tmp = tmp %>%
  mutate(model = map(data, mod_fxn)) %>%
  mutate(summ = map(model, tidy)) %>%
  unnest(summ)
# kable(tmp, caption = "Effects of Sex, Wave & Diet Days on Phenotypes")
~~~
{: .language-r}


~~~
tmp %>%
  filter(term != "(Intercept)") %>%
  mutate(neg.log.p = -log10(p.value)) %>%
  ggplot(aes(term, neg.log.p)) +
    geom_point() +
    facet_wrap(~phenotype) +
    labs(title = "Significance of Sex, Wave & Diet Days on Phenotypes") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
rm(tmp)
~~~
{: .language-r}

<img src="../fig/rmd-02-sex_diet_wave_effects-1.png" title="plot of chunk sex_diet_wave_effects" alt="plot of chunk sex_diet_wave_effects" width="612" style="display: block; margin: auto;" />

### Weight vs. Food Intake


~~~
pheno_clin_log %>%
  select(mouse, sex, food_ave:weight_10wk) %>%
  gather(phenotype, value, -mouse, -sex, -food_ave) %>%
  separate(phenotype, c("phenotype", "week")) %>%
  mutate(week = factor(week, levels = c("2wk", "6wk", "10wk"))) %>%
  ggplot(aes(food_ave, value, color = sex)) +
    geom_point() +
    geom_smooth(method = "lm") +
    labs(title = "Food Intake vs. Body Weight", y = "log(Body Weight)") + 
    facet_wrap(~week)
~~~
{: .language-r}



~~~
`geom_smooth()` using formula 'y ~ x'
~~~
{: .output}

<img src="../fig/rmd-02-bw_vs_food-1.png" title="plot of chunk bw_vs_food" alt="plot of chunk bw_vs_food" width="612" style="display: block; margin: auto;" />


~~~
model_fxn = function(df) { lm(value ~ sex*food_ave, data = df) }
tmp = pheno_clin_log %>%
  select(mouse, sex, food_ave:weight_10wk) %>%
  gather(phenotype, value, -mouse, -sex, -food_ave) %>%
  separate(phenotype, c("phenotype", "week")) %>%
  mutate(week = factor(week, levels = c("2wk", "6wk", "10wk"))) %>%
  group_by(week) %>%
  nest() %>%
  mutate(model = map(data, model_fxn)) %>%
  mutate(summ = map(model, tidy)) %>%
  unnest(summ) %>%
  filter(term != "(Intercept)") %>%
  mutate(p.adj = p.adjust(p.value))
# kable(tmp, caption = "Effects of Sex and Food Intake on Body Weight")
~~~
{: .language-r}

### Correlation Plots

Females


~~~
tmp = pheno_clin_log %>% 
        filter(sex == "F") %>%
        select(num_islets:weight_10wk)
tmp = cor(tmp, use = "pairwise")
corrplot.mixed(tmp, upper = "ellipse", lower = "number", 
               main = "Female Clinical Phenotype Correlation")
corrplot.mixed(tmp, upper = "ellipse", lower = "number", 
               main = "Female Clinical Phenotype Correlation")
~~~
{: .language-r}

<img src="../fig/rmd-02-female_corr_plot-1.png" title="plot of chunk female_corr_plot" alt="plot of chunk female_corr_plot" width="1080" style="display: block; margin: auto;" />

Males


~~~
tmp = pheno_clin_log %>% 
        filter(sex == "M") %>%
        select(num_islets:weight_10wk)
tmp = cor(tmp, use = "pairwise")
corrplot.mixed(tmp, upper = "ellipse", lower = "number", 
               main = "Male Clinical Phenotype Correlation")
corrplot.mixed(tmp, upper = "ellipse", lower = "number", 
               main = "Male Clinical Phenotype Correlation")
~~~
{: .language-r}

<img src="../fig/rmd-02-male_corr_plot-1.png" title="plot of chunk male_corr_plot" alt="plot of chunk male_corr_plot" width="1080" style="display: block; margin: auto;" />

~~~
knitr::include_graphics("../fig/rmd-02-female_corr_plot-1.png")
~~~
{: .language-r}

<img src="../fig/rmd-02-female_corr_plot-1.png" title="plot of chunk pheno_female_corrplot" alt="plot of chunk pheno_female_corrplot" width="100%" style="display: block; margin: auto;" />

~~~
knitr::include_graphics("../fig/rmd-02-male_corr_plot-1.png")
~~~
{: .language-r}

<img src="../fig/rmd-02-male_corr_plot-1.png" title="plot of chunk pheno_male_corrplot" alt="plot of chunk pheno_male_corrplot" width="100%" style="display: block; margin: auto;" />
### Founder Allele Frequency

Load in the genoprobs and markers.


~~~
genoprobs = readRDS("../data/genotypes/attie_DO500_genoprobs_v5.rds")
markers = readRDS("../data/marker_grid_0.02cM_plus.rds")
K = calc_kinship(probs = genoprobs, type = "loco", cores = 4)
saveRDS(K, file = "../data/kinship")
map = map_df_to_list(map = markers, pos_column = "pos")
~~~
{: .language-r}


~~~
map_fxn = function(g, m) {
  retval = apply(g, 2:3, mean) %>%
             t() %>%
             data.frame() %>%
             mutate(pos = m) %>%
             gather(founder, prop, -pos)
  return(retval)
}
allele.freq = map2(genoprobs, map, map_fxn)
allele.freq = map2(allele.freq, 1:length(allele.freq), function(af, chr) { mutate(af, chr = chr) })
tmp = allele.freq[[1]]
for(i in 2:length(allele.freq)) {
  tmp = rbind(tmp, allele.freq[[i]])
}
allele.freq = data.frame(tmp)
rm(tmp)

cc = CCcolors
names(cc) = LETTERS[1:8]
ggplot(allele.freq, aes(pos, prop, color = founder)) +
  geom_line() +
  scale_color_manual(values = cc) +
  facet_grid(founder~chr, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.spacing = unit(0.1, "lines")) +
  labs(title = "Attie Founder Allele Proportions")
~~~
{: .language-r}

<img src="../fig/rmd-02-founder_allele_freq-1.png" title="plot of chunk founder_allele_freq" alt="plot of chunk founder_allele_freq" width="864" style="display: block; margin: auto;" />


### QTL Scans


~~~
rownames(pheno_clin_log) = pheno_clin_log$mouse
covar = model.matrix(~sex + DOwave, data = pheno_clin_log)

qtl = scan1(genoprobs = genoprobs, pheno = pheno_clin_log[,12:31, drop = FALSE], kinship = K, addcovar = covar, cores = 2)
~~~
{: .language-r}

### QTL plots


~~~
for(i in 1:ncol(qtl)) {
  plot_scan1(x = qtl, map = map, lodcolumn = i, main = colnames(qtl)[i])
  abline(h = 6, col = 2, lwd = 2)
}
~~~
{: .language-r}

<img src="../fig/rmd-02-qtl_plots-1.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-2.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-3.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-4.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-5.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-6.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-7.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-8.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-9.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-10.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-11.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-12.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-13.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-14.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-15.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-16.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-17.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-18.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-19.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" /><img src="../fig/rmd-02-qtl_plots-20.png" title="plot of chunk qtl_plots" alt="plot of chunk qtl_plots" width="612" style="display: block; margin: auto;" />

### QTL Peaks


~~~
lod_threshold = 6
peaks = find_peaks(scan1_output = qtl, map = map, threshold = lod_threshold, peakdrop = 4, prob = 0.95)
kable(peaks %>% select (-lodindex) %>% arrange(chr, pos), caption = "Phenotype QTL Peaks with LOD >= 6")
~~~
{: .language-r}



Table: Phenotype QTL Peaks with LOD >= 6

|lodcolumn     |chr |        pos|       lod|     ci_lo|      ci_hi|
|:-------------|:---|----------:|---------:|---------:|----------:|
|Ins_14wk      |1   | 135.481336|  6.411892|   3.00000| 136.156743|
|TG_14wk       |2   |  59.584104|  6.941041|  57.06789|  60.941372|
|weight_2wk    |2   | 102.703260|  6.850378|  27.19212| 105.158368|
|TG_10wk       |2   | 107.367424|  6.335945|  98.11262| 109.629183|
|Glu_6wk       |3   |   5.843185|  7.390397|   3.00000|   9.907714|
|WPIC          |4   |  59.732535|  6.010250|  59.16445|  87.661023|
|Ins_per_islet |4   |  59.736655|  6.599126|  59.46325|  63.541039|
|TG_10wk       |5   |  31.962653|  6.358020|  28.82120| 141.946043|
|Ins_14wk      |5   | 151.520204|  6.746210| 149.08909| 151.833620|
|Glu_10wk      |8   |  48.435138|  6.040469|  24.25070| 120.579268|
|Ins_10wk      |10  | 115.899454|  6.049395|  78.43255| 117.979203|
|weight_10wk   |11  |  11.726117|  6.595914|  10.41106|  16.881123|
|weight_6wk    |11  |  11.826027|  6.321523|  10.41106|  16.881123|
|Glu_tAUC      |11  |  71.782692|  6.185468|  61.06009|  94.220647|
|Ins_tAUC      |11  |  83.594665| 11.258841|  83.58553|  84.954440|
|HOMA_B_0min   |11  |  84.954440| 10.246737|  80.76487|  84.954440|
|HOMA_B_0min   |11  |  94.142728|  8.401139|  91.74957|  94.478235|
|HOMA_IR_0min  |12  |  75.147681|  6.636235|  72.84618|  79.844595|
|HOMA_B_0min   |12  |  75.168727|  6.538015|  73.60843|  77.028816|
|Ins_10wk      |12  |  79.450744|  6.310718|  77.15929|  97.135996|
|Ins_14wk      |13  |  38.450528|  6.174213|  34.80443| 109.389360|
|HOMA_IR_0min  |13  |  97.462245|  6.783359|  96.43627|  98.394465|
|TG_10wk       |14  |  55.802312|  6.120587|  46.98388|  86.148708|
|HOMA_IR_0min  |15  |  72.661699|  6.423236|  72.23535|  77.426124|
|Ins_14wk      |17  |  31.648772|  7.570437|  31.34581|  34.961498|
|weight_6wk    |17  |  34.565685|  6.831975|  31.46779|  44.523656|
|weight_10wk   |17  |  34.660915|  6.942955|  31.34289|  44.044398|
|HOMA_IR_0min  |17  |  72.177507|  6.587342|  71.48575|  87.735040|
|HOMA_B_0min   |18  |  47.130148|  6.099928|  15.39198|  47.422415|
|num_islets    |18  |  76.826830|  6.561935|  74.47827|  78.195167|
|food_ave      |19  |  22.638956|  6.266982|  21.84092|  60.056498|
|weight_2wk    |X   |  45.989364|  7.315103|  39.46486|  47.608349|
|weight_6wk    |X   |  45.989364|  6.373361|  39.23834|  49.934828|
|TG_6wk        |X   | 135.388462|  6.441733|  37.92750| 139.946194|



~~~
write_csv(peaks, file = "../data/pheno_clin_QTL_peaks.csv")
~~~
{: .language-r}

### QTL Peaks Figure


~~~
peaks = peaks %>%
          arrange(lodcolumn)
pdf("../fig/phenotype_qtl.pdf", width = 10, height = 8)
plot_peaks(peaks, map, col = c("blue","red"), lwd = 3, tick.height = 0.8, gap = 0, main = "LOD > 6")
box()
dev.off()
~~~
{: .language-r}



~~~
quartz_off_screen 
                2 
~~~
{: .output}



~~~
plot_peaks(peaks, map, col = c("blue","red"), lwd = 3, tick_height = 0.8, gap = 0, main = "LOD > 6")
box()
~~~
{: .language-r}

<img src="../fig/rmd-02-qtl_peaks_figure-1.png" title="plot of chunk qtl_peaks_figure" alt="plot of chunk qtl_peaks_figure" width="612" style="display: block; margin: auto;" />

~~~
#ggplot_peaks(peaks, map, col = c("blue","red"), legend.title = "LOD > 6")
~~~
{: .language-r}


~~~
#pdf("../fig/pheno_qtl_heatmap.pdf", width = 12, height = 8)
#qtl_heatmap(qtl = qtl, map = map, low.thr = 3.5)
#dev.off()
#qtl_heatmap(qtl = qtl, map = map, low.thr = 3.5)
~~~
{: .language-r}


## Gene Expression Phenotypes

### summarized data, matrices w/sample annotation, exp data, gene annotation, 
### like bioconductor summarized experiment etc = whole set of rectangles
