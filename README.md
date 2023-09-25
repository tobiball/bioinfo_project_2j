# bioinfo_project_2j

Command generate ROC plot for one score - *example polyphen*

`python3 skeleton_script_create_roc_plot.py -ibench data/HGVS_2020_small_benchmark.tsv -ipred data/vep/HGVS_2020_small_polyphen_scores.tsv -color output/ROCplot_HGVS_2020_small_polyphen.png`

Command to generate ROC plot for all scores

`python3 skeleton_script_create_roc_plot.py -ipred data/vep/HGVS_2020_small_polyphen_scores.tsv -ipred data/vep/HGVS_2020_small_sift_scores.tsv -ipred data/HGVS_2020_small_baseline_scores.tsv -ibench data/HGVS_2020_benchmark.tsv -o output/ROCplot_all.png
`
