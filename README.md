# mBALvTA_comparison


**ABSTRACT:**
Accurate and informative microbiologic testing is essential for guiding diagnosis and management of pneumonia in critically ill patients. Sampling of tracheal aspirate (TA) is less invasive compared to mini-bronchoalveolar lavage (mBAL) and is now recommended as a frontline diagnostic approach in mechanically ventilated patients, despite the historical belief that TA was suboptimal due to contamination from oral microbes. Advancements in metagenomic next generation sequencing (mNGS) now permit assessment of airway microbiota without a need for culture, and as such provide an opportunity to examine differences between mBAL and TA at a resolution previously unachievable. Here, we engaged shotgun mNGS to quantitatively assess the airway microbiome in matched mBAL and TA specimens from a prospective cohort of critically ill adults. We observed moderate differences between sample types across all patients (Pearson correlation of 0.72, 95% CI: 0.68 – 0.76), however we found significant compositional similarity in patients with bacterial pneumonia, whose microbial communities were characterized by a dominant pathogen (Pearson correlation of 0.92, 95% CI: 0.88 – 0.95). In addition, we found that both mBAL and TA were similar in terms of microbial burden, abundance of oropharyngeal taxa, and microbial diversity. Our findings suggest that TA sampling provides a similar assessment of airway microbiota as more invasive testing by mBAL, and that this similarity is most significant in the setting of bacterial pneumonia.

This repository contains the code and analysis files for comparison of mBAL and TA sample types by DNA-sequencing.

### mbal_v_ta.Rmd
Primary markdown file containing the functions for analysis.


### mbal_v_ta.md
Output from mbal_v_ta.Rmd with all analyses configured to show results in the github interface.


### /data/032118/
**all_mBALvTA.txt** : contains the filenames used in the analysis

**background_model.txt** : dummy file containing the background model used when pulling the report files from the pipeline

**tavmbal_metadata_noviruses.tsv** : contains the metadata for the files - all viral infections have been removed from the analysis at this point


**/BM_4** : contains the raw microbial counts data per patient - output from DeRisi lab pipeline

**merged_genusrpm.csv** : contains the genus-level collapsed microbe counts; a single matrix containing all microbe counts for all patients in the study.

**mBAL-XXX.report.csv** : there is a single report.csv file produced per sample run through the DeRisi Lab pipeline. These include the RPM values that were grouped to produce the merged file.
