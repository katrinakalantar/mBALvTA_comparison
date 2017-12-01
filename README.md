# mBALvTA_comparison


**ABSTRACT:**
Pneumonia causes more deaths each year in the United States than any other type of infectious disease.1 The ability to accurately detect etiologic pathogens and distinguish them from insignificant colonizing or commensal microbes is essential for guiding optimal antimicrobial treatments. In patients requiring mechanical ventilation, samples obtained by less invasive tracheal aspirate (TA) have historically been considered inferior to those obtained from mini-bronchial alveolar lavage/telescoping catheter (mBAL) sampling due to suspicion for inevitable contamination from oropharyngeal microbiota. This idea has been challenged, however, by studies demonstrating a lack of clinically significant differences between sample types and a shift towards broader acceptance of TA sampling, as reflected by recent updates in clinical practice guidelines.4 Despite the broad potential implications of this diagnostic sampling approach transition, few studies have employed an unbiased molecular approach to evaluate the microbial differences between mBAL and TA specimens.

This repository contains the code and analysis files for comparison of mBAL and TA sample types by DNA-sequencing.

### mbal_v_ta.Rmd
Primary markdown file containing the functions for analysis.

### /data
**all_mBALvTA.txt** : contains the filenames used in the analysis

**merged_genusrpmphylum.csv** : contains the phylum-level collapsed microbe counts

**phylum_proportions.csv** : contains the normalized phylum proportions for all DNA-seq samples in the analysis

**tavmbal_metadata_noviruses.tsv** : contains the metadata for the files - all viral infections have been removed from the analysis at this point

**/BM_4** : contains the raw microbial counts data per patient - output from DeRisi lab pipeline

### /outputs

all outputs from the markdown script

