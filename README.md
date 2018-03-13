# mBALvTA_comparison


**ABSTRACT:**
Accurate and informative microbiologic testing is essential for guiding diagnosis and management of pneumonia in critically ill patients. Sampling of tracheal aspirate (TA) is less invasive compared to mini-bronchoalveolar lavage (mBAL) and is now recommended as a frontline diagnostic approach, despite concerns of inferiority due to increased potential for oropharyngeal contamination.  Advancements in metagenomic next generation sequencing (mNGS) now permit assessment of airway microbiota at a resolution previously unachievable using culture-based methods. Here we leverage mNGS to quantitatively assess microbial composition in paired mBAL and TA specimens and find moderate differences that resolve in the setting of bacterial pneumonia. 

This repository contains the code and analysis files for comparison of mBAL and TA sample types by DNA-sequencing.

### mbal_v_ta.Rmd
Primary markdown file containing the functions for analysis.

### /data
**all_mBALvTA.txt** : contains the filenames used in the analysis

**merged_genusrpmphylum.csv** : contains the phylum-level collapsed microbe counts

**phylum_proportions.csv** : contains the normalized phylum proportions for all DNA-seq samples in the analysis

**tavmbal_metadata_noviruses.tsv** : contains the metadata for the files - all viral infections have been removed from the analysis at this point

**/BM_4** : contains the raw microbial counts data per patient - output from DeRisi lab pipeline


