# Gene Set Enrichment Analysis
## Shiny Application

Comparing pathway activity states of NCI-60 Genome Scale Metabolic Models with corresponding transcriptomic data using GSEA.

In order to describe a correlation between Angiogenesis and the activity of CMP-N-Acetylneuraminate/Gangliosid pathways in Carciogenesis and Metastation transcriptomic data from the [NCI60- Cellines](https://discover.nci.nih.gov/cellminer/home.do) was used. Pathway activities were determined from the Genome Scale Metabolic Models from [Yizhal et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4238051/). (Integrated transcriptomic data + growth rates, Validation of models using Metabolite Secretion/Uptake).

A GSEA was performed comparing cellines with active/non active CMP-N-Acetylneuraminate activity. As well as a linear combination test, taking the continous (flux) phenotypes into account [Dinu et al.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-212).

No signficant relation between CMP-N-Acetylneuraminate/Gangliosid pathway activity and Angiogenesis was found. There were some indication for a correlation of Gangliosid pathway activity and expression of Epithelial Mesenchymal Transition related genes.
