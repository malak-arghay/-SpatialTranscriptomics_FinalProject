# FINAL DEGREE PROJECT : Mapping Y chromosome Activity across the Tumor Landscape in Bladder Cancer patients: a Spatial Transcriptomics Approach 

Welcome to the repository for my final degree project in Spatial Transcriptomics. This study maps mosaic loss of the Y chromosome (LOY) in 33 muscle-invasive bladder cancer samples from the DUTRENEO trial using Visium spatial transcriptomics. We designed a 43-gene ChrY expression signature and computed raw RNA count–based Y_sum per spot. By integrating these data with inferred clonal lineages and consensus molecular subtypes, we show that LOY is heterogeneous, not random, forming clusters at invasive fronts and within specific clones , especially in Basal-Squamous and Stroma-Rich subtypes—and might correlate with chemoresistance. 
Interestingly, LOY also coincides with tumor regions enriched in immunosuppressive features, suggesting a potential role in immune evasion. These insights establish LOY as a spatially and clonally driven modifier of tumor evolution, immune interaction, and treatment response


## Folders and Contents

| Folder name            | Description                                                                 |
|------------------------|-----------------------------------------------------------------------------|
| `Principal analysis`  | Main pipeline scripts to compute Y chromosome scores (Y_signature, Y_sum) and prepare data for Cancer and TME |
| `Molecular subtypes`  | Analyses relating ChrY loss to consensus molecular subtypes                 |
| `Clones`              | Scripts for generating boxplots of ChrY activity grouped by clone           |
| `Treatment response`  | Code to assess the link between ChrY loss and treatment response            |




## Contact

For questions, feedback, or collaboration, feel free to contact me at:  
**malak.arghay@alum.esci.upf.edu**
