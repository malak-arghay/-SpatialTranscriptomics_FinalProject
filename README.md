# FINAL DEGREE PROJECT : Mapping Y chromosome Activity across the Tumor Landscape in Bladder Cancer patients: a Spatial Transcriptomics Approach 

Welcome to the repository for my final degree project in Spatial Transcriptomics. This study maps mosaic loss of the Y chromosome (LOY) in 33 muscle-invasive bladder cancer samples from the DUTRENEO trial using Visium spatial transcriptomics. We designed a 43-gene ChrY expression signature and computed raw RNA count–based Y_sum per spot. By integrating these data with inferred clonal lineages and consensus molecular subtypes, we show that LOY is heterogeneous, not random, forming clusters at invasive fronts and within specific clones—especially in Basal-Squamous and Stroma-Rich subtypes—and might correlate with chemoresistance. Interestingly, LOY also coincides with tumor regions enriched in immunosuppressive features, suggesting a potential role in immune evasion. These insights establish LOY as a spatially and clonally driven modifier of tumor evolution, immune interaction, and treatment response


## Folders and Contents

(All folders are in the spatial_transcriptomics directory.)

- **Preprocessing/**: Scripts for Seurat object cleaning and metadata annotation.
- **ChrY_Signature/**: Code to compute raw ChrY expression signatures from 43 curated genes.
- **ChrY_Loss_Thresholding/**: Threshold testing and generation of binary LOY maps per spot.
- **Clonal_LOY_Analysis/**: Scripts to map LOY status onto clones.
- **Subtype_And_Response/**: Analysis of LOY by consensus subtype and treatment response.
- **Visualization/**: ggplot2 and Seurat plots used in the final report and figures.

## Contact

For questions, feedback, or collaboration, feel free to contact me at:  
**malak.arghay@alum.esci.upf.edu**
