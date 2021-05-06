Discovering cellular programs of intrinsic and extrinsic drivers of metabolic traits in adipocytes using LipocyteProfile
S. Laber and S. Strobel et al.

LipocyteProfiler is an unbiased high-throughput microscopy profiling assay that uses a combination of cellular stains, making it amenable to large-scale profiling of lipid-accumulating cell types. The following repository stores main analysis and their visualisations of the LipocyteProfiler pipeline using LipocytePainting data and RNA-seq data to survey diverse cellular mechanisms by generating context-, process-, and allele-specific morphological and cellular profiles, that offer insights into metabolic disease risk mechanisms. 


LipocyteProfiles represent morphological and cellular phenotypes and are generated from 3,005 morphological and cellular features that map to three cellular compartments (Cell, Cytoplasm, Nucleus) across four channels differentiating the organelles, namely DNA (Hoechst), Mito (MitoTracker Red which stains mitochondria), AGP (actin, golgi, plasma membrane; stained with Phalloidin (F-actin cytoskeleton) and Wheat Germ Agglutinin ( golgi and plasma membranes)), and Lipid (BODIPY, which stains neutral lipids, multiplexed with SYTO14, which stains nucleoli and cytoplasmic RNA). We showed that LipocyteProfiles can be used to survey diverse cellular mechanisms of cell types, polygenic-risk of metabolic disease and allelic-risk of common complex diseases.
 
LipocyteProfiles of white adipocytes compared to brown adipocytes (generation of LP profiles - different cell types.R; Figure 1g)
LipocyteProfiles of subcutaneous vs. visceral adipocytes across differentiation (generation of UMAP LP profile.R; Figure 2b) 
LipocyteProfiles of primary human hepatocytes with several drug treatments (generation of LP profiles in hepatocytes  with drug treatment.R; Figure; 4f+h) 
LipocyteProfiles of tail ends of polygenic-risk of insulin resistance (generation of LP profiles - high vs low polygenic risk.R; Figure 5b)  
LipocyteProfiles of polygenic-risk for lipodystrophy (generation of LP profiles - linear regression model.R; Figure 6b)
LipocyteProfiles of allelic risk of 18q21.33 and 2p23.3 locus (generation of LP profiles - risk vs non-risk haplotype single loci.R,  Figure 7; generation of LP profiles - risk vs non-risk haplotype.R, link BCL2)

Visualisation tools: 
Heatmap (visualisation of LP profiles -  heatmap.R)
Pie Chart (visualisation of LP profiles - pie chart.R)
Barplot (visualisation of LP profiles -  barplot.R)

Imaging based high-throughput profiling data of samples derived from different donors are pruned to be confounded by both technical and biological variables. We assessed the effect of both extrinsic and intrinsic variance on LipocyteProfiler features by performing:

Variance component analysis across all data collected on LipocyteProfiles of 65 donor-derived differentiating AMSCs. We assessed intrinsic genetic variation compared to the contribution of other possible confounding factors such as batch, adipose depot, T2D status, age, sex, BMI, cell density, month/year of sampling, and passage numbers. (LP variance explained across data set.R; Figure S3)
To obtain a measure of batch-to-batch variance associated with our experimental set-up, we differentiated hWAT, hBAT and SGBS preadipocytes in three independent experiments and tested batch effects by visualizing the data using a Principle component analysis and quantifying it using a Kolmogorov-Smirnov test implemented in the “BEclear” R package. (BEclear.R; PCA hWAT hBAT SGBS.R; Figure S2a)
We performed a k-nearest neighbour (knn) supervised machine learning algorithm implemented in the “class” R package to investigate the accuracy of predicting biological and technical variation. For this analysis the data set, consisting of 3 different cell types (hWAT, hBAT, SGBS) distributed on the 96-well plate, imaged at 4 days of differentiation, was split into equally balanced testing (n=18) and training (n=56) sets. Accuracy of the classification model was predicted based on three different categories cell type, batch and column of the 96-well plate. (knn-analysis of hWAT hBAT SGBS.R; Figure S2)

Networks of LipocyteProfiles and RNA-seq data link gene sets with morphological and cellular features, capturing a broad range of cell activity and identifying relevant cellular processes. We generated those networks using a linear regression model across 2,760 LipocyteProfiler features and expression of 52,170 genes across differentiation in both adipocytes depot. (LMM_gene_LP.R; Figure 3). Using this network we interrogated association of expression of a specific gene to LP features (LP gene profile - features extraction.R Figure 3c) or identified transcriptional pathways of a specific LP feature and correlated genes. (Enrichr pathway analysis.R;  Figure 3b)

Visualisation tools: 
Network using “igraph” R package (visualisation of LMM network.R; Figure 3b) 
Heatmap of LipocyteProfile of genes (visualisation of LMM network.R; Figure 3c)

By linking cellular programs and LipocyteProfiles we generated maps of cellular programs revealing morphological phenotypes and transcriptomic states within adipocytes. Hallmark genes of apoptosis, ROS, OXPHOS and thermogenesis were correlated with LipocyteProfiler data. (LMM_cellularprogramm_LP.R; BCL2 Figure xx)

Visualisation tools: 
Heatmap of genes (code;  BCL2 Figure xx)



