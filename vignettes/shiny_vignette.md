### Introduction


### About
The initial "About" tab provides an introduction to how Escort works. It serves to evaluate sets of analysis choices, including feature selection, dimension reduction, and trajectory inference methods, along with their hyperparameters.
[![screen shot of the entry page](vignettes/shiny_about.png)](vignettes/shiny_about.png)


### Step1
To initiate the analysis assessing the data's support for trajectory existence, users can click on the "Step 1" button in the left column sidebar. The user inputs cleaned raw and normalized single-cell RNA-seq data in preformatted CSV files and then click “Import” button. Escort checks if the datasets share the same dimensions and showing the number of genes and cells in the dataset; if qualified, it proceeds to assess trajectory feasibility. Two scenarios deemed inappropriate for trajectory fitting include cells from biologically distinct clusters and sets of cells with insufficient heterogeneity. The outcomes are displayed in boxes on the left, while UMAP and t-SNE representations visualize the datasets on the right. If the dataset fails the evaluation, such as detecting distinct clusters, Escort conducts differential expression analysis using the edgeR R package to identify the top 30 differentially expressed genes between pairwise clusters. Additionally, if the dataset is deemed homogeneous, Escort lists highly variable genes and performs Gene Ontology (GO) enrichment analysis for better understanding biological processes.
[![screen shot of step1](vignettes/shiny_step1.png)](vignettes/shiny_step1.png)


### Downstream Analysis
Users can employ Escort for dimension reduction by selecting a technique (e.g., UMAP, MDS, TSNE) and specifying the number of highly variable genes. By default, Escort fits a trajectory using Slingshot. Users can visualize the embedding and trajectory on the left and download the embedding and trajectory information into a .rds file by clicking “Download Trajectory” button.
[![screen shot of dr](vignettes/shiny_dr.png)](vignettes/shiny_dr.png)

### Step2:
In step 2, users upload the .rds files generated in the downstream analysis. Escort evaluates embeddings based on inter-cellular relationships, preservation of similarity relationships, and cell density, presenting results in tables. 
[![screen shot of step2](vignettes/shiny_step2.png)](vignettes/shiny_step2.png)

### Step3
Utilizing the embedding-specific trajectory, Escort estimates the proportion of cells with ambiguous projections along the trajectory, displaying the results in a table.
[![screen shot of step3](vignettes/shiny_step3.png)](vignettes/shiny_step3.png)


### Conclusion
Based on the three-step evaluation, Escort establishes a realistic and standardized benchmark scoring system ranging from negative infinity to two. Higher scores indicate better performance. Embeddings with a score greater than zero are recommended, while those with a score less than or equal to zero are considered non-recommended. Users can easily identify optimal embedding choices for constructing a trajectory through this scoring system. The conclusive results, featuring embeddings and trajectories, are presented in a table sorted by score. Additionally, a visual representation of the top 6 highest-scoring embeddings and trajectories is plotted below the table for enhanced insight.
[![screen shot of conclusion](vignettes/shiny_conclusion.png)](vignettes/shiny_conclusion.png)