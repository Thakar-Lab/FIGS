Fuzzy Inference of Gene-Sets (FIGS)


![figs](https://cloud.githubusercontent.com/assets/21067499/19664277/a6f15b0e-9a0d-11e6-95fb-4376af6808ff.png)



FIGS implements optimized fuzzy c-means clustering algorithm to produce functionally related genes (gene-sets) from gene expression data. This package optimizes initial cluster centroids with Ward's hierarchical clustering method that produces robust and stable clustering solution. FIGS GUI interface makes it easy to upload data and set the fuzzy c-means clustering parameters. The package also implements different criteria for associating genes to the clusters (see manual for details). Additionally, the users can upload their list of pathways for enrichment of fuzzy gene-sets. The results (fuzzy gene-sets) are stored in excel files and are also displayed in tabular form. The enrichment results and overlap between the gene-sets are displayed as interactive circular graphs.         

This interactive package is developed by Atif Khan Ph.D., Dejan Katanic and Juilee Thakar Ph.D. from University of Rochester Medical Center, Rochester NY. For futher information regarding FIGS, please contact: 

 
Dr. Juilee Thakar 

601 Elmwood Avenue,

Rochester, NY 14642.

Phone: 5852766925

Email: juilee_thakar@urmc.rochester.edu

https://www.urmc.rochester.edu/labs/thakar-lab/


Other projects of Thakar Lab can be accessed at: http://www.bio-networks.com/ 

 
 ================================
 
 
Files required for running FIGS package:
 
Run_FIGS.m 
FIGS_FCM.m 
HypGeoTest.m
Process_FCM.m
FIGS_INIT_FCM.m 
Display_Circular_Graph.m


Prerequisites for running FIGS: Matlab 
