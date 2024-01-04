# Co-occurrence-network
Cooccurrence networks of 
- rhizosphere fungal communities and
- phyllosphere fungal communities

were constructed and compared for 

- asymptomatic samples and
- symptomatic samples

with the NetCoMi package (network construction and comparison for microbiome data) (version 1.1.0) in R version 4.2.0 using the taxonomic profiles at generic level (Peschel et al. 2021).

Covariates from field was confounded in the network analysis because of the limitation of each covariate for plant sampling. 

To simplify the networks, with the top 100 most frequent ASVs were used. 

Correlation network analysis was computed with the Sparse Correlations for Compositional (SparCC) method. 

Eigenvector centrality was used for scaling node size and picking hub or keystone taxa. 

Cluster patterns were detected by the 'cluster_fast greedy' algorithm (Clauset et al. 2004) and nodes with the same colors indicated that they were from the same cluster. 

Networks were compared using 1,000 permutations to estimate statistical differences and the P values were adjusted with Benjamini-Hochberg correction (Benjamini and Hochberg 2000). 

Similarities in community structure or clustering between networks were measured using the adjusted Rand index (ARI) and the Jaccard index (Gates et al. 2019).
