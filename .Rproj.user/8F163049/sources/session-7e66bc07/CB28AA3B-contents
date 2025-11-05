Task 1:

When the clusters are far apart, the gap statistic correctly returns the true
number. As the side length decreases, the algorithm starts merging clusters.
From the figure, we can see that for 2D and 3D, it mergres almost immediately.
For 4D and 5D, it merges at side length = 3. And for 6D, it merges at side
length = 2.

In higher dimensions, the clusters stay more distinct, so the method holds up
longer.

Task 2:

The plot did not stay flat at 4 for all large radii. The method sometimes
reported 5 clusters instead of 4. This is because I kept the same number of 
points per shell while the shell radius increased. The outer shells became
very sparse, so the similarity graph actually broke a shell into more than one
connected component, and spectral clustering treated those as separate clusters.
When the radius got smaller the opposite happened so the algorithm merged
everything and the gap statistic picked 1 cluster. Only around radius = 1 did
the algorithm find 4.

If the d_threshold were to be lowered to 0.8, then the over-merging at small
radii would happen later. If the d_threshold were to be raised to 1.2, I would
probably get 1 big cluster even earlier.