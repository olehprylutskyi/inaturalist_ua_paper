# State-of-the-art of iNaturalist as a source of data on Ukrainian Funga
Repository contains data and script to reproduce all analysis and visualizations from the manuscript "State-of-the-art of iNaturalist as a source of data on Ukrainian Funga", submitted to the Plant Introduction journal.

![iNaturalist fungal observations as of January 1st, 2024]()

# Data
We downloaded an archive of all observation, available on iNaturalist as of December 31, 2023, and meet following criteria: (i) lie within official state boundary of Ukraine, and (ii) belong to the Fungi Kingdom (iconic taxon name = “Fungi”, https://www.inaturalist.org/taxa/47170-Fungi). No additional filtering was applied. Total number of downloaded records was 62,255.
State boundary of Ukraine, as well as its administrative division, acquired using rgeoboundaries package v. 1.2.9 (Dicko, 2023)

# Code
To reproduce analyses, download the repository into a local computer and run `script_inat_ua.R`.

# Requirements