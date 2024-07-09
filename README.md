# State-of-the-art of iNaturalist as a source of data on Ukrainian Funga
Repository contains data and script to reproduce all analysis and visualizations from the manuscript "State-of-the-art of iNaturalist as a source of data on Ukrainian Funga".

This study evaluates the utility of iNaturalist, a citizen science platform, for mycological research in Ukraine. With over 69,000 fungal observations documented since its adoption, iNaturalist presents a significant potential for supplementing professional mycological surveys. This work is novel in critically analysing the spatial distribution, taxonomic coverage, and identification accuracy of iNaturalist fungal data, with a particular focus on social processes underlying citizen scientists' activity in Ukraine over the last five years. We explored observers' activity and compared iNaturalist records with curated checklists and additional data sources. Despite inherent biases, our analysis shows that iNaturalist effectively documents conspicuous macroscopic fungi, often surpassing traditional sources in record numbers. However, limitations such as identification uncertainty and lack of voucher specimens persist, necessitating greater professional engagement. Enhanced collaboration between amateur and professional mycologists could unlock iNaturalist’s full potential, making it a robust tool for biodiversity monitoring and research in Ukraine.

![Taxonomic credibility of iNaturalist fungal data from Ukraine as of June 1st, 2024](https://github.com/olehprylutskyi/inaturalist_ua_paper/blob/main/figures/fig1_taxonomy.png)

![Spatial distribution of iNaturalist fungal data from Ukraine](https://github.com/olehprylutskyi/inaturalist_ua_paper/blob/main/figures/fig2_distribution.png)

![The social facet of fungal part of iNaturalist activity in Ukraine](https://github.com/olehprylutskyi/inaturalist_ua_paper/blob/main/figures/fig3_social.png)

# Data
All data are stored into the `data` subfolder.

We downloaded an archive of all observations, available on iNaturalist as of 1 June 2024, and meet the following criteria: (i) lie within the official state boundary of Ukraine, (ii) belong to the Fungi Kingdom, and (iii) verifiable (accompanied by photos). Total number of downloaded records was 69,211.

The state boundary of Ukraine, as well as its administrative division, was acquired using the rgeoboundaries package v. 1.2.9 (Dicko 2023). Population density data for Ukraine were acquired from the Gridded Population of the World, Version 4 (GPWv4): Population Density, Revision 11 (https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11). Ukrainian war events data were downloaded from ACLED's Ukraine Conflict Monitor (https://acleddata.com/ukraine-conflict-monitor/).

# Code
To reproduce analyses, download the repository into a local computer and run `inat_ua.Rmd` either as a single script or chunk by chunk. Knitted HTML version inat_ua.html can be viewed with any web-browser.

# Requirements
R version 4.3 or above, packages `sf`, `tidyverse`, `rgeoboundaries`, `iNEXT`, `knitr`, `rgbif`, `terra`, `ggspatial`, `gridExtra`.