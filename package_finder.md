---
title: "package finder"
author: "Daniel Fuller"
date: "04/03/2020"
output:
      html_document:
        keep_md: true
---




```r
library(tidyverse)
```

```
## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──
```

```
## ✔ ggplot2 3.2.1     ✔ purrr   0.3.3
## ✔ tibble  2.1.3     ✔ dplyr   0.8.3
## ✔ tidyr   1.0.0     ✔ stringr 1.4.0
## ✔ readr   1.3.1     ✔ forcats 0.4.0
```

```
## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
library(stringr)
library(dplyr)
library(kableExtra)
```

```
## 
## Attaching package: 'kableExtra'
```

```
## The following object is masked from 'package:dplyr':
## 
##     group_rows
```

```r
library(cranly)
```
## Creating a task view for physical activity researchers

First step is use the `cranly` package to see what is out there. I found this tutorial to help https://rviews.rstudio.com/2018/05/31/exploring-r-packages/. Now I have a tidy dataframe with all the packages. 


```r
package_db <- clean_CRAN_db(tools::CRAN_package_db())
```

## Search dataframe for census and canada

I have used the search terms:

- can
- canada
- census


```r
search <-
    c("canada",
      "census", 
      "cchs"
    )

package_db$search <-
    str_extract_all(package_db$description, paste(search, collapse = "|"))
canada_data_packages <- filter(package_db, search != "character(0)")
```

From this search there are 10 packages are designed to analyse some type of physical activity data.




```r
package_names <- canada_data_packages$package
```

## Packages' description


```r
canada_data_packages %>% select(package, description) %>%
   kable() %>%
   kable_styling("striped", full_width = F) %>% 
   column_spec(2, width_max = 100)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> package </th>
   <th style="text-align:left;"> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> acs </td>
   <td style="text-align:left;max-width: 100; "> Provides a general toolkit for downloading, managing,
  analyzing, and presenting data from the U.S. Census
  (&lt;https://www.census.gov/data/developers/data-sets.html&gt;), including
  SF1 (Decennial short-form), SF3 (Decennial long-form), and the
  American Community Survey (ACS).  Confidence intervals provided with
  ACS data are converted to standard errors to be bundled with
  estimates in complex acs objects.  Package provides new methods to
  conduct standard operations on acs objects and present/plot data in
  statistically appropriate ways. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AmostraBrasil </td>
   <td style="text-align:left;max-width: 100; "> Generates samples or complete list of Brazilian IBGE (Instituto Brasileiro de Geografia e Estatistica, see
  &lt;http://www.ibge.gov.br/&gt; for more information) census
    households, geocoding it by Google Maps. The package connects IBGE site and
    downloads maps and census data. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> apyramid </td>
   <td style="text-align:left;max-width: 100; "> Provides a quick method for visualizing non-aggregated line-list
    or aggregated census data stratified by age and one or two categorical
    variables (e.g. gender and health status) with any number of values. It
    returns a 'ggplot' object, allowing the user to further customize the
    output. This package is part of the 'R4Epis' project 
    &lt;https://r4epis.netlify.com&gt;. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cancensus </td>
   <td style="text-align:left;max-width: 100; "> Integrated, convenient, and uniform access to Canadian
    Census data and geography retrieved using the 'CensusMapper' API. This package produces analysis-ready 
    tidy data frames and spatial data in multiple formats, as well as convenience functions
    for working with Census variables, variable hierarchies, and region selection. API
    keys are freely available with free registration at &lt;https://censusmapper.ca/api&gt;.
    Census data and boundary geometries are reproduced and distributed on an "as
    is" basis with the permission of Statistics Canada (Statistics Canada 2001; 2006;
    2011; 2016). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cchs </td>
   <td style="text-align:left;max-width: 100; "> Contains a function, also called 'cchs', that calculates Estimator III of Borgan et al (2000), &lt;DOI:10.1023/A:1009661900674&gt;. This estimator is for fitting a Cox proportional hazards model to data from a case-cohort study where the subcohort was selected by stratified simple random sampling. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Census2016 </td>
   <td style="text-align:left;max-width: 100; "> Contains selected variables from the time series profiles for statistical areas level 2 from the 2006, 2011, and 2016 censuses of population and housing, Australia. Also provides methods for viewing the questions asked for convenience during analysis. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> censusxy </td>
   <td style="text-align:left;max-width: 100; "> Provides access to the U.S. Census Bureau's API for batch geocoding 
    American street addresses (&lt;https://geocoding.geo.census.gov/geocoder&gt;).
    The package offers a batch solution for address geocoding, as opposed to geocoding
    a single address at a time. It has also been developed specifically with large 
    data sets in mind - only unique addresses are passed to the API for geocoding. 
    If a data set exceeds 1,000 unique addresses, it will be automatically subset 
    into appropriately sized API calls, geocoded, and then put back together so that 
    a single object is returned. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DDM </td>
   <td style="text-align:left;max-width: 100; "> A set of three two-census methods to the estimate the degree of death registration coverage for a population. Implemented methods include the Generalized Growth Balance method (GGB), the Synthetic Extinct Generation method (SEG), and a hybrid of the two, GGB-SEG. Each method offers automatic estimation, but users may also specify exact parameters or use a graphical interface to guess parameters in the traditional way if desired. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gridsample </td>
   <td style="text-align:left;max-width: 100; "> Multi-stage cluster surveys of households are commonly performed by
  governments and programmes to monitor population-level demographic, social,
  economic, and health outcomes. Generally, communities are sampled from
  subpopulations (strata) in a first stage, and then households are listed and
  sampled in a second stage. In this typical two-stage design, sampled
  communities are the Primary Sampling Units (PSUs) and households are the
  Secondary Sampling Units (SSUs). Census data typically serve as the sample
  frame from which PSUs are selected. However, if census data are outdated
  inaccurate, or too geographically course, gridded population data (such as
  &lt;http://www.worldpop.org.uk&gt;) can be used as a sample frame instead.
  GridSample (&lt;doi:10.1186/s12942-017-0098-4&gt;) generates PSUs from
  gridded population data according to user-specified complex survey design
  characteristics and household sample size. In gridded population sampling,
  like census sampling, PSUs are selected within each stratum using a
  serpentine sampling method, and can be oversampled in urban or rural areas to
  ensure a minimum sample size in each of these important sub-domains.
  Furthermore, because grid cells are uniform in size and shape, gridded
  population sampling allows for samples to be representative of both the
  population and of space, which is not possible with a census sample frame. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hipread </td>
   <td style="text-align:left;max-width: 100; "> Read hierarchical fixed width files like those commonly used by 
    many census data providers. Also allows for reading of data in chunks,
    and reading 'gzipped' files without storing the full file in memory. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> idbr </td>
   <td style="text-align:left;max-width: 100; "> Use R to make requests to the US Census Bureau's International Data Base API.
             Results are returned as R data frames.  For more information about the IDB API, visit
             &lt;http://www.census.gov/data/developers/data-sets/international-database.html&gt;. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> inca </td>
   <td style="text-align:left;max-width: 100; "> Specific functions are provided for rounding real weights to integers and performing an integer programming algorithm for calibration problems. They are useful for census-weights adjustments, or for performing linear regression with integer parameters. This research was supported in part by the U.S. Department of Agriculture, National Agriculture Statistics Service. The findings and conclusions in this publication are those of the authors and should not be construed to represent any official USDA, or US Government determination or policy. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ipfr </td>
   <td style="text-align:left;max-width: 100; "> Performs iterative proportional updating given a seed table and
  an arbitrary number of marginal distributions. This is commonly used in
  population synthesis, survey raking, matrix rebalancing, and other
  applications. For example, a household survey may be weighted to match the
  known distribution of households by size from the census. An origin/
  destination trip matrix might be balanced to match traffic counts.
  The approach used by this package is based on a paper from
  Arizona State University (Ye, Xin, et. al. (2009)
  &lt;http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.537.723&amp;rep=rep1&amp;type=pdf&gt;).
  Some enhancements have been made to their work including primary and 
  secondary target balance/importance, general marginal agreement, and weight 
  restriction. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ipumsr </td>
   <td style="text-align:left;max-width: 100; "> An easy way to import census, survey and geographic data provided by 'IPUMS'
    into R plus tools to help use the associated metadata to make analysis easier. 'IPUMS'
    data describing 1.4 billion individuals drawn from over 750 censuses and surveys is
    available free of charge from our website &lt;https://ipums.org&gt;. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PakPC2017 </td>
   <td style="text-align:left;max-width: 100; "> Provides data sets and functions for exploration of Pakistan Population Census 2017 (&lt;http://www.pbscensus.gov.pk/&gt;). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ppcSpatial </td>
   <td style="text-align:left;max-width: 100; "> Spatial Analysis for exploration of Pakistan Population Census 2017 (&lt;http://www.pbscensus.gov.pk/&gt;). It uses data from R package 'PakPC2017'. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RapidPolygonLookup </td>
   <td style="text-align:left;max-width: 100; "> Facilitates efficient polygon search using kd trees.
    Coordinate level spatial data can be aggregated to higher geographical
    identities like census blocks, ZIP codes or police district boundaries.
    This process requires mapping each point in the given data set to a
    particular identity of the desired geographical hierarchy. Unless efficient
    data structures are used, this can be a daunting task. The operation
    point.in.polygon() from the package sp is computationally expensive.
    Here, we exploit kd-trees as efficient nearest neighbor search algorithm
    to dramatically reduce the effective number of polygons being searched. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> revengc </td>
   <td style="text-align:left;max-width: 100; "> Decoupled (e.g. separate averages) and censored (e.g. &gt; 100 species) variables are continually reported by many well-established organizations (e.g. World Health Organization (WHO), Centers for Disease Control and Prevention (CDC), World Bank, and various national censuses).  The challenge therefore is to infer what the original data could have been given summarized information.  We present an R package that reverse engineers decoupled and/or censored count data with two main functions.  The cnbinom.pars function estimates the average and dispersion parameter of a censored univariate frequency table.  The rec function reverse engineers summarized data into an uncensored bivariate table of probabilities. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> satscanMapper </td>
   <td style="text-align:left;max-width: 100; "> Supports the generation of maps based on the results from 
     'SaTScan' (TM) cluster analysis.
     The package handles mapping of Spatial and Spatial-Time analysis using
     the discrete Poisson, Bernoulli, and exponential models of case data generating
     cluster and location ('GIS') records containing observed, expected and observed/expected
     ratio for U. S. states (and DC), counties or census tracts of individual 
     states based on the U. S. 'FIPS' codes for state, county and census tracts 
     (locations) using 2000 or 2010 Census areas, 'FIPS' codes, and boundary data.
     'satscanMapper' uses the 'SeerMapper' package for the boundary data and 
     mapping of locations.  Not all of the 'SaTScan' (TM) analysis and models generate
     the observed, expected and observed/expected ratio values for the clusters and 
     locations.
     The user can map the observed/expected ratios for locations 
     (states, counties, or census tracts) for each cluster with a p-value less than 0.05 
     or a user specified p-value.  
     The locations are categorized and colored based on either the cluster's Observed/Expected 
     ratio or the locations' Observed/Expected ratio. 
     The place names are provided for each census tract using data from 'NCI', the 'HUD' crossover 
     tables (Tract to Zip code) as of December, 2013, the USPS Zip code 5 database for 1999, 
     and manual look ups on the USPS.gov web site. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SeerMapper </td>
   <td style="text-align:left;max-width: 100; "> Provides an easy way to map seer registry area rate data on a U. S, map.  
   The U. S. data may be mapped at the state, U. S. NCI Seer Register, state/county 
   or census tract level. The function can categorize the data into "n" quantiles, where "n" is 3 to 11 or
   the caller can specify a cut point list for the categorizes.  
   The caller can also provide the data and the comparison operation to request
   hatching over any areas.  The default operation and value are &gt; 0.05 (p-values).
   The location id provided in the data determines the geographic level of the mapping.
   If states, state/counties or census tracts are being mapped, the location ids 
   used must be the U.S. FIPS codes for states (2 digits), state/counties (5 digits)
   or state/county/census tracts (11 digits). If the location id references the U.S. Seer Registry 
   areas, the Seer Registry area identifier used to link the data to the geographical 
   areas, then the location id is the Seer Registry name or abbreviation.
   Additional parameters are used to provide control over the drawing of the boundaries
   at the data's boundary level and higher levels.
   The package uses modified boundary data from the 2000 and 2010 U. S. Census to reduce the 
   storage requirements and improve drawing speed.  
   The 'SeerMapper' package contains the U. S. Census 2000 and 2010 boundary data
   for the regional, state, Seer Registry, and county levels.  Six supplement packages 
   contain the census tract boundary data (see manual for more details.) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SeerMapper2010East </td>
   <td style="text-align:left;max-width: 100; "> Provides supplemental 2010 census tract boundary package for 23 states
   without Seer Registries that are east of the Mississippi river 
   for use with the 'SeerMapper' package.  
   The data contained in this 
   package is derived from U. S. Census data and is in public domain. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SeerMapper2010Regs </td>
   <td style="text-align:left;max-width: 100; "> Provides  supplemental 2010 census tract boundaries of the 15 states 
   containing Seer Registries for use with the 'SeerMapper' package.
   The data contained in this 
   package is derived from U. S. 2010 Census data and is in public domain. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SeerMapper2010West </td>
   <td style="text-align:left;max-width: 100; "> Provides supplemental 2010 census tract boundaries for the 14 states
   without Seer Registries that are west of the Mississippi river 
   for use with the 'SeerMapper' package.
   The data contained in this 
   package is derived from U. S. 2010 Census data and is in public domain. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SeerMapperEast </td>
   <td style="text-align:left;max-width: 100; "> Provides supplemental 2000 census tract boundaries for the 23 states
   without Seer Registries that are east of the Mississippi river 
   for use with the 'SeerMapper' package.  
   The data contained in this 
   package is derived from U. S. Census data and is in the public domain. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SeerMapperRegs </td>
   <td style="text-align:left;max-width: 100; "> Provides supplemental 2000 census tract boundaries for the 15 states
   containing Seer Registries for use with the 'SeerMapper' package.  
   The data contained in this 
   package is derived from U. S. Census data and is in the public domain. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SeerMapperWest </td>
   <td style="text-align:left;max-width: 100; "> Provides supplemental 2000 census tract boundaries for the 14 states
   without Seer Registries that are west of the Mississippi river
   for use with the 'SeerMapper' package.  
   The data contained in this 
   package is derived from U. S. Census data and is in the public domain. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sms </td>
   <td style="text-align:left;max-width: 100; "> Produce small area population estimates by fitting census data to
    survey data. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SpatialVS </td>
   <td style="text-align:left;max-width: 100; "> Perform variable selection for the spatial Poisson regression model under the adaptive elastic net penalty. Spatial count data with covariates is the input. We use a spatial Poisson regression model to link the spatial counts and covariates. For maximization of the likelihood under adaptive elastic net penalty, we implemented the penalized quasi-likelihood (PQL) and the approximate penalized loglikelihood (APL) methods. The proposed methods can automatically select important covariates, while adjusting for possible spatial correlations among the responses. More details are available in Xie et al. (2018, &lt;arXiv:1809.06418&gt;). The package also contains the Lyme disease dataset, which consists of the disease case data from 2006 to 2011, and demographic data and land cover data in Virginia. The Lyme disease case data were collected by the Virginia Department of Health. The demographic data (e.g., population density, median income, and average age) are from the 2010 census. Land cover data were obtained from the Multi-Resolution Land Cover Consortium for 2006. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TexMix </td>
   <td style="text-align:left;max-width: 100; "> A collection of functions and data - mostly from Texas - is provided. 
  These are used as teaching tools for geo-spatial data analytics courses in the 
  GISciences program at The University of Texas at Dallas. In addition, several 
  vignettes illustrate geo-spatial data analytics practices, such as relative 
  risk kernel density estimations based on food store locations within Dallas 
  County or the identification of homogenous and spatially contiguous market areas 
  built on socio-economic, demographic and infrastructure census information. The 
  spatial resolution of the data-sets ranges from 1623 food store locations, over 
  geo-referenced areal data of 258 Texan counties, to 529 census tracts as well 
  as 1669 block groups in Dallas County. Cartographic, specialized regression and 
  data handling functions are provided. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyqwi </td>
   <td style="text-align:left;max-width: 100; "> The purpose of this package is to access the 
    United States Census Bureau's Quarterly Workforce Indicator data. Additionally, 
    the data will be retrieved in a tidy format for further manipulation with full variable
    descriptions added if desired. Information about the United States Census Bureau's 
    Quarterly Workforce Indicator is available at 
    &lt;https://www.census.gov/data/developers/data-sets/qwi.html&gt;. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyUSDA </td>
   <td style="text-align:left;max-width: 100; "> Provides a consistent API to pull United States
    Department of Agriculture census and survey data from the National
    Agricultural Statistics Service (NASS) QuickStats service
    &lt;https://quickstats.nass.usda.gov&gt;. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tigris </td>
   <td style="text-align:left;max-width: 100; "> Download TIGER/Line shapefiles from the United States Census Bureau 
    (&lt;https://www.census.gov/geo/maps-data/data/tiger-line.html&gt;) and load into R as 'SpatialDataFrame' or 'sf' objects. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> totalcensus </td>
   <td style="text-align:left;max-width: 100; "> Download summary files from Census Bureau &lt;https://www2.census.gov/&gt; 
    and extract data, in particular high resolution data at 
    block, block group, and tract level, from decennial census and 
    American Community Survey 1-year and 5-year estimates. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UScancer </td>
   <td style="text-align:left;max-width: 100; "> This package contains functions to read cancer data from SEER (http://seer.cancer.gov/) and IARC (http://www.iarc.fr) to create datasets at the county level based on US census information. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UScensus2010 </td>
   <td style="text-align:left;max-width: 100; "> US Census 2010 shape files and additional demographic data
        from the SF1 100 percent files. This package contains a number
        of helper functions for the UScensus2010blk,
        UScensus2010blkgrp, UScensus2010tract, UScensus2010cdp
        packages. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> valetr </td>
   <td style="text-align:left;max-width: 100; "> Interface to Bank of Canada's 'Valet' API (&lt;https://www.bankofcanada.ca/valet/docs&gt;). Please read the API terms and conditions: &lt;https://www.bankofcanada.ca/terms/&gt;. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x12 </td>
   <td style="text-align:left;max-width: 100; "> The 'X13-ARIMA-SEATS' &lt;https://www.census.gov/srd/www/x13as/&gt; methodology and software is a widely used software and developed by the US Census Bureau. It can be accessed from 'R' with this package and 'X13-ARIMA-SEATS' binaries are provided by the 'R' package 'x13binary'. </td>
  </tr>
</tbody>
</table>
