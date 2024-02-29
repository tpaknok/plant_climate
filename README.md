# plant_climate
Using GBIF data to understand habitat usage of plants under different climatic conditinos. Data are too large and thus are not uploaded. See Script 01 for code used to download data from GBIF.

Script 02 are for SQL operations to select useful data only. Data are too large to be processed by R. 

It is possible to use more SQL operations (e.g., extract unique species names first, combine synonyms into the same species and remove species with insufficent records) to further clean the data before cleaning the coordinates. However, this is not computationally efficient because too many names are needed for fuzzy matching.

Thus, only a simple SQL operation is written, and the data is imported into the R environment for coordinate cleaning (requires a R function to do it)


