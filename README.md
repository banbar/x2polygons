This Python package allows the calculation of distance between two polygons. Polygons may have thematic attributes (e.g. name) to add further description.

The available geometric distance functions between polygons are:
* Hausdorff distance
* Chamfer distance
* PoLis distance
* Turn function distance

The avaiable thematic distance functions are:
* Levenshtein distance


For Windows:
1. Follow the steps described [here](https://stackoverflow.com/a/58943939/1959766) to download geopandas.
2. Resolve the *geos_c.dll file missing* error by running *conda install shapely* as described [here](https://github.com/Toblerity/Shapely/issues/1032).
