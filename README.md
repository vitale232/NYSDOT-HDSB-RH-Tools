# NYSDOT-HDSB-RH-Tools
Utilities to accompany NYSDOT's Esri Roads and Highways implementation (ArcGIS 10.5.1). The tools have three general categories:

|Category|Description|
|--------|-----------|
|Calibration|QA Tools to ensure routes are calibrated properly|
|Database|Tools to generate new Route IDs and DOT IDs from the database sequences|
|Signed Routes|Tools to work with gapped routes and dual-carriageway roads|


# Included Tools

## Calibration
#### Compare Calibrated Length and GIS Length
This tool checks for valid calibration of an Esri Roads and Highways route. It assumes a valid calibration uses the geographic length to calibrate the measure values, calibrates the route in the "drawing direction", and that gapped routes use the "Add Euclidean Distance" method of gap calibration. The tool will add four new fields to the output feature class's attribute table:

1. GIS_Length - This field will provide the calculated 2 dimensional length of the route's shape field. The length is only valid if the data are in a projected coordinate system that can be represented by a Cartesian Plane (e.g. Universal Transverse Mercator). If the route is multiple parts, the euclidean distance of the gap(s) between parts will be added to the GIS_Length. This value represents the hypothetical calibrated length of the route or routes.
2. M_Length - This field provides the maximum "m-value" stored on the route vertices. This value is the current calibrated length of the route in the R&H platform.
3. length_diff - This value is a simple difference of the M_length and the GIS_Length. It is calculated as follows: `length_diff = M_Length - GIS_Length`
4. part_count - This field will display the total number of polylines that make up a route. It is simply the `partCount` attribute of the shape field.

#### Generate Equal Interval Points
This tool creates two point feature classes. One point feature class will represent equal interval points while traversing the geographic length of the route. For every given `interval` of 2 dimensional geographic length, a point will be produced. The other feature class will be euqal interval points while traversing the "m" values of the line. Each point represents a given length `interval` between the previously calculated length.

The output of this tool is most useful when used in tandom. If the two point feature calsses align with one another, the route should be calibrated correctly. If the points start to diverge from one another, you are likely in an area that led to a miscalibration of the route. The centerline feature class and routes in this area should be examined closesly.

*Note - This tool has a known bug on gapped routes. The points do not correctly account for the gap length. The first part of the line feature will have reliable results, but parts 2+ will not.*

## Database
#### Generate New DOT_ID and/or ROUTE_ID
This tool will look in your ArcCatalog "Database Connections" folder for a file matching either `dev_elrs_ad_Lockroot.sde` or `dev_elrs_ad_lockroot.sde`, which connect to the R&H development environment. The database in the R&H dev environment is set up with two sequences to track new DOT_IDs and ROUTE_IDs. This tool calls those sequences to generate new values, and presents them to the user in the output messages. There are input options to generate either a DOT_ID, ROUTE_ID, or both.

# Signed Routes
#### Calculate Euclidean Distance of Gapped Routes
This tool only works on the selected routes, and calculates the Euclidean Distance of any gaps that are present. The output is presented in a table.

#### Identify Primary and Reverse Divergence
This tool essentially ensures that a dual carriageway geometry (i.e. a reverse route) is required for the input DOT_ID and COUNTY_ORDER. If the output geometry from this tool is empty, the primary and reverse route represented by the input DOT ID and COUNTY_ORDER are fully overlapping. If there are feature geometries in the output feature class, they represent the areas that are truly divided roadways.

