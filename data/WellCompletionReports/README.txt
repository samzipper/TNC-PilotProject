Navarro_WellCompletionReports shapefile

Well completion reports downloaded 2018-11-23 from CA DWR Well Completion Report Map Application
https://dwr.maps.arcgis.com/apps/webappviewer/index.html
and trimmed to Navarro River Watershed + Adjacent HUC12s using QGIS.

Some potentially important fields might be:
-Decimal latitude = latitude at centroid of township/range/section quandrangle
-Decimal longitude = latitude at centroid of township/range/section quandrangle
-Township
-Range
-Section
-Planned use/former use = use of well
-Top of Perforated Interval = depth to top of well screen [ft]
-Bottom of Perforated Interval = depth to bottom of well screen [ft]

Navarro_WellCompletionReports_BedrockInfo.csv = wells from Navarro_WellCompletionReports.shp file that have both the screen depth
and a PDF well completion report available online, with some information manually entered from the PDFs:
-WeatheredDepthToBedrock_ft = depth of first layer from driller log that corresponds to bedrock (sandstone, rock, shale, etc) (9999 = did not reach bedrock)
-DepthToBedrock_ft = depth of first layer from driller log that corresponds to bedrock (sandstone, rock, shale, etc), ignoring weathered bedrock (9999 = did not reach bedrock)
-BedrockType_WellScreen = description of bedrock materials in screened interval
(Note: I screened these while looking at the PDFs to eliminate ones with errors, for example where the shapefile well screen depth did not match the PDF)

Raw well completion reports for all of Mendocino County are stored on GSAS
Z:\2.active_projects\Zipper\1.Spatial_data\regional\NavarroRiver\use_withdrawl_abstraction_use\1original\WellCompletionReports