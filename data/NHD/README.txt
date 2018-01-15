NHD: Data from US National Hydrologic Dataset

HYD = National Hydrography Dataset: https://nhd.usgs.gov/data.html
    Contents:
        NHDFlowline_SouthForkEel.shp = NHDFlowline for South Fork Eel watershed
        NHDFlowline_SouthForkEel+Adjacent = NHDFlowline for South Fork Eel watershed + adjacent watersheds on all sides
    These are all extracted from the NHD_H_1801_HU4_Shape dataset downloaded from ( http://prd-tnm.s3-website-us-west-2.amazonaws.com/?prefix=StagedProducts/Hydrography/NHD/HU4/HighResolution/Shape/ )
    which is stored on GSAS at ( Z:\2.active_projects\Zipper\1.Spatial_data\regional\SouthForkEel\riv_river_network_streamflow\1original\NHD_H_1801_HU4_Shape )

WBD = Watershed Boundary Dataset (also known as HUC boundaries): https://nhd.usgs.gov/wbd.html
    Contents:
        WBDHU6_NorthCaliforniaCoast.shp = HUC6 watershed containing South Fork Eel
        WBDHU8_SouthForkEel.shp = South Fork Eel watershed
        WBDHU8_SouthForkEel+Adjacent.shp = South Fork Eel watershed + adjacent watersheds on all sides
        WBDHU10_SouthForkEel.shp = Subwatersheds within South Fork Eel
        WBDHU12_SouthForkEel.shp = Sub-subwatersheds within South Fork Eel
    These are all extracted from the WBD_18_HU2_Shape dataset downloaded from ( ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/Hydrography/WBD/HU2/Shape/ )
    which is stored on GSAS at ( Z:\2.active_projects\Zipper\1.Spatial_data\regional\SouthForkEel\ws_watersheds\1original )