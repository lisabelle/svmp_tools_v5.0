import arcpy
import numpy as np

"""
svmpUtils.py

This module contains shared functions and constants that are used
with the SVMP Site Transect Analysis Tools

This script was designed and developed by Allison Bailey, Sound GIS.
Developed in ArcGIS Pro provided Python 3.7.11, NumPy 1.2.3, Pandas 1.20.1
ArcGIS Pro 2.9.5
6/30/2023
"""

#---------------------------- CONSTANTS -----------------------------------#
# Columns
vegcodeCol = 'veg_code'
visityearCol = 'visit_year'
sampselCol = 'samp_sel'
sampstatCol = 'samp_status'
studycodeCol = 'study_code'
sitevisitidCol = 'site_visit_id'
sitecodeCol = 'site_code'
sampidCol = 'site_samp_id'
datesampCol = 'date_samp_start'
transectidCol = 'transect_id'
surveyidCol = 'survey_id'
surveystatCol = 'survey_status'
maxdepflagCol = 'tran_maxd_qual'
mindepflagCol = 'tran_mind_qual'
datetimesampCol = 'date_time_samp'
depInterpCol = 'depth_interp'
videoCol = 'video'
ptidCol = 'ID'

#Tables
vegcodesTbl = 'veg_codes'
sitevisitsTbl = 'site_visits'
sitesamplesTbl = 'site_samples'
studyassociationsTbl = 'study_associations'
transectsTbl = 'transects'
surveysTbl = 'surveys'
vegoccurTbl = 'veg_occur'
segmentsTbl = 'segments'

#Feature Classes
samppolyFC = 'samp_polygons'

NULL_DEPTH = -9999
NULL_VEG = -9999
NULL_VIDEO = -9999
NULL_VAR = -9999

## Transect_results table fields
transect_results_fields = [
    "tran_results_id",
    "transect_id",
    "veg_code",
    "tran_len_ft",
    "veg_len_ft",
    "tran_veg_frac",
    "tran_veg_maxd_ft",
    "tran_maxd_ft",
    "tran_veg_mind_ft",
    "tran_mind_ft",
    "tran_maxd_qual",
    "tran_mind_qual",
    "site_results_id"
]

# site_results table fields
site_results_fields = [
    "site_results_id",
    "site_samp_id",
    "veg_code",
    "veg_area_n_tran",
    "veg_frac",
    "samp_area_ft2",
    "veg_area_ft2",
    "veg_area_se_ft2",
    "veg_mind_n_tran",
    "veg_mind_mean_ft",
    "veg_mind_deepest_ft",
    "veg_mind_shallowest_ft",
    "veg_mind_se_ft",
    "veg_maxd_n_tran",
    "veg_maxd_mean_ft",
    "veg_maxd_deepest_ft",
    "veg_maxd_shallowest_ft",
    "veg_maxd_se_ft"
]


# Zero and no data values for site_results table
site_results_zero = [
    0,
    0,
    0,
    0,
    0,
    0,
    -9999,
    -9999,
    -9999,
    -9999,
    0,
    -9999,
    -9999,
    -9999,
    -9999,
]



# Spatial Reference
sr = arcpy.SpatialReference(2927) # NAD_1983_HARN_StatePlane_Washington_South_FIPS_4602_Fee

# Filter for Sample Status -- only these are used for statistics calcs
sampstat4stats = ["sampled","exception"]


def unique_values(table, field):
    # Get list of all unique values in a field
    # search cursor wrapped in list generator creating list of all values
    values = (row[0] for row in arcpy.da.SearchCursor(table, (field)))
    # pass list into set to get only unique values and return the result
    return sorted(set(values))

def unique_values_np(table, field):
    data = arcpy.da.TableToNumPyArray(table, [field])
    return np.unique(data[field])

def tables_fcs_list(gdb):
    # Get list of tables and feature classes in a geodatabase
    tables = []
    fcs = []
    # -------- Initial Code --------------
    # #Save initial state of workspace
    # initial_ws = arcpy.env.workspace
    #
    # #change workspace to input geodatabase
    # arcpy.env.workspace = gdb
    # tables = arcpy.ListTables()
    # fcs = arcpy.ListFeatureClasses()
    # #reset workspace back to original
    # arcpy.env.workspace = initial_ws
    #----------------------------------------------

    # --------  New approach ------------------------------------
    #------  Increased speed by ~38 seconds ------------------------------
    for path, path_names, data_names in arcpy.da.Walk(gdb, "FeatureClass"):
        for data_name in data_names:
            fcs.append(data_name)

    for path, path_names, data_names in arcpy.da.Walk(gdb, "Table"):
        for data_name in data_names:
            tables.append(data_name)
    #------  Increased speed by ~38 seconds ------------------------------

    return {"tables":tables,"fcs":fcs,"both":tables + fcs}


def fieldExists(dataset, field_name):
    if field_name in [field.name for field in arcpy.ListFields(dataset)]:
        return True



def stdDev(sample):
    """ Calculates Standard Deviation of a sample
    :param sample: list of values that make up the sample
    :return: standard deviation
    """
    N = len(sample)  # number of samples
    # Can't calculate if number if samples is 0 or 1
    mean = float(sum(sample)) / N  # mean of the samples
    sum_sqdif = 0  # initialize sum of squared differences
    # Calculate sum of squared differences
    for val in sample:
        sqdif = (val - mean) ** 2
        sum_sqdif = ((val - mean) ** 2) + sum_sqdif
    # Standard Deviation
    s = ((1 / (float(N) - 1)) * sum_sqdif) ** 0.5
    return s


def variance(stdDev):
    """ Calculates variance of a sample
    :param stdDev:  standard deviation
    :return: standard deviation
    """
    var = stdDev ** 2
    return var


def stdErr(s, N):
    """ Calculate Standard Error of a sample
    :param s: Standard Deviation
    :param N: number of samples
    :return: standard error
    """
    SE = s / (float(N) ** 0.5)
    return SE


def ci95(SE):
    """ Calculate 95% confidence interval
    :param SE: standard error
    :return: 95% confidence interval
    """
    confint95 = 1.96 * SE
    return confint95