
import time
import datetime
import os
import csv
import fnmatch
import re
import numpy as np
import pandas as pd
import arcpy

"""
csv2pt.py

Convert a CSV file to a point feature class using pandas and NumPy

This script was designed and developed by Allison Bailey, Sound GIS.
All updates after June 2018 by Lisa Ferrier.
Developed in ArcGIS Pro provided Python 3.7.11, NumPy 1.20.1, Pandas 1.2.3
ArcGIS Pro 2.9.5
7/12/2023
"""

__version__= '5.0'

# import svmpUtils as utils

t0 = time.time()

# ------------- VARIABLES Related to Source ASCII data files ------------- #
# Source ASCII data are provided by MRC in a standard comma-delimited format
# The first row contains the column names,
# but the order of columns within the file may vary

# Columns needed from source csv file
sourceSiteCol = 'Site'  # site code e.g. core001
sourceTrkCol = 'trk'  # column to identify a track/transect
sourceDateCol = 'date'  # column name for date
sourceTimeCol = 'time'  # column for time
sourceDepthObsCol = 'BSdepth'  # observed depth
sourceDepthInterpCol = 'BSdepth_interp'  # interpolated depth
sourceVideoCol = 'video'  # video data quality (0,1)
sourceLatCol = 'latitude'  # column name for latitude
sourceLonCol = 'lon'  # column name for longitude

# Columns required to be input csv file
REQD_COLUMNS =[
    sourceSiteCol,
    sourceTrkCol,
    sourceDateCol,
    sourceTimeCol,
    sourceDepthObsCol,
    sourceDepthInterpCol,
    sourceVideoCol,
    sourceLatCol,
    sourceLonCol,
]

# ------------- VARIABLES Related to Output Point Feature Classes ------------- #
ptidCol = 'ID'
surveyidCol = 'survey_id'  # Unique survey ID = site_code + date of transect + transect number as string
siteCol = 'site_code' # site code e.g. core001
tranCol = 'tran_num' # column to identify a track/transect
dateCol = 'date_svy' # Date of survey
timeCol = 'time_svy' # Time of survey
datetimeCol = 'date_time_samp'  # Column with time stamp (date/time) of point location
depObsCol = 'depth_obs' # observed depth
depInterpCol = 'depth_interp' # interpolated depth
videoCol = 'video'  # Column for video data quality (0,1)
latCol = sourceLatCol # column name for latitude
lonCol = sourceLonCol # column name for longitude
vegcodeCol = 'veg_code'
nativesgCol = 'nativesg'


# Columns that are used to calculate nativesg attribute
NATIVESG_CODES = [
            'Zm',
            'Phyllo',
            'undiff',
            ]

NULL_DEPTH = -9999
NULL_VEG = -9999
NULL_VIDEO = -9999

# Convert Degree/minute format to Decimal Degrees
# Assumes input in format ##d##.####' Dir  (i.e. 48d33.8342' N)
def dm2dd(coordDegMin):
    [deg,minDir] = coordDegMin.split('d')
    [minute,Dir] = minDir.split("' ")
    coordDD = int(deg) + (float(minute) / 60)
    if Dir in ['W','S']:
        coordDD = coordDD * (-1)
    return coordDD

class CsvPath(object):
    """ Represents a directory path for a single site

    Properties:
    sitecode -- the code for the site
    csvdir -- the full path to the directory -- concatenation of the base directory and site code
    search_pattern -- the pattern to match for transect data files:   sitecode_YYYY_##_TD.csv
    valid -- flag to indicate if path is valid -- directory exists and their are transect files
    dir_exist -- flag to indicate if the directory exists
    files_exist -- flag to indicate if files matching the search pattern exist in the directory
    tdfiles -- list of transect data files within the directory

    """

    def __init__(self, sitecode, basedir):
        self.sitecode = sitecode
        self.csvdir = os.path.normpath(os.path.join(basedir, sitecode))
        self.search_pattern = self.sitecode + '_*_*_TD.csv'

    @property
    def valid(self):
        """Binary attribute indicating the existence of the csv path"""
        if self.dir_exists and self.files_exist:
            return True
        else:
            return False

    @property
    def dir_exists(self):
        """Binary attribute indicating the existence of the csv path"""
        if os.path.exists(self.csvdir):
            return True
        else:
            return False
    @property
    def files_exist(self):
        if self.tdfiles:
            return True
        else:
            return False

    @property
    def tdfiles(self):
        """List of transect data files matching the specified pattern within the directory"""
        _tdfiles = []
        if self.dir_exists:
            files = os.listdir(self.csvdir)
            _tdfiles = fnmatch.filter(files, self.search_pattern)
        return _tdfiles


class CsvSource(object):
    """ Represents a csv source data file for a single site visit (site, year, group)

    Properties
    file_path -- full path to the csv file
    valid -- binary to flag the overall validity of the csv file
    file_exists -- binary attribute indicating the existence of the csv file
    rows -- all rows from the csv file as extracted via DictReader
    all_columns -- list of all columns in the source csv file
    base_columns -- list of required columns that exist in the source csv file
    missing_columns -- list of required columns that are missing from the csv
    veg_columns -- list of vegetation columns in source csv file
    columns -- combined list of base columns and veg columns in source csv file
    lat_errors -- list of csv rows with latitude values that don't match the expected pattern
    lon_errors -- list of csv rows with longitude values that don't match the expected pattern
    dataframe -- a pandas dataframe of the csv source data

    """

    def __init__(self, file_path, veg_codes):
        self.file_path = file_path
        self.reqd_columns = REQD_COLUMNS
        self.sourceLatCol = sourceLatCol
        self.sourceLonCol = sourceLonCol
        self.sourceTimeCol = sourceTimeCol
        self.sourceDateCol = sourceDateCol
        self.veg_codes = veg_codes

    @property
    def valid(self):
        """ property to flag the overall validity of the csv file, specifically,
            File must exist
            No missing base columns
            At least one vegetation column
            No malformed latitude or longitude values
            No malformed date or time values
        """
        if self.file_exists and len(self.missing_columns) == 0 and len(self.veg_columns) > 0 and \
                len(self.lat_errors) == 0 and len(self.lon_errors) == 0 and len(self.time_errors) == 0 and len(self.date_errors) == 0:
            return True
        else:
            return False

    @property
    def file_exists(self):
        """Binary attribute indicating the existence of the csv file"""
        if os.path.exists(self.file_path):
            return True
        else:
            return False

    @property
    def rows(self):
        """a text-based csv reader with all rows from the source csv file"""
        _rows = []
        try:
            csv_file = open(self.file_path, newline="")
            csv_rows = csv.DictReader(csv_file)
            for row in csv_rows:
                _rows.append(row)
            csv_file.close()
            return _rows
        except:
            return _rows

    @property
    def all_columns(self):
        """list of all columns in the source csv file"""
        try:
            csv_file = open(self.file_path, newline="")
            csv_rows = csv.DictReader(csv_file)
            _all_columns = csv_rows.fieldnames
            csv_file.close()
            return _all_columns
        except:
            return []

    @property
    def base_columns(self):
        """list of required columns that are present in the source csv"""
        _base_columns = set(self.all_columns).intersection(set(self.reqd_columns))
        return list(_base_columns)

    @property
    def missing_columns(self):
        """list of required columns that are missing from the source csv"""
        _missing_columns = set(self.reqd_columns).difference(set(self.all_columns))
        return list(_missing_columns)

    @property
    def veg_columns(self):
        """list of vegetation columns in source csv file"""
        _veg_columns = set(self.veg_codes).intersection(set(self.all_columns))
        return list(_veg_columns)

    @property
    def columns(self):
        """combined list of base and vegetation columns in source csv file"""
        _columns = self.base_columns + self.veg_columns
        return _columns

    @property
    def lat_errors(self):
        """list of rows in source csv file with erroneous latitude values"""
        try:
            _lat_errors = self._validate_latlon(self.sourceLatCol)
            return _lat_errors
        except:
            return None

    @property
    def lon_errors(self):
        """list of rows in source csv file with erroneous longitude values"""
        try:
            _lon_errors = self._validate_latlon(self.sourceLonCol)
            return _lon_errors
        except:
            return None

    @property
    def time_errors(self):
        """list of rows in source csv file with erroneous time values"""
        try:
            _time_errors = self._validate_time(self.sourceTimeCol)
            return _time_errors
        except:
            return None

    @property
    def date_errors(self):
        """list of rows in source csv file with erroneous date values"""
        try:
            _date_errors = self._validate_date(self.sourceDateCol)
            return _date_errors
        except:
            return None

    def _validate_latlon(self, col):
        """Function to validate lat or lon values based on expected pattern"""
        # Assumes input in format ##d##.####' Dir or ##d##' Dir
        # latlon_pattern = "[0-9]+[d][0-9]+[.][0-9]+\' [nsewNSEW]"
        latlon_pattern = "[0-9]+[d][0-9]+[.]?[0-9]+\' [nsewNSEW]"
        error_rows = [] # initialize list of rows with errors
        # Loop through data and validate lat/long value
        for i, row in enumerate(self.rows):
            csv_row = i + 1
            coord = row[col]
            if not re.search(latlon_pattern, coord):
                error_rows.append(csv_row)
        #
        # """
        # ALTERNATE UNUSED APPROACH, but kept as an example for reference
        # Pandas approach: Creates a True/False series showing records that match the pattern
        # The tilde (~) in front of series inverts the True/False values,
        # which is the same as records that DO NOT match the pattern
        # """
        # validrows = self.dataframe[col].str.contains(latlong_pattern)
        # errors = self.dataframe[~validrows]

        return error_rows

    def _validate_time(self, col):
        """
        Function to validate time values based on expected patter
        Input can be in 12-hour (AM/PM) or 24-hour format
        :param col:
        :return:
        """
        error_rows = [] # initialize list of rows with errors
        # Loop through data and validate time values
        for i, row in enumerate(self.rows):
            csv_row = i + 1
            time_of_survey = row[col]
            time24hr_pattern = "^(2[0-3]|[01]?[0-9]):([0-5]?[0-9]):([0-5]?[0-9])$"
            time12hr_pattern = "^(1[0-2]|0?[1-9]):([0-5]?[0-9]):([0-5]?[0-9])( ?[AP]M)?$"

            if "M" in time_of_survey:
                if not re.search(time12hr_pattern, time_of_survey):
                    error_rows.append(csv_row)
            else:
                if not re.search(time24hr_pattern, time_of_survey):
                    error_rows.append(csv_row)
        return error_rows

    def _validate_date(self, col):
        """
        Function to validate date values based on expected pattern
        Input expected to be in this format:  m/d/yyyy or mm/dd/yyyy
        :param col:
        :return:
        """
        error_rows = [] # initialize list of rows with errors
        # Loop through data and validate time values
        for i, row in enumerate(self.rows):
            csv_row = i + 1
            date_of_survey = row[col]
            try:
                [m, d, y] = date_of_survey.split('/')
                testdate = datetime.date(int(y), int(m), int(d))
            except:
                error_rows.append(csv_row)
        return error_rows

    @property
    def dataframe(self):
        if self.valid:

            # usecols=self.columns
            print(self.columns)
            _dataframe = pd.read_csv(self.file_path,
                                     memory_map=True,
                                     converters={
                                         sourceLatCol: dm2dd,
                                         sourceLonCol: dm2dd,
                                     },
                                     skipinitialspace=True,
                                     parse_dates={
                                         datetimeCol: [sourceDateCol, sourceTimeCol]
                                     },
                                     keep_date_col=True
                                     )
            print(_dataframe.columns)
            return _dataframe
        else:
            return None


class CsvData(object):
    """ Represents the data for a single site visit

    Properties:
    csv_source -- the input CsvSource object
    source_columns -- the relevant columns in the csv source that are needed for data processing
    source_veg_columns -- list of vegetation columns in source csv file
    veg_columns -- vegetation columns -- may get updated if nativesg is added
    df -- pandas dataframe of the csv source data
    nparray -- structured NumPy array created from the pandas dataframe

    """

    def __init__(self, csv_source):
        # Get some properties from the csv_source object
        self.csv_source = csv_source
        self.source_columns = self.csv_source.columns
        self.veg_columns = self.csv_source.veg_columns
        self.source_veg_columns = self.csv_source.veg_columns
        self.df = self.csv_source.dataframe
        # Initialize warnings property
        self.warnings = False
        # process the data
        self._process_data()
        # validate the data
        self._validate_data()


    @property
    def base_dtype(self):
        """" NumPy Data types for base columns """
        _base_dtype = [
                        (ptidCol, '<i8'),
                        (surveyidCol, 'S25'),
                        (siteCol, 'S10'),
                        (tranCol, '<i4'),
                        (datetimeCol, '<M8[us]'),
                        (dateCol, 'S12'),
                        (timeCol, 'S11'),
                        (latCol, '<f8'),
                        (lonCol, '<f8'),
                        (depObsCol, '<f8'),
                        (depInterpCol, '<f8'),
                        (videoCol, '<i4'),
                    ]
        return _base_dtype

    @property
    def veg_dtype(self):
        """" NumPy Data types for veg columns """
        _veg_dtype = []
        for veg in self.veg_columns:
            _veg_dtype.append((veg,'<i4'))
        return _veg_dtype

    @property
    def nparray(self):

        """" NumPy Array converted from the pandas dataframe """
        # Convert pandas data frame to numpy record array (a type of structured array)
        # NOTE 7/10/2023
        # before we can convert to nparray, need to force the column order.
        # this was not an issue in Arc10.6 but is in Arc10.8 and in ArcPro.
        # this step should occur in a separate function but struggling to implement, will refine at a later date.
        _base_col_order = [i[0] for i in self.base_dtype]
        _veg_col_order = [i[0] for i in self.veg_dtype]
        _col_order = _base_col_order + _veg_col_order

        _nparray = self.df[_col_order].to_records(index=False)

        # Change datatypes to ones compatible with ArcPy
        # This will also subset and re-order columns according to order in array
        dtypes = self.base_dtype + self.veg_dtype
        # Create NumPy data types object
        td_dtype = np.dtype(dtypes)

        # Apply data types to NumPy array
        _nparray = _nparray.astype(td_dtype)

        return _nparray

    def _process_data(self):
        """Run a series of functions to set or update values in the csv data"""
        # Rename columns to match final feature class
        self._rename_columns()
        # Add point ID column
        self._add_pointid()
        # Sort rows by transect id and timestamp
        self._sort_rows()
        # Fill Null records with a value
        self._fill_nulls()
        # Set site_code to lower case
        self._lower_site_code()
        # Create survey_id
        self._calc_survey_id()
        # Calculate nativesg column if at least one of the veg columns is a Native seagrass type
        if len(set(self.veg_columns).intersection(set(NATIVESG_CODES))) > 0:
            self.nativesg_columns = list(set(self.veg_columns).intersection(set(NATIVESG_CODES)))
            self._calc_nativesg()


    def _calc_nativesg(self):
        # Calculate the nativesg column based on values in the individual seagrass observations
        # Rows with any NULLs for native seagrass types
        df_null = (self.df[self.nativesg_columns] == NULL_VEG).any(axis=1)
        # Rows with any PRESENT for native seagrass types
        df_present = (self.df[self.nativesg_columns] == 1).any(axis=1)
        # Rows with all ABSENT for native seagrass types
        df_absent = (self.df[self.nativesg_columns] == 0).all(axis=1)

        # This Gets Maximum value of all veg columns - works when all input veg columns are non-null and video = 1
        # Not needed, but saved for reference
        # self.df[nativesgCol] = self.df[self.nativesg_columns].max(axis=1)

        # If _any_ native seagrass observations are present, set nativesg column to present (1)
        self.df.loc[df_present, nativesgCol] = 1
        # If _all_ native seagrass observations are absent, set nativesg column to absent (0)
        self.df.loc[df_absent, nativesgCol] = 0
        # If _any_ native seagrass observations are unknown, and there are NO present observations,
        #  set nativesg column to NULL/unknown (-9999)
        # Note:  in pandas ~ gets the inverse of a boolean dataframe
        self.df.loc[df_null & ~df_present, nativesgCol] = NULL_VEG
        # Finally, If video = 0 (bad quality) or null, set nativesg column to null/unknown (-9999)
        self.df.loc[self.df[videoCol] == 0, nativesgCol] = NULL_VEG
        self.df.loc[self.df[videoCol] == NULL_VIDEO, nativesgCol] = NULL_VEG

        # Add the native seagrass column to the veg columns list
        self.veg_columns.append(nativesgCol)

    @property
    def video_zero(self):
        # Validate video field.  Max of video column must be > 0
        if self.df[videoCol].max() < 1:
            return True
        else:
            return False

    @property
    def veg1_video0(self):
        # rows where video=0 but any veg column = 1
        df_video0 = self.df[videoCol] == 0
        df_veg1 = (self.df[self.source_veg_columns] == 1).any(axis=1)
        _veg1_video0 = self.df.loc[df_video0 & df_veg1, ptidCol].tolist()
        return _veg1_video0

    @property
    def dupe_ts(self):
        # Rows with duplicate time stamps
        _dupe_ts = self.df[self.df.duplicated([datetimeCol])][ptidCol].tolist()
        return _dupe_ts

    @property
    def null_veg(self):
        # Filter for Rows with null vegetation codes
        df_null = (self.df[self.source_veg_columns] == NULL_VEG).any(axis=1)
        # Make a list of row ids with null vegetation values
        _null_veg = self.df.loc[df_null][ptidCol].tolist()
        return _null_veg

    @property
    def null_video(self):
        # Filter for Rows with null video codes
        df_null = self.df[videoCol] == NULL_VIDEO
        # Make a list of row ids with null vegetation values
        _null_video = self.df.loc[df_null][ptidCol].tolist()
        return _null_video

    @property
    def veg_gt1(self):
        # Filter for rows with vegetation values > 1
        df_gt1 = (self.df[self.source_veg_columns] > 1).any(axis=1)
        # Make a list of row ids with null vegetation values
        _veg_gt1 = self.df.loc[df_gt1][ptidCol].tolist()
        return _veg_gt1

    @property
    def video_gt1(self):
        # Filter for Rows with video values > 1
        df_gt1 = self.df[videoCol] > 1
        # Make a list of row ids with null vegetation values
        _video_gt1 = self.df.loc[df_gt1][ptidCol].tolist()
        return _video_gt1

    @property
    def transect_video0(self):
        # Transects where max of video column < 1
        # df_max = self.df.groupby(tranCol)[videoCol].max()
        grouped = self.df.groupby(tranCol)
        df_max = grouped.filter(lambda x: x[videoCol].sum() == 0)

        # print df_max.loc[tranCol > 1]
        return df_max

    def _validate_data(self):
        # Check validation properties
        if self.video_zero or self.null_video or self.null_veg or self.video_gt1 \
                or self.veg_gt1 or self.veg1_video0 or self.dupe_ts:
            self.warnings = True

    def _rename_columns(self):
        # Rename some input columns to match feature class
        col_in2out = {
            sourceSiteCol: siteCol,
            sourceTrkCol: tranCol,
            sourceDateCol: dateCol,
            sourceTimeCol: timeCol,
            sourceDepthObsCol: depObsCol,
            sourceDepthInterpCol: depInterpCol,
        }
        self.df.rename(columns=col_in2out, inplace=True )

    def _add_pointid(self):
        # Add point ID column
        self.df[ptidCol] = self.df.index + 1

    def _fill_nulls(self):
        # Change null depth values to a placeholder integer value
        self.df.fillna({
            depObsCol: NULL_DEPTH,
            depInterpCol: NULL_DEPTH,
            videoCol: NULL_VIDEO,
        }, inplace=True)
        # Change veg nulls to placeholder integer value
        for col in self.veg_columns:
            self.df[col].fillna(NULL_VEG, inplace=True)

    def _lower_site_code(self):
        # Set site_code to lower case
        self.df[siteCol] = self.df[siteCol].str.lower()

    def _calc_survey_id(self):
        # Get the minimum date stamp for each transect and apply it back to the original data frame
        transect_mindate = self.df.groupby(tranCol)[datetimeCol].transform(lambda x: x.min())
        # Create survey_id
        # concatenate survey_id, date as string and transect number
        self.df[surveyidCol] = self.df[siteCol] + "_" + \
                          transect_mindate.dt.year.map(str) + \
                          transect_mindate.dt.month.map(lambda x: '{:02d}'.format(x)) + \
                          transect_mindate.dt.day.map(lambda x: '{:02d}'.format(x)) + "_" + \
                          self.df[tranCol].map(lambda x: '{:02d}'.format(x))

    def _sort_rows(self):
        # Sort the dataframe by transect number and time stamp
        self.df.sort_values([tranCol, datetimeCol], inplace=True)

class PointFC(object):
    """ Represents the a point feature class for a single site visit

    Properties:
    td_nparray -- The transect point data stored as a NumPy Array
    output_fc -- full path to the output feature class
    source_sr -- spatial reference of the source csv data
    output_sr -- spatial reference of the final output feature class

    """

    def __init__(self, td_nparray, output_fc):
        self.td_nparray = td_nparray # Transect data in a structured NumPy array
        self.output_fc = output_fc # Full path to output feature class
        self.source_sr = arcpy.SpatialReference(4326)  # Geographic, WGS-84
        self.output_sr = arcpy.SpatialReference(2927) # NAD_1983_HARN_StatePlane_Washington_South_FIPS_4602_Feet

    def create_fc(self):
        if arcpy.Exists(self.output_fc):
            arcpy.Delete_management(in_data=self.output_fc)
        # Create temporary lat/long version of point using "in_memory" location
        temp_fc = "in_memory/temp"
        arcpy.da.NumPyArrayToFeatureClass(self.td_nparray, temp_fc, [lonCol, latCol], self.source_sr)
        # Project to output coordinate system
        arcpy.Project_management(temp_fc, self.output_fc, self.output_sr)
        # Remove the temporary feature class
        arcpy.Delete_management(temp_fc)


class VegCodes(object):
    """ Represents the vegetation codes in the database

    Properties:
    vegcode_table -- geodatabase table with the vegetation codes
    veg_list -- list of vegetation codes

    """

    def __init__(self, vegcode_table):
         self.vegcode_table = vegcode_table

    @property
    def veg_list(self):
        """list of vegetation codes in the specified table"""
        # Convert veg_codes table to NumPy Array
        veg_table = arcpy.da.TableToNumPyArray(self.vegcode_table, (vegcodeCol), skip_nulls=True)
        # pull values out of veg_code column into a list
        _veg_list = [str(v) for v in veg_table[vegcodeCol].tolist()] # to deal with unicode issues
        del veg_table
        return _veg_list


class SiteVisit(object):

    def __init__(self, sitecode, yr, group):
        self.sitecode = sitecode
        self.yr = yr
        self.group = group
        self.fc = "_".join([sitecode, yr, group, 'transect','pt'])


class LogFile(object):
    """ Represents the a Log File used for error and validation reporting

    Properties:
    log_header -- Column headers for log file, varies according to log file type

    """

    def __init__(self, log_dir, log_type='csv2ptErrorLog'):
        self.log_dir = log_dir
        self.log_type = log_type
        self.log_file = os.path.join(self.log_dir, timeStamped(log_type))
        self.fh = None

    @property
    def log_header(self):
        if self.log_type == 'csv2ptErrorLog':
            return 'Directory,FileName,ErrorType,Details\n'
        if self.log_type == 'csv2ptWarningLog':
            return 'Directory,FileName,WarningType,Details\n'

    def open_log(self):
        # Open the log file for writing
        self.fh = open(self.log_file, 'w')
        # Write column headers to output file
        self.fh.write(self.log_header)

    def close_log(self):
        # Close the log file and set the file handle back to none
        self.fh.close()
        self.fh = None

    def write_csverr(self, csvsource):

        csv_dir = os.path.normpath(os.path.dirname(csvsource.file_path))
        csv_file = os.path.basename(csvsource.file_path)

        # Open the log file if it is not already open
        if not self.fh:
            self.open_log()

        if not csvsource.file_exists:
            err_type = "File Does Not Exist"
            details = ""
            self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")
        else:
            if csvsource.missing_columns:
                err_type = "Missing Columns"
                details = ';'.join(csvsource.missing_columns)
                self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")
                # File has no vegetation columns
            if csvsource.veg_columns == []:
                err_type = "No Vegetation Columns"
                details = ""
                self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")
                # File has latitude format errors
            if csvsource.lat_errors:
                err_type = "Bad Latitude Values"
                details = 'Rows: ' + ';'.join(str(r) for r in csvsource.lat_errors)
                self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")
                # file has longitude format errors
            if csvsource.lon_errors:
                err_type = "Bad Longitude Values"
                details = 'Rows: ' + ';'.join(str(r) for r in csvsource.lon_errors)
                self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")
            if csvsource.time_errors:
                err_type = "Bad Time Values"
                details = 'Rows: ' + ';'.join(str(r) for r in csvsource.time_errors)
                self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")
            if csvsource.date_errors:
                err_type = "Bad Date Values"
                details = 'Rows: ' + ';'.join(str(r) for r in csvsource.date_errors)
                self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")

    def write_direrr(self, csv_dir):

        # Open the log file if it is not already open
        if not self.fh:
            self.open_log()

        if not csv_dir.dir_exists:
            err_type = "Directory Does Not Exist"
            self.fh.write(",".join((csv_dir.csvdir, "", err_type, "")) + "\n")
        elif not csv_dir.files_exist:
            err_type = "No Transect Files Found"
            self.fh.write(",".join((csv_dir.csvdir, csv_dir.search_pattern, err_type, "")) + "\n")


    def write_datawarn(self, csvdata):

        csv_dir = os.path.normpath(os.path.dirname(csvdata.csv_source.file_path))
        csv_file = os.path.basename(csvdata.csv_source.file_path)

        # Open the log file if it is not already open
        if not self.fh:
            self.open_log()

        if csvdata.video_zero:
            err_type = "All Video < 1"
            details = ""
            self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")
        if csvdata.null_video:
            err_type = "Null Video Values"
            # print(', '.join(str(x) for x in list_of_ints))
            details = 'Rows: ' + ';'.join(str(i) for i in csvdata.null_video)
            self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")
        if csvdata.null_veg:
            err_type = "Null Veg Values"
            details = 'Rows: ' + ';'.join(str(i) for i in csvdata.null_veg)
            self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")
        if csvdata.video_gt1:
            err_type = "Video Values > 1"
            details = 'Rows: ' + ';'.join(str(i) for i in csvdata.video_gt1)
            self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")
        if csvdata.veg_gt1:
            err_type = "Veg Values > 1"
            details = 'Rows: ' + ';'.join(str(i) for i in csvdata.veg_gt1)
            self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")
        if csvdata.veg1_video0:
            err_type = "Veg = 1 and Video = 0"
            details = 'Rows: ' + ';'.join(str(i) for i in csvdata.veg1_video0)
            self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")
        if csvdata.dupe_ts:
            err_type = "Duplicate Time Stamps"
            details = 'Rows: ' + ';'.join(str(i) for i in csvdata.dupe_ts)
            self.fh.write(",".join((csv_dir, csv_file, err_type, details)) + "\n")


def timeStamped(fname, fmt='{fname}_%Y%m%d_%H%M%S.csv'):
    # Create time stamped filename
    return datetime.datetime.now().strftime(fmt).format(fname=fname)

def make_sitelist(sites_file):
    # Get all lines from input file without leading or trailing whitespace
    # and omit lines that are just whitespace
    site_list = [line.strip() for line in open(sites_file,'r') if not line.isspace()]
    return site_list

# General message accumulator
def msg(msg):
    arcpy.AddMessage(msg)


def main(in_dir, sites_file, vegcode_table, out_gdb, err_dir):
    # Main function to run code

    # Generate list of sites from text file
    site_codes = make_sitelist(sites_file)

    # Get list of vegetation codes available
    vegCodes = VegCodes(vegcode_table)

    # initiate Log File objects
    error_log = LogFile(err_dir,'csv2ptErrorLog')
    warning_log = LogFile(err_dir,'csv2ptWarningLog')

    # Loop through all of the sites in the site list
    for site in site_codes:
        msg("----- Processing site: {0} -----".format(site))

        # Locate and validate directory.  Get list of transect data files
        csvDir = CsvPath(site, in_dir)

        # Process valid directories
        if csvDir.valid:
            # # If there are matching transect files in the directory, process them
            # if csvDir.tdfiles:
            # Process all transect data files in the directory
            for tdfile in csvDir.tdfiles:
                # Site Visit object
                [sitecode, yr, group] = os.path.basename(tdfile).split('_')[0:3]
                site_visit = SiteVisit(sitecode, yr, group)
                # Create an csv source object
                csvSource = CsvSource(os.path.join(csvDir.csvdir,tdfile), vegCodes.veg_list)

                # If the source CSV file is valid, convert to point feature class
                if csvSource.valid:
                    fc_path = os.path.join(out_gdb, site_visit.fc)
                    transectData = CsvData(csvSource)
                    # Write validation warnings to log file
                    if transectData.warnings:
                        msg("Data Validation Warnings.\nWriting to log file: {1}".format(csvSource.file_path, warning_log.log_file))
                        warning_log.write_datawarn(transectData)
                    msg("Creating Point feature class {0}".format(fc_path))
                    ptFC = PointFC(transectData.nparray, fc_path)
                    ptFC.create_fc()
                else:
                    # Log invalid csv source files to error log
                    msg("Invalid source csv file {0}.\nWriting to error log file: {1}".format(csvSource.file_path, error_log.log_file))
                    error_log.write_csverr(csvSource)
        # Log Invalid directories to Error Log
        else:
            # write a line to the log file about the error
            msg("Directory, {0}, does not exist, "
                "or has no files matching the pattern: {1}_YYYY_##_TD.csv .\n"
                "Writing to error log file: {2}".format(csvDir.csvdir, csvDir.sitecode, error_log.log_file))
            error_log.write_direrr(csvDir)

    # Close the log files if they were opened
    if error_log.fh:
        error_log.close_log()
    if warning_log.fh:
        warning_log.close_log()


if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True, report=True)

    t0 = time.time()

    # Input parameter 1:  Parent directory for site data folders and input csv files
    in_dir = "C:/Users/lfer490/Documents/projects/svmp/2021_2022_site_processing/field_data/2021/site_folders/"

    # Input parameter 2:  Text file with list of sites to process
    sites_file = os.path.join(in_dir, "update.txt")

    # Input parameter 3: Table with vegetation codes
    vegcode_table = "C:/Users/lfer490/Documents/svmp_data/masterSVMPdb_20230502.gdb/veg_codes"

    # Input parameter 4: Output Geodatabase to store point feature classes
    out_gdb = "C:/Users/lfer490/Documents/svmp_data/transect_pt_data_2021.gdb"

    # Input parameter 5: Error Log directory
    # err_dir = in_dir
    err_dir = "C:/Users/lfer490/Documents/svmp_data/transect_pt_logs/"

    main(in_dir, sites_file, vegcode_table, out_gdb, err_dir)

    t1 = time.time()
    elapsed_time = str(t1-t0)
    print(f"Total time elapsed is: {elapsed_time}.")