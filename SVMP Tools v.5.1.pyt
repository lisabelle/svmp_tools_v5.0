import arcpy
import os
import svmpUtils as utils
import sys
# Not working, yet.....
# tool_path = os.path.dirname(os.path.realpath(__file__))
# script_path = os.path.join(tool_path, "scripts")
# sys.path.append(script_path)


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "SVMP Tools v5.0"
        self.alias = "svmp50"

        import svmpUtils as utils

        # List of tool classes associated with this toolbox
        self.tools = [TransectDatatoPtFC, TransectAndSiteStatistics]


class TransectDatatoPtFC(object):

    def __init__(self):
        """Tool to convert SVMP video survey point files (csv) to point feature classes"""
        self.label = "(1) Transect Data to Point Feature Class"
        self.description = "This tool converts SVMP video survey files from csv format to point feature classes"
        self.canRunInBackground = True

    def getParameterInfo(self):
        """Define parameter definitions"""

        # Input parameter 1:  Parent directory for site data folders and input csv files
        in_dir = arcpy.Parameter(
            displayName="Input Data Parent Directory",
            name="in_dir",
            datatype="Folder",
            parameterType="Required",
            direction="Input"
        )
        # Input parameter 2:  Text file with list of sites to process
        sites_file = arcpy.Parameter(
            displayName="List of Sites file",
            name="sites_file",
            datatype="File",
            parameterType="Required",
            direction="Input"
        )
        # Input parameter 3: Table with vegetation codes
        vegcode_table = arcpy.Parameter(
            displayName="Vegetation Code Table",
            name="vegcode_table",
            datatype="Table",
            parameterType="Required",
            direction="Input"
        )

        out_gdb = arcpy.Parameter(
            displayName="Output Geodatabase",
            name="out_gdb",
            datatype="Workspace",
            parameterType="Required",
            direction="Input"
        )
        out_gdb.filter.list = ['Local Database','Remote Database']

        err_dir = arcpy.Parameter(
            displayName="Output Error Log Directory",
            name="err_dir",
            datatype="Folder",
            parameterType="Required",
            direction="Input"
        )

        # Default values  -- Change or remove as needed
        # in_dir.value = 'C:/Users/lfer490/Documents/projects/svmp/2021_2022_site_processing/field_data/2022/site_folders/'
        # sites_file.value = os.path.join(in_dir.value, "sites2process_all.txt")
        # vegcode_table.value = "C:/Users/lfer490/Documents/svmp_data/masterSVMPdb_20230502.gdb/veg_codes"
        # out_gdb.value = "C:/Users/lfer490/Documents/svmp_data/transect_pt_data_2022.gdb"
        # err_dir.value = "C:/Users/lfer490/Documents/svmp_data/transect_pt_logs/"

        params = [in_dir, sites_file, vegcode_table, out_gdb, err_dir]

        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""

        if parameters[2].value:
            vegcode_table = parameters[2].value
            vegcode_field = 'veg_code'
            # table = arcpy.Describe(vegcode_path).baseName
            field_names = [f.name for f in arcpy.ListFields(vegcode_table)]
            if vegcode_field not in field_names:
                errtext = f"[SVMP ERROR]: The selected table, {vegcode_table}, has no field {vegcode_field}."
                errtext += "\nChoose a different table."
                parameters[2].setErrorMessage(errtext)
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        import csv2pt

        # Input parameter 1:  Parent directory for site data folders and input csv files
        in_dir = parameters[0].valueAsText # 'C:/Users/lfer490/Documents/projects/svmp/2021_2022_site_processing/field_data/2022/site_folders/'

        # Input parameter 2:  Text file with list of sites to process
        sites_file = parameters[1].valueAsText # 'C:/Users/lfer490/Documents/projects/svmp/2021_2022_site_processing/field_data/2022/site_folders/update.txt'

        # Input parameter 3: Table with vegetation codes
        vegcode_table = parameters[2].valueAsText # "C:/Users/lfer490/Documents/svmp_data/masterSVMPdb_20230502.gdb/veg_codes"

        # Input parameter 4: Output Geodatabase to store point feature classes
        out_gdb = parameters[3].valueAsText # "C:/Users/lfer490/Documents/svmp_data/transect_pt_data_2022.gdb"

        # Input parameter 5: Error Log directory
        # err_dir = parameters[4].valueAsText  # in_dir
        err_dir = parameters[4].valueAsText # "C:/Users/lfer490/Documents/svmp_data/transect_pt_logs/"

        # Call the main function to process the csv point data
        csv2pt.main(in_dir, sites_file, vegcode_table, out_gdb, err_dir)

        return


class TransectAndSiteStatistics(object):

    def __init__(self):
        """Tool to calculate SVMP transect and site statistics from transect point features """
        self.label = "(2) Calculate Transect and Site Statistics"
        self.description = "This tool calculates SVMP transect and site statistics from transect point features"
        self.canRunInBackground = True

        # Input sources for parameter lists derived from database tables:
        # index = index of parameter (as returned from getParameterInfo)
        # table = SVMP table containing the column used for parameter list
        # field = column in the table that includes the items for the parameter list
        # reverse = boolean indicate reverse sort option
        self.svmpgdb_idx = 1  # parameter index for master SVMP geodatabase
        self.parameter_inputs = {
            "survey year": {
                "index": 3,
                "table": utils.sitevisitsTbl,
                "field": utils.visityearCol,
                "reverse": True,
            },
            "veg type": {
                "index": 4,
                "table": utils.vegcodesTbl,
                "field": utils.vegcodeCol,
                "reverse": False,
            },
            "study code": {
                "index": 6,
                "table": utils.studyassociationsTbl,
                "field": utils.studycodeCol,
                "reverse": False,
            },
            "sample selection": {
                "index": 7,
                "table": utils.sitesamplesTbl,
                "field": utils.sampselCol,
                "reverse": False,
            }
        }


    def getParameterInfo(self):
        """Define parameter definitions"""
        # Input parameter 1:  Geodatabase with Transect Point Feature Class(es)
        transect_gdb = arcpy.Parameter(
            displayName="Transect Point Geodatabase",
            name="transect_gdb",
            datatype="Workspace",
            parameterType="Required",
            direction="Input"
        )
        transect_gdb.filter.list = ['Local Database','Remote Database']

        # Input parameter 2:  SVMP Geodatabase with Tables needed for selecting correct transects
        svmp_gdb = arcpy.Parameter(
            displayName="SVMP Core Geodatabase",
            name="svmp_gdb",
            datatype="Workspace",
            parameterType="Required",
            direction="Input"
        )
        svmp_gdb.filter.list = ['Local Database','Remote Database']

        # Input parameter 3: Site Statistics Geodatabase with Template results tables
        stats_gdb = arcpy.Parameter(
            displayName="Site Statistics Database",
            name="stats_db",
            datatype="Workspace",
            parameterType="Required",
            direction="Input"
        )
        stats_gdb.filter.list = ['Local Database','Remote Database']

        # Input parameter 4: Survey Year to be Processed
        survey_year = arcpy.Parameter(
            displayName="Survey Year",
            name="survey_year",
            datatype="String",
            parameterType="Required",
            direction="Input"
        )
        survey_year.enabled = False  # Disabled until value in svmp_gdb

        # Input parameter 5: Vegetation Type to be Processed
        veg_code = arcpy.Parameter(
            displayName="Vegetation Type",
            name="veg_code",
            datatype="String",
            parameterType="Required",
            direction="Input"
        )
        veg_code.enabled = False # Disabled until value in svmp_gdb

        # Input parameter 6: Optional List of Sites file
        sites_file = arcpy.Parameter(
            displayName = "List of Sites File",
            name = "sites_file",
            datatype="File",
            parameterType="Optional",
            direction="Input",
        )

        # Input parameter 7: Study or Studies to Be Processed
        study = arcpy.Parameter(
            displayName="Study",
            name="study",
            datatype="String",
            parameterType="Optional",
            direction="Input",
            multiValue=True,
            category="Optional - Choose Study"
        )
        study.filter.type = "ValueList"
        study.enabled = False # Disabled until value in svmp_gdb

        # Input parameter 8: Vegetation Type to be Processed
        samp_sel = arcpy.Parameter(
            displayName="Sample Selection Method",
            name="samp_sel",
            datatype="String",
            parameterType="Optional",
            direction="Input",
            multiValue=True,
            category="Optional - Sample Selection Method"
        )
        samp_sel.filter.type = "ValueList"
        samp_sel.enabled = False # Disabled until value in svmp_gdb

        # Default values  -- Change or remove as needed
        # transect_gdb.value = "C:/Users/lfer490/Documents/svmp_data/transect_pt_data_2022.gdb"
        # svmp_gdb.value = "C:/Users/lfer490/Documents/svmp_data/masterSVMPdb_20230502.gdb"
        # stats_gdb.value = "C:/Users/lfer490/Documents/svmp_data/svmp_sitesdb.gdb"
        # sites_file.value = 'C:/Users/lfer490/Documents/projects/svmp/2021_2022_site_processing/field_data/2022/site_folders/update.txt'
        # survey_year.value = 2022
        # veg_code.value = "veg"

        params = [transect_gdb, svmp_gdb, stats_gdb, survey_year, veg_code, sites_file, study, samp_sel]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        # Run this section if the SVMP GDB parameter is not blank
        if parameters[self.svmpgdb_idx].altered:
            svmp_gdb = str(parameters[self.svmpgdb_idx].value)  # Full pathname of SVMP geodatabase
            svmp_gdb_base = os.path.basename(svmp_gdb)  # Geodatabase name without path
            # Loop through all parameters that use SVMP GDB to set filters
            for param, input in self.parameter_inputs.items():
                index = input["index"]  # parameter index
                table = input["table"]   # Table name for the parameter filter table
                table_path = os.path.normpath(os.path.join(svmp_gdb,table))  # Full path to survey year table
                if arcpy.Exists(table_path):
                    field = input["field"] # Field used to get filter list
                    if utils.fieldExists(table_path, field):
                        if parameters[index].value and "[SVMP ERROR]" in parameters[index].valueAsText:
                            parameters[index].value = "" # Reset parameter value if there was previous error
                        values_list = utils.unique_values(table_path, field) # List of unique values for parameter filter list
                        parameters[index].filter.list = sorted(values_list, reverse=input["reverse"])
                        parameters[index].enabled = True  # Enable the parameter
                    else:
                        field_error = "[SVMP ERROR]: Field {0} is not present in Table {1} in {2}".format(field, table, svmp_gdb_base)
                        parameters[index].value = field_error
                        parameters[index].enabled = False
                        parameters[index].filter.list = []
                        # self.parameter_inputs["survey year"]["error"] = field_error
                else:
                    table_error = "[SVMP ERROR]: Table {0} is not present in {1}".format(table, svmp_gdb_base)
                    parameters[index].value = table_error
                    parameters[index].enabled = False
                    parameters[index].filter.list = []
                    # This works in updateParamters, but the value for "error" in the parameters_input dictionary does
                    #   not get propagated to updateMessages
                    # self.parameter_inputs["survey year"]["error"] = table_error
                    # parameters[self.parameter_inputs["veg type"]["index"]].value = self.parameter_inputs["survey year"]["error"]
        # If SVMP GDB parameter is blank
        else:
            for param, input in self.parameter_inputs.items():
                index = input["index"]
                parameters[index].value = ""
                parameters[index].filter.list = []
                parameters[index].enabled = False

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""

        # Loop through all parameters that require tables from the SVMP GDB
        for param, input in self.parameter_inputs.items():
            index = input["index"]
            if parameters[index].value:
                # If there is a missing table or field
                if "[SVMP ERROR]" in parameters[index].valueAsText:
                    if not parameters[self.svmpgdb_idx].hasError():
                        # Set an error on the SVMP geodatabase
                        parameters[self.svmpgdb_idx].setErrorMessage("Database is missing required tables or fields.  Select a new GDB")

        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        import statsdb

        # Input parameter 1:  Geodatabase with individual transect point data -- REQUIRED
        transect_gdb = parameters[0].valueAsText  #

        # Input parameter 2:  SVMP Geodatabase with Base Tables -- REQUIRED
        svmp_gdb = parameters[1].valueAsText

        # Input parameter 3: Site Statistics Geodatabase with Template results tables -- REQUIRED
        stats_gdb = parameters[2].valueAsText

        # Input parameter 4: Survey Year to be Processed -- REQUIRED
        survey_year = parameters[3].valueAsText

        # Input parameter 5: Vegetation Type to be Processed -- REQUIRED
        veg_code = parameters[4].valueAsText

        # Input parameter 6: List of Sites file -- OPTIONAL
        sites_file = parameters[5].valueAsText

        # Input parameter 7: Study or Studies to Be Processed -- OPTIONAL
        study = parameters[6].valueAsText

        # Input parameter 8: Vegetation Type to be Processed -- OPTIONAL
        samp_sel = parameters[7].valueAsText

        # Call the main function to process the csv point data
        statsdb.main(transect_gdb, svmp_gdb, stats_gdb, survey_year, veg_code, sites_file, study, samp_sel)

        return
