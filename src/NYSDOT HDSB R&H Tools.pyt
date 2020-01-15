from collections import namedtuple
from datetime import datetime
import math
import os
import re
import time

import arcpy


class QueryMustSelectTwoRoutes(Exception):
    pass


class NoGapsOnInputRoutes(Exception):
    pass


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (.pyt file is name of Toolbox)"""
        # Metadata for ArcGIS
        self.label = 'NYSDOT R&H Tools'
        self.alias = 'nysdotrh'

        # List of the tools available in this Toolbox based on class names
        self.tools = [
            GenerateRouteIdAndDotId, CheckCalibration, GenerateIntervalPoints,
            FindPrimaryAndReverseRouteDivergence, GapLength,
        ]


class ReverseDescriptionDeventFields(object):
    def __init__(self):
        self.label = 'Reverse Begin and End Descriptions'
        self.description = (
            'This tool is to be executed after using the Reverse Route tool from the Roads and Highways Toolbar. ' +
            'You must either "Save Edits" or "Apply Updates" before using this tool!' +
            'It takes flips the Begin and End descriptions in the feature class to correspond with the edit.'
        )
        self.category = 'Editing'
        self.canRunInBackground = False

    def getParameterInfo(self):

        return


class GenerateRouteIdAndDotId(object):
    def __init__(self):
        self.label = 'Generate New DOT_ID and/or ROUTE_ID'
        self.description = (
            'Get the next DOT_ID and/or ROUTE_ID from the database sequence.'
        )
        self.category = 'Database'
        self.canRunInBackground = False

    def getParameterInfo(self):
        sde_filepath_param = arcpy.Parameter(
            displayName='Database Connection File (SDE file)',
            name='sde_filepath',
            datatype='DEWorkspace',
            parameterType='Required',
            direction='Input'
        )

        route_id_flag_param = arcpy.Parameter(
            displayName='Get new ROUTE_ID',
            name='route_id_flag',
            datatype='GPBoolean',
            parameterType='Optional',
            direction='Input'
        )

        dot_id_flag_param = arcpy.Parameter(
            displayName='Get a new DOT_ID',
            name='dot_id_flag',
            datatype='GPBoolean',
            parameterType='Optional',
            direction='Input'
        )

        route_id_flag_param.value = 'true'
        dot_id_flag_param.value = 'true'

        if arcpy.Exists(r'Database Connections\dev_elrs_ad_Lockroot.sde'):
            sde_filepath_param.value = r'Database Connections\dev_elrs_ad_Lockroot.sde'
        if arcpy.Exists(r'Database Connections\dev_elrs_ad_lockroot.sde'):
            sde_filepath_param.value = r'Database Connections\dev_elrs_ad_Lockroot.sde'

        return [
            sde_filepath_param, route_id_flag_param, dot_id_flag_param
        ]

    def execute(self, parameters, messages):
        sde_filepath = parameters[0].valueAsText
        route_id_flag = parameters[1].valueAsText
        dot_id_flag = parameters[2].valueAsText

        if route_id_flag == 'true':
            route_id_flag = True
        else:
            route_id_flag = False

        if dot_id_flag == 'true':
            dot_id_flag = True
        else:
            dot_id_flag = False

        route_id_sql = "SELECT NEXT VALUE FOR ELRS.S_LRSN_MILEPOINT_ROUTE_ID AS ROUTE_ID"
        dot_id_sql = "SELECT NEXT VALUE FOR ELRS.S_RIS_ROADWAY_REGISTRY AS DOT_ID"

        connection = arcpy.ArcSDESQLExecute(sde_filepath)

        if route_id_flag:
            route_id = connection.execute(route_id_sql)
            messages.AddMessage ("###############################")
            messages.AddMessage ("ROUTE_ID: {:.0f}".format(route_id))
            messages.AddMessage ("###############################")

        if dot_id_flag:
            dot_id = connection.execute(dot_id_sql)
            messages.AddMessage ("###############################")
            messages.AddMessage ("DOT_ID: {:.0f}".format(dot_id))
            messages.AddMessage ("###############################")

        return True


class CheckCalibration(object):
    def __init__(self):
        """Define the Tool's class (Tool name is name of class)"""
        self.label = 'Compare Calibrated Length and GIS Length'
        self.description = (
            'Calculate max M value and GIS length and compare')
        self.category = 'Calibration'
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define the parameters definitions of the Tool. These are the inputs
        required to run the tool. These definitions set up the user input screen
        visible in ArcGIS."""
        param0 = arcpy.Parameter(
            displayName='Input Routes',
            name='milepoint_path',
            datatype='GPFeatureLayer',
            parameterType='Required',
            direction='Input')
        param0.filter.list = ['Line']

        param1 = arcpy.Parameter(
            displayName='Output Routes',
            name='output_path',
            datatype='GPFeatureLayer',
            parameterType='Required',
            direction='Output')

        param2 = arcpy.Parameter(
            displayName='Output Length Field Name',
            name='length_field',
            datatype='Field',
            parameterType='Required',
            direction='Input')

        param3 = arcpy.Parameter(
            displayName='Output Measure Field Name',
            name='m_field',
            datatype='Field',
            parameterType='Required',
            direction='Input')

        param4 = arcpy.Parameter(
            displayName='Output Difference Field Name',
            name='diff_field',
            datatype='Field',
            parameterType='Required',
            direction='Input')

        param5 = arcpy.Parameter(
            displayName='Output Part Count Field Name',
            name='part_count_field',
            datatype='Field',
            parameterType='Required',
            direction='Input')

        param6 = arcpy.Parameter(
            displayName='Decimal Places of Output Difference Field',
            name='decimal_places',
            datatype='GPLong',
            parameterType='Required',
            direction='Input')

        param7 = arcpy.Parameter(
            displayName='Use selected features',
            name='selected_features_flag',
            datatype='GPBoolean',
            parameterType='Required',
            direction='Input')

        # Set sensible default values
        param2.value = 'GIS_length'
        param3.value = 'M_length'
        param4.value = 'length_diff'
        param5.value = 'part_count'
        param6.value = 3
        param7.value = 'true'

        params = [param0, param1, param2, param3, param4, param5, param6, param7]

        return params

    def isLicensed(self):
        return True

    def execute(self, parameters, messages):
        """Retrieve parameters as Python variables and execute tool"""
        # Retrieve parameters and input them to the GIS work method
        input_path = parameters[0].valueAsText
        output_path = parameters[1].valueAsText
        length_field = parameters[2].valueAsText
        m_field = parameters[3].valueAsText
        diff_field = parameters[4].valueAsText
        part_count_field = parameters[5].valueAsText
        decimal_places = parameters[6].value
        selected_features_flag = parameters[7].value

        self.check_calibration(input_path, output_path, length_field,
                               m_field, diff_field, part_count_field,
                               selected_features_flag,
                               decimal_places=decimal_places)
        return True

    def check_calibration(self, input_path, output_path,
                          length_field, m_field, diff_field,
                          part_count_field, selected_features_flag,
                          decimal_places=3):
        """Manage input/output and execute GIS work"""
        # Set up string of date/time. This will be appended to temporary file
        # names. If the tool fails to execute, the in-memory data from the
        # failed execution will still exist. Appending this to filenames ensures
        # the tool will not error out if re-run in the same ArcGIS session
        time_chars = str(datetime.now())
        time_chars = time_chars.replace(' ', '').replace('-', '').replace(
            ':', '').replace('.', '').strip('0')

        # Apply selected features flag and create in_memory layer of input
        selected_count = int(arcpy.GetCount_management(input_path)[0])
        if (not selected_features_flag or
           (selected_features_flag and selected_count == 0)):
            arcpy.SelectLayerByAttribute_management(
                input_path, 'CLEAR_SELECTION')
        input_features = 'in_memory/in_layer{}'.format(time_chars)
        arcpy.CopyFeatures_management(input_path, input_features)
        arcpy.MakeFeatureLayer_management(input_features, 'input_layer')

        # Add relevant fields to in_memory layer
        self.shape_length('input_layer', length_field)
        self.max_m('input_layer', m_field)
        self.calc_diff('input_layer', m_field, length_field,
                       diff_field, decimal_places=decimal_places)
        self.part_count('input_layer', part_count_field)

        # Create copy of output on disk
        arcpy.CopyFeatures_management('input_layer', output_path)

        return True


    def shape_length(self, input_path, new_field):
        """Calculate GIS length and add it to table. Input must be meters
        (requires projected geographic reference system in meters e.g. UTM),
        output will be miles."""
        # Generate new field if required
        field_names = [field.name for field in arcpy.ListFields(input_path)]
        if not new_field in field_names:
            arcpy.AddField_management(input_path, new_field, 'DOUBLE')
        # Calculate GIS length for each feature. If multipart feature, use
        # Euclidean distance to account for length between "gaps"
        with arcpy.da.UpdateCursor(input_path, ['SHAPE@', new_field]) as cursor:
            for row in cursor:
                shape, length = row
                if shape.partCount > 1:
                    length = self.length_between_points(shape)
                elif shape.partCount == 1:
                    length = shape.length * 0.000621371
                else:
                    raise ValueError('0 Parts in feature. No length calculated.')
                cursor.updateRow([shape, length])
        return True

    def length_between_points(self, shape, conversion_factor=0.000621371):
        """Calculate GIS length of shape, including Euclidean distance
        between parts. Input should be UTM projection, output will be miles
        unless conversion factor is changed."""
        # Initialize variables to store values
        length = 0.
        point_list = []

        # Get the ArcPy point objects of each point along the feature
        part_count = shape.partCount
        for i in range(part_count):
            array = shape.getPart(i)
            for point in array:
                point_list.append(point)

        # Calculate the distance between each set of points. If x or y
        # coordinates are identical, use subtraction. Otherwise, use the
        # Pythagorean Theorem
        point_count = len(point_list)
        for i, ii in zip(range(point_count), range(point_count)[1:]):
            p1, p2 = point_list[i], point_list[ii]
            if (p1.X - p2.X) == 0.:
                segment_l = abs(p1.Y - p2.Y) * conversion_factor
            elif (p1.Y - p2.Y) == 0.:
                segment_l = abs(p1.X - p2.X) * conversion_factor
            else:
                a = abs(float(p1.X - p2.X))
                b = abs(float(p1.Y - p2.Y))
                c = math.sqrt(a**2 + b**2)
                segment_l = c * conversion_factor
            length += segment_l

        return length

    def max_m(self, input_path, new_field):
        """Calculate the maximum calibrated measure value and append the table"""
        # Add field if needed
        field_names = [field.name for field in arcpy.ListFields(input_path)]
        if not new_field in field_names:
            arcpy.AddField_management(input_path, new_field, 'DOUBLE')

        # Get M value from the extent of the feature
        with arcpy.da.UpdateCursor(input_path, ['SHAPE@', new_field]) as cursor:
            for row in cursor:
                shape, m_field = row
                m_field = shape.extent.MMax
                cursor.updateRow([shape, m_field])

        return True


    def calc_diff(self, input_path, m_field, l_field, new_field, decimal_places=3):
        """Subract the GIS length from the max M and add to table"""
        # Add new field if needed
        field_names = [field.name for field in arcpy.ListFields(input_path)]
        if not new_field in field_names:
            arcpy.AddField_management(input_path, new_field, 'DOUBLE')

        # Retrieve m/l values and subtract
        fields = [m_field, l_field, new_field]
        with arcpy.da.UpdateCursor(input_path, fields) as cursor:
            for row in cursor:
                m, length = row[0], row[1]
                try:
                    diff = round(float(m) - float(length), decimal_places)
                except:
                    diff = None
                cursor.updateRow([m, length, diff])

        return True

    def part_count(self, input_path, new_field):
        """Add part count to table"""
        # Add new field if needed
        field_names = [field.name for field in arcpy.ListFields(input_path)]
        if not new_field in field_names:
            arcpy.AddField_management(input_path, new_field, 'SHORT')

        # Add part count to table
        fields = [new_field, 'SHAPE@']
        with arcpy.da.UpdateCursor(input_path, fields) as cursor:
            for row in cursor:
                part_count, shape = row
                part_count = shape.partCount
                cursor.updateRow([part_count, shape])

        return True


class GenerateIntervalPoints(object):
    def __init__(self):
        """Define the Tool's class (Tool name is name of class)"""
        self.label = 'Generate Equal Interval Points'
        self.description= (
            'Generate equal points at equal interval of calibrated and GIS Length')
        self.category = 'Calibration'
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define the parameters definitions of the Tool. These are the inputs
        required to run the tool. These definitions set up the user input screen
        visible in ArcGIS."""
        param0 = arcpy.Parameter(
            displayName='Input Routes',
            name='input_path',
            datatype='GPFeatureLayer',
            parameterType='Required',
            direction='Input')
        param0.filter.list = ['Line']

        param1 = arcpy.Parameter(
            displayName='Route ID Field Name',
            name='route_id',
            datatype='Field',
            parameterType='Required',
            direction='Input')

        param2 = arcpy.Parameter(
            displayName='Length interval between points (mi)',
            name='interval',
            datatype='GPString',
            parameterType='Required',
            direction='Input')

        param3 = arcpy.Parameter(
            displayName='Output M Points',
            name='m_points',
            datatype='GPFeatureLayer',
            parameterType='Required',
            direction='Output')

        param4 = arcpy.Parameter(
            displayName='Output Length Points',
            name='l_points',
            datatype='GPFeatureLayer',
            parameterType='Required',
            direction='Output')

        param5 = arcpy.Parameter(
            displayName='Use selected features',
            name='selected_features_flag',
            datatype='GPBoolean',
            parameterType='Required',
            direction='Input')

        # Set logical defaults
        param1.parameterDependencies = [param0.name]
        param1.value = 'ROUTE_ID'
        param2.value = '0.10'
        param5.value = 'false'

        # Set categories
        # param3.category = 'Outputs'
        # param4.category = 'Outputs'

        params = [param0, param1, param2, param3, param4, param5]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return True

    def execute(self, parameters, messages):
        """Execute the tool. Called on run from ArcMap session"""
        # Retrieve and manipulate parameters as Python variables
        input_path = parameters[0].valueAsText
        route_id = parameters[1].valueAsText
        interval = float(parameters[2].valueAsText)
        m_path = parameters[3].valueAsText
        l_path = parameters[4].valueAsText
        selected_features_flag = parameters[5].value

        # Set up string of date/time. This will be appended to temporary file
        # names. If the tool fails to execute, the in-memory data from the
        # failed execution will still exist. Appending this to filenames ensures
        # the tool will not error out if re-run in the same ArcGIS session
        time_chars = str(datetime.now())
        time_chars = time_chars.replace(' ', '').replace('-', '').replace(
            ':', '').replace('.', '').strip('0')

        # Apply selected features flag
        if not selected_features_flag:
            arcpy.SelectLayerByAttribute_management(
                input_path, 'CLEAR_SELECTION')
        input_features = 'in_memory/in_layer{}'.format(time_chars)
        arcpy.CopyFeatures_management(input_path, input_features)
        arcpy.MakeFeatureLayer_management(input_features, 'input_layer')

        # Do the GIS Work
        self.generate_points('input_layer', route_id, interval,
                             m_path, l_path, time_chars)

##        self.add_data_to_layer(m_path)
##        self.add_data_to_layer(l_path)

    def add_data_to_layer(self, filepath):
        """Add output data to active ArcMap dataframe"""
        arcpy.MakeFeatureLayer_management(filepath, 'view_layer')
        mxd = arcpy.mapping.MapDocument('CURRENT')
        df = arcpy.mapping.ListDataFrames(mxd, 'Layers')[0]
        arcpy.mapping.AddLayer(df, 'view_layer', 'TOP')
        arcpy.RefreshActiveView()
        arcpy.RefreshTOC()
        return


    def generate_points(self, input_layer, route_id, interval,
                        m_path, l_path, time_chars):
        """Generate points along the feature at an equal interval of
        GIS length and calibrated length"""
        # Retrieve the M values from the extent
        m_values = []
        with arcpy.da.SearchCursor(input_layer, ['SHAPE@', route_id]) as cursor:
            for row in cursor:
                shape, id_ = row
                m_values.append([int(id_), float(shape.extent.MMax)])

        # Create a table of data including route_id and measure value for each
        # interval of the input feature
        table_data= []
        for row in m_values:
            id_, max_m = row
            placement_m = interval
            if max_m < interval:
                continue
            else:
                id_intervals = [[id_, placement_m]]
            while placement_m < max_m:
                placement_m += interval
                id_intervals.append([id_, round(placement_m, 3)])
            table_data.append(id_intervals)

        # Create an in_memory table of the route_id/measure value
        table_name = 'table_{}'.format(time_chars)
        table_path = 'in_memory/{}'.format(table_name)
        arcpy.CreateTable_management('in_memory', table_name)
        arcpy.AddField_management(table_path, route_id, 'LONG')
        arcpy.AddField_management(table_path, 'm_value', 'DOUBLE')
        with arcpy.da.InsertCursor(table_path, [route_id, 'm_value']) as cursor:
            for id_intervals in table_data:
                for row in id_intervals:
                    cursor.insertRow(row)

        # Create an event layer of the m value table along the input features
        # and save it to disk as the m_path
        arcpy.MakeRouteEventLayer_lr(
            input_layer, route_id, table_path,
            '{} POINT m_value'.format(route_id), 'point_layer',
            add_error_field='ERROR_FIELD')
        arcpy.CopyFeatures_management('point_layer', m_path)

        # Create a point at each equal distance along the GIS length and save it
        # to disk as the l_path
        # arcpy.GeneratePointsAlongLines_management(
        #     input_layer, l_path,'DISTANCE',
        #     Distance='{} miles'.format(interval))
        describe = arcpy.Describe(m_path)
        spatial_info = namedtuple('spatial_info', 'spatialReference extent')
        sp_info = spatial_info(spatialReference=describe.spatialReference,
                               extent=describe.extent)

        input_fc = 'in_memory/layer_{}'.format(time_chars)
        arcpy.CopyFeatures_management(input_layer, input_fc)
        self.create_points_from_lines(input_fc, l_path,
                                      sp_info.spatialReference, route_id,
                                      dist=interval)

        # add a field for the length value to l_path. If the GIS ID is the same
        # as the previous iteration, add the interval to the value and write
        # to the table. If the GIS ID is different, start back at the interval
        # length. Assumes the features are ordered in ascending GIS length
        arcpy.AddField_management(l_path, 'l_value', 'DOUBLE')
        with arcpy.da.UpdateCursor(l_path, [route_id, 'l_value']) as cursor:
            last_id = 'none'
            l_value = interval
            for row in cursor:
                this_id = row[0]
                if this_id == last_id:
                    l_value += interval
                else:
                    l_value = interval
                last_id = this_id
                cursor.updateRow([this_id, l_value])
        return True

    def create_points_from_lines(self, input_fc, output_fc, spatial_ref,
                                 route_id, percent=False, dist=True,
                                 add_end_points=False):
        """Convert line features to feature class of points. This function
        was taken from the GeneratePointsAlongLines_management script in
        ArcGIS 10.5.1 and reworked into this toolbox so that this toolbox
        has 10.3 support

        :param input_fc: Input line features
        :param output_fc: Output point feature class
        :param spatial_ref: The spatial reference of the input
        :param percent: If creating points by percentage (or distance)
        :param dist: The distance used to create points (if percentage == False).
        The distance should be in units of the input (see convert_units)
        :param add_end_points: If True, an extra point will be added from start
        and end point of each line feature
        :return: None
        """
        if percent:
            is_percentage = True
        else:
            is_percentage = False

        ## Hardcoded conversion from miles to meters
        if dist:
            dist *= 1609.347218694438

        # Create output feature class
        arcpy.CreateFeatureclass_management(
            os.path.dirname(output_fc),
            os.path.basename(output_fc),
            geometry_type="POINT",
            spatial_reference=spatial_ref)

        # Add a field to transfer FID from input
        fid_name = 'ORIG_FID'
        arcpy.AddField_management(output_fc, fid_name, 'LONG')

        # Create new points based on input lines
        in_fields = ['SHAPE@', 'OID@']
        out_fields = ['SHAPE@', fid_name]

        with arcpy.da.SearchCursor(input_fc, in_fields) as search_cursor:
            with arcpy.da.InsertCursor(output_fc, out_fields) as insert_cursor:
                for row in search_cursor:
                    line = row[0]

                    if line:  # if null geometry--skip
                        if line.type == 'polygon':
                            line = line.boundary()

                        if add_end_points:
                            insert_cursor.insertRow([line.firstPoint, row[1]])

                        increment = (percent or dist)
                        cur_length = increment

                        if is_percentage:
                            max_position = 1.0
                        else:
                            max_position = line.length

                        while cur_length < max_position:
                            new_point = line.positionAlongLine(cur_length,
                                                               is_percentage)
                            insert_cursor.insertRow([new_point, row[1]])
                            cur_length += increment

                        if add_end_points:
                            end_point = line.positionAlongLine(1, True)
                            insert_cursor.insertRow([end_point, row[1]])

            try:
                arcpy.JoinField_management(output_fc,
                                           fid_name,
                                           input_fc,
                                           arcpy.Describe(input_fc).OIDFieldName,
                                           fields=[route_id])
            except arcpy.ExecuteError:
                # In unlikely event that JoinField fails, proceed regardless,
                # as spatial and join field are already complete
                pass

        return


class GapLength(object):
    def __init__(self):
        self.label = 'Calculate Euclidean Distance of Gapped Routes'
        self.description = 'Finds Euclidean distance of gapped routes'
        self.category = 'Signed Routes'
        self.canRunInBackground = False

    def getParameterInfo(self):
        param0 = arcpy.Parameter(
            displayName='Input Routes',
            name='input_path',
            datatype='GPFeatureLayer',
            parameterType='Required',
            direction='Input')
        param0.filter.list = ['Line']

        param1 = arcpy.Parameter(
            displayName='Route ID',
            name='rid_field',
            datatype='Field',
            parameterType='Required',
            direction='Input')

        param2 = arcpy.Parameter(
            displayName='Output Table',
            name='output_table',
            datatype='DETable',
            parameterType='Required',
            direction='Output')

        param1.value = 'ROUTE_ID'

        params = [param0, param1, param2]

        return params

    def isLicensed(self):
        return True

    def execute(self, parameters, messages):
        input_path = parameters[0].valueAsText
        rid_field = parameters[1].valueAsText
        output_table = parameters[2].valueAsText

        arcpy.CreateTable_management(os.path.dirname(output_table),
                                     os.path.basename(output_table))
        arcpy.AddField_management(output_table, rid_field, 'TEXT')
        arcpy.AddField_management(output_table, 'GAP_ID', 'TEXT')
        arcpy.AddField_management(output_table, 'GAP_LENGTH', 'DOUBLE')

        print('\nFinding route gaps and their lengths in:\n {}'.format(input_path))
        arcpy.AddMessage('\nFinding route gaps and their lengths in:\n {}'.format(input_path))
        output_list = []
        with arcpy.da.SearchCursor(input_path, [rid_field, 'SHAPE@']) as cursor:
            for row in cursor:
                rid, polyline = row
                part_count = polyline.partCount
                if part_count > 1:
                    points = []
                    for i in range(part_count):
                        points.append(
                            (polyline.getPart(i)[0],
                             polyline.getPart(i)[-1]))
                    gap_count = 1
                    for first_part, last_part in zip(points, points[1:]):
                        gap_start = first_part[1]
                        gap_end = last_part[0]
                        length = self.gap_dist(gap_start, gap_end)
                        if length == 0.:
                            continue
                        output_list.append(
                            [rid,
                             '{0}:{1}'.format(gap_count, gap_count + 1),
                             length
                            ]
                        )
                        gap_count += 1

        gap_count = len(output_list)
        if gap_count == 0:
            raise NoGapsOnInputRoutes(
                'No gaps were found on the input routes.')

        print('\nWriting to output table:\n {}'.format(output_table))
        arcpy.AddMessage('\nWriting to output table:\n {}'.format(output_table))
        insert_fields = [rid_field, 'GAP_ID', 'GAP_LENGTH']
        with arcpy.da.InsertCursor(output_table, insert_fields) as out_cursor:
            for output in output_list:
                out_cursor.insertRow(output)

        return True

    def gap_dist(self, p1, p2, conversion_factor=0.000621371):
        if (p1.X - p2.X) == 0.:
            segment_l = abs(p1.Y - p2.Y) * conversion_factor
        elif (p1.Y - p2.Y) == 0.:
            segment_l = abs(p1.X - p2.X) * conversion_factor
        else:
            a = abs(float(p1.X - p2.X))
            b = abs(float(p1.Y - p2.Y))
            c = math.sqrt(a**2 + b**2)
            segment_l = c * conversion_factor
        return segment_l


class FindPrimaryAndReverseRouteDivergence(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = 'Identify Primary and Reverse Route Divergence'
        self.description = (
            'This tool will return features in the output feature class ' +
            'where the primary and reverse of a route split onto separate ' +
            'carriageways or an empty geometry if the primary/reverse ' +
            'take the same path.'
        )
        self.category = 'Signed Routes'
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        in_routes = arcpy.Parameter(
            displayName='Input Routes',
            name='milepoint_path',
            datatype='GPFeatureLayer',
            parameterType='Required',
            direction='Input')
        in_routes.filter.list = ['Line']

        dot_id_field_name = arcpy.Parameter(
            displayName='DOT ID Field Name',
            name='dot_id_field',
            datatype='GPString',
            parameterType='Required',
            direction='Input')

        in_dot_id = arcpy.Parameter(
            displayName='DOT ID (6 digits)',
            name='dot_id',
            datatype='GPString',
            parameterType='Required',
            direction='Input')

        county_order_field_name = arcpy.Parameter(
            displayName='County Order Field Name',
            name='county_order_field',
            datatype='GPString',
            parameterType='Required',
            direction='Input')

        in_county_order = arcpy.Parameter(
            displayName='County Order (2 digits)',
            name='county_order',
            datatype='GPString',
            parameterType='Required',
            direction='Input')

        direction_field_name = arcpy.Parameter(
            displayName='Direction Field Name',
            name='direction_field',
            datatype='GPString',
            parameterType='Required',
            direction='Input')

        out_fc_path = arcpy.Parameter(
            displayName='Output Feature Class',
            name='output_path',
            datatype='GPFeatureLayer',
            parameterType='Required',
            direction='Output')

        # dot_id_field_name.value = 'DOT_ID'
        # county_order_field_name.value = 'COUNTY_ORDER'
        # direction_field_name.value = 'DIRECTION'

        dot_id_field_name.filter.type = 'ValueList'
        county_order_field_name.filter.type = 'ValueList'
        direction_field_name.filter.type = 'ValueList'
        # dot_id_field_name.filter.list = ['DOT_ID']

        # params = [
        #     in_routes, dot_id_field_name, in_dot_id,
        #     county_order_field_name, in_county_order,
        #     direction_field_name, out_fc_path,
        # ]

        params = [
            in_routes, in_dot_id, in_county_order,
            dot_id_field_name, county_order_field_name, direction_field_name,
            out_fc_path,
        ]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        field_indices = [3, 4, 5]
        if parameters[0].value:
            field_names = [
                field.name for field in arcpy.ListFields(parameters[0].valueAsText)
                if field.name not in ['OBJECTID', 'SHAPE', 'SHAPE.LEN']
            ]

            for idx in field_indices:
                parameters[idx].filter.list = sorted(field_names)
            # parameters[1].filter.list = sorted(field_names)
            # parameters[3].filter.list = sorted(field_names)
            # parameters[5].filter.list = sorted(field_names)

            default_values = ['DOT_ID', 'COUNTY_ORDER', 'DIRECTION']

            for idx, default_name in zip(field_indices, default_values):
                if not parameters[idx].altered:
                    if default_name in field_names:
                        parameters[idx].value = default_name
        else:
            parameters[1].filter.list = []
            parameters[3].filter.list = []
            parameters[5].filter.list = []
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        arcpy.env.workspace = 'in_memory'

        input_routes = parameters[0].valueAsText
        dot_id = parameters[1].valueAsText
        county_order = parameters[2].valueAsText
        dot_id_field = parameters[3].valueAsText
        county_order_field = parameters[4].valueAsText
        direction_field = parameters[5].valueAsText
        output_path = parameters[6].valueAsText

        # Progressor step 0
        where_clause = "({0} = '{1}') AND ({2} = '{3}')".format(
            dot_id_field, dot_id, county_order_field, county_order)
        arcpy.AddMessage('Selecting routes where: "{}"'.format(where_clause))

        route_layer = arcpy.MakeFeatureLayer_management(
            input_routes, 'route_layer')
        route_layer = arcpy.SelectLayerByAttribute_management(
            route_layer, 'NEW_SELECTION', where_clause)

        selection_count = int(arcpy.GetCount_management(route_layer)[0])
        if selection_count > 2:
            raise QueryMustSelectTwoRoutes(
                'The input parameters led to more than 2 routes being selected')
        elif selection_count < 2:
            raise QueryMustSelectTwoRoutes(
                'The input parameters led to less than 2 routes being selected. ' +
                'Did you clear your selection on the input routes?')

        # Progressor step 1
        arcpy.AddMessage('\nSplitting feature class in memory')
        for layer in arcpy.ListFeatureClasses():
            arcpy.Delete_management(os.path.join('in_memory', layer))
        arcpy.SplitByAttributes_analysis (
            route_layer, 'in_memory', [direction_field])


        feature_classes = sorted(
            [os.path.join('in_memory', feature_class)
             for feature_class in arcpy.ListFeatureClasses()],
             reverse=True)
        arcpy.AddMessage('in mem: {}'.format(feature_classes))

        # Progressor step 2
        arcpy.AddMessage('\nGenerating output feature class')
        erase_path = os.path.join('in_memory', 'erased')
        arcpy.Erase_analysis(
            os.path.join(arcpy.env.workspace, feature_classes[0]),
            os.path.join(arcpy.env.workspace, feature_classes[1]),
            erase_path)

        # Progressor step 3
        arcpy.AddMessage('\nSplitting multipart features on output')
        arcpy.MultipartToSinglepart_management(
            erase_path, output_path)

        arcpy.AddMessage('\nSuccessfully saved output to:\n{}\n'.format(output_path))

        output_count = int(arcpy.GetCount_management(output_path)[0])
        if output_count == 0:
            arcpy.AddWarning(
                '000117: Warning empty output generated.\n')

        return
