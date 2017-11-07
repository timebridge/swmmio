#!/usr/bin/env python
#coding:utf-8
import re
import os
from time import ctime
import pandas as pd
from .utils import functions, spatial
import glob
import math
import geojson
from .utils import text as txt
from .utils.dataframes import create_dataframeINP, create_dataframeRPT, get_link_coords
from definitions import *

class Model(object):

	#Class representing a complete SWMM model incorporating its INP and RPT
	#files and data

	def __init__(self, in_file_path):

		#can init with a directory containing files, or the specific inp file
		"""
		initialize a swmmio.Model object by pointing it to a directory containing
		a single INP (and optionally an RPT file with matching filename) or by
		pointing it directly to an .inp file.
		"""

		inp_path = None
		if os.path.isdir(in_file_path):
			#a directory was passed in
			#print 'is dir = {}'.format(in_file_path)
			inps_in_dir = glob.glob1(in_file_path, "*.inp")
			if len(inps_in_dir) == 1:
				#there is only one INP in this directory -> good.
				inp_path = os.path.join(in_file_path, inps_in_dir[0])
				#print 'only 1 inp found = {}'.format(inp_path)

		elif os.path.splitext(in_file_path)[1] == '.inp':
			#an inp was passed in
			inp_path = in_file_path
			#print 'is inp path = {}'.format(in_file_path)

		if inp_path:
			wd = os.path.dirname(inp_path) #working dir
			name = os.path.splitext(os.path.basename(inp_path))[0] #basename
			self.name = name
			self.inp = inp(inp_path) #inp object
			self.rpt = None #until we can confirm it initializes properly
			#slots to hold processed data
			self.organized_node_data = None
			self.organized_conduit_data = None
			self.bbox = None #to remember how the model data was clipped
			self.scenario = self._get_scenario()

			#try to initialize a companion RPT object
			rpt_path = os.path.join(wd, name + '.rpt')
			if os.path.exists(rpt_path):
				try:
					self.rpt = rpt(rpt_path)
				except:
					print '{}.rpt failed to initialize'.format(name)

			self._nodes_df = None
			self._conduits_df = None
                        self._orifices_df = None
			self._subcatchments_df = None

	def rpt_is_valid(self , verbose=False):
		"""Return true if the .rpt file exists and has a revision date more
		recent than the .inp file. If the inp has an modified date later than
		the rpt, assume that the rpt should be regenerated"""

		if self.rpt is None:
			if verbose:
				print '{} does not have an rpt file'.format(self.name)
			return False


		#check if the rpt has ERRORS output from SWMM
		with open (self.rpt.path) as f:
			#jump to 500 bytes before the end of file
		    f.seek(self.rpt.file_size - 500)
		    for line in f:
		        spl = line.split()
		        if len(spl) > 0 and spl[0]=='ERROR':
					#return false at first "ERROR" occurence
		            return False

		rpt_mod_time = os.path.getmtime(self.rpt.path)
		inp_mod_time = os.path.getmtime(self.inp.path)

		if verbose:
			print "{}.rpt: modified {}".format(self.name, ctime(rpt_mod_time))
			print "{}.inp: modified {}".format(self.name, ctime(inp_mod_time))

		if inp_mod_time > rpt_mod_time:
			#inp datetime modified greater than rpt datetime modified
			return False
		else:
			return True

	def to_map(self, filename=None, inproj='epsg:2272'):

		conds = self.conduits()
		nodes = self.nodes()
		try:
			import pyproj
		except ImportError:
			raise ImportError('pyproj module needed. get this package here: https://pypi.python.org/pypi/pyproj')

		#SET UP THE TO AND FROM COORDINATE PROJECTION
		pa_plane = pyproj.Proj(init=inproj, preserve_units=True)
		wgs = pyproj.Proj(proj='longlat', datum='WGS84', ellps='WGS84') #google maps, etc

		#get center point
		c = ((nodes.X.max() + nodes.X.min())/2 , (nodes.Y.max() + nodes.Y.min())/2)
		c = pyproj.transform(pa_plane, wgs, c[0], c[1])
		bbox = [(nodes.X.min(), nodes.Y.min()),
				(nodes.X.max(), nodes.Y.max())]
		bbox = [pyproj.transform(pa_plane, wgs, *xy) for xy in bbox]


		geo_conduits = spatial.write_geojson(conds)
		geo_nodes = spatial.write_geojson(nodes, geomtype='point')

		if filename is None:
			filename = os.path.join(self.inp.dir, self.inp.name + '.html')

		with open(BETTER_BASEMAP_PATH, 'r') as bm:
			with open(filename, 'wb') as newmap:
				for line in bm:
					if '//INSERT GEOJSON HERE ~~~~~' in line:
						newmap.write('conduits = {};\n'.format(geojson.dumps(geo_conduits)))
						newmap.write('nodes = {};\n'.format(geojson.dumps(geo_nodes)))
						newmap.write('parcels = {};\n'.format(0))

					if 'center: [-75.148946, 39.921685],' in line:
						newmap.write('center:[{}, {}],\n'.format(c[0], c[1]))
					if '//INSERT BBOX HERE' in line:
						 newmap.write('map.fitBounds([[{}, {}], [{}, {}]]);\n'.format(bbox[0][0], bbox[0][1], bbox[1][0], bbox[1][1]))

					else:
						newmap.write(line)

	def _get_scenario(self):
		"""get a descrition of the model scenario by reading the raingage data"""
		rg = create_dataframeINP(self.inp.path, '[RAINGAGES]')
		storms = rg.DataSourceName.unique()
		if len(storms) > 1:
			return ', '.join(storms[:3]) + '...'
		else:
			return '{}'.format(storms[0])

	def conduits(self):

		"""
		collect all useful and available data related model conduits and
		organize in one dataframe.
		"""

		#check if this has been done already and return that data accordingly
		if self._conduits_df is not None:
			return self._conduits_df

		#parse out the main objects of this model
		inp = self.inp
		rpt = self.rpt

		#create dataframes of relevant sections from the INP
		conduits_df = create_dataframeINP(inp.path, "[CONDUITS]", comment_cols=False)
		xsections_df = create_dataframeINP(inp.path, "[XSECTIONS]", comment_cols=False)
		conduits_df = conduits_df.join(xsections_df)
		coords_df = create_dataframeINP(inp.path, "[COORDINATES]").drop_duplicates()

		if rpt:
			#create a dictionary holding data from an rpt file, if provided
			link_flow_df = create_dataframeRPT(rpt.path, "Link Flow Summary")
			conduits_df = conduits_df.join(link_flow_df)

		#add conduit coordinates
		#the xys.map() junk is to unpack a nested list
		verts = create_dataframeINP(inp.path, '[VERTICES]')
		xys = conduits_df.apply(lambda r: get_link_coords(r,coords_df,verts), axis=1)
		df = conduits_df.assign(coords=xys.map(lambda x: x[0]))

		#add conduit up/down inverts and calculate slope
		elevs = self.nodes()[['InvertElev']]
		df = pd.merge(df, elevs, left_on='InletNode', right_index=True, how='left')
		df = df.rename(index=str, columns={"InvertElev": "InletNodeInvert"})
		df = pd.merge(df, elevs, left_on='OutletNode', right_index=True, how='left')
		df = df.rename(index=str, columns={"InvertElev": "OutletNodeInvert"})
		df['UpstreamInvert'] = df.InletNodeInvert + df.InletOffset
		df['DownstreamInvert'] = df.OutletNodeInvert + df.OutletOffset
		df['SlopeFtPerFt'] = (df.UpstreamInvert - df.DownstreamInvert) / df.Length

		self._conduits_df = df

		return df

	def orifices(self):

		"""
		collect all useful and available data related model orifices and
		organize in one dataframe.
		"""

		#check if this has been done already and return that data accordingly
		if self._orifices_df is not None:
			return self._orifices_df

		#parse out the main objects of this model
		inp = self.inp
		rpt = self.rpt

		#create dataframes of relevant sections from the INP
		orifices_df = create_dataframeINP(inp.path, "[ORIFICES]", comment_cols=False)
		xsections_df = create_dataframeINP(inp.path, "[XSECTIONS]", comment_cols=False)
		orifices_df = conduits_df.join(xsections_df)
		coords_df = create_dataframeINP(inp.path, "[COORDINATES]").drop_duplicates()

		#add conduit coordinates
		#the xys.map() junk is to unpack a nested list
                verts = create_dataframeINP(inp.path, '[VERTICES]')
		xys = conduits_df.apply(lambda r: get_link_coords(r,coords_df,verts), axis=1)
		df = conduits_df.assign(coords=xys.map(lambda x: x[0]))

		self._conduits_df = df

		return df

	def nodes(self, bbox=None, subset=None):

		"""
		collect all useful and available data related model nodes and organize
		in one dataframe.
		"""

		#check if this has been done already and return that data accordingly
		if self._nodes_df is not None and bbox==self.bbox:
			return self._nodes_df

		#parse out the main objects of this model
		inp = self.inp
		rpt = self.rpt

		#create dataframes of relevant sections from the INP
		juncs_df = create_dataframeINP(inp.path, "[JUNCTIONS]")
		outfalls_df = create_dataframeINP(inp.path, "[OUTFALLS]")
		storage_df = create_dataframeINP(inp.path, "[STORAGE]")
		coords_df = create_dataframeINP(inp.path, "[COORDINATES]")

		#concatenate the DFs and keep only relevant cols
		all_nodes = pd.concat([juncs_df, outfalls_df, storage_df])
		cols =['InvertElev', 'MaxDepth', 'SurchargeDepth', 'PondedArea']
		all_nodes = all_nodes[cols]

		if rpt:
			#add results data if a rpt file was found
			depth_summ = create_dataframeRPT(rpt.path, "Node Depth Summary")
			flood_summ = create_dataframeRPT(rpt.path, "Node Flooding Summary")

			#join the rpt data (index on depth df, suffixes for common cols)
			rpt_df = depth_summ.join(flood_summ,lsuffix='_depth',rsuffix='_flood')
			all_nodes = all_nodes.join(rpt_df) #join to the all_nodes df

		all_nodes = all_nodes.join(coords_df[['X', 'Y']])
		def nodexy(row):
			if math.isnan(row.X) or math.isnan(row.Y):
				return None
			else:
				return [(row.X, row.Y)]

		xys = all_nodes.apply(lambda r: nodexy(r), axis=1)
		all_nodes = all_nodes.assign(coords = xys)

		self._nodes_df = all_nodes

		return all_nodes

	def subcatchments(self):
		"""
		collect all useful and available data related subcatchments and organize
		in one dataframe.
		"""
		subs = create_dataframeINP(self.inp.path, "[SUBCATCHMENTS]")
		subs = subs.drop([';', 'Comment', 'Origin'], axis=1)

		if self.rpt:
			flw = create_dataframeRPT(self.rpt.path, 'Subcatchment Runoff Summary')
			subs = subs.join(flw)

			#more accurate runoff calculations
			subs['RunoffAcFt'] = subs.TotalRunoffIn/ 12.0 * subs.Area
			subs['RunoffMGAccurate'] = subs.RunoffAcFt / 3.06888785

		self._subcatchments_df = subs

		return subs


	def node(self, node, conduit=None):

		"""
		DEPRECIATED/NOT SUPPORTED: organizeNodeData()

		method for provide information about specific model elements
		returns a node object given its ID"""
		if not self.organized_node_data:
			self.organized_node_data = self.organizeNodeData()

		n = self.organized_node_data['node_objects'][node]
		subcats_inp = self.inp.createDictionary("[SUBCATCHMENTS]")
		subcats_rpt = self.rpt.createDictionary('Subcatchment Runoff Summary')

		n.nodes_upstream = functions.trace_from_node(self, node, mode='up')['nodes']
		n.subcats_direct = [k for k,v in subcats_inp.items() if v[1]==node]
		n.subcats_upstream = [k for k,v in subcats_inp.items() if v[1] in n.nodes_upstream]


		n.drainage_area_direct = sum([float(x) for x in [v[2] for k,v in subcats_inp.items() if k in n.subcats_direct]])
		n.drainage_area_upstream = sum([float(x) for x in [v[2] for k,v in subcats_inp.items() if k in n.subcats_upstream]])

		n.runoff_upstream_mg = sum([float(x) for x in [v[5] for k,v
								in subcats_rpt.items() if k in n.subcats_upstream]])
		n.runoff_upstream_cf = n.runoff_upstream_mg*1000000/7.48
		return n


	def export_to_shapefile(self, shpdir, prj=None):
		"""
		export the model data into a shapefile. element_type dictates which type
		of data will be included.

		default projection is PA State Plane - untested on other cases
		"""

		#CREATE THE CONDUIT shp
		conds = self.conduits()
		conds_path = os.path.join(shpdir, self.inp.name + '_conduits.shp')
		spatial.write_shapefile(conds, conds_path, prj=prj)

		#CREATE THE NODE shp
		nodes = self.nodes()
		nodes_path = os.path.join(shpdir, self.inp.name + '_nodes.shp')
		spatial.write_shapefile(nodes, nodes_path, geomtype='point', prj=prj)


class SWMMIOFile(object):

	defaultSection = "Link Flow Summary"

	def __init__(self, file_path):

		#file name and path variables
		self.path = file_path
		self.name = os.path.splitext(os.path.basename(file_path))[0]
		self.dir = os.path.dirname(file_path)
		self.file_size = os.path.getsize(file_path)


	def findByteRangeOfSection(self, startStr):

		#returns the start and end "byte" location of substrings in a text file

		with open(self.path) as f:
			start = None
			end = None
			l = 0 #line bytes index
			for line in f:

				#if start and len(line) <= 3 and (l - start) > 100:
				if start and line.strip() == "" and (l - start) > 100:
					#LOGIC ^ if start exists (was found) and the current line length is 3 or
					#less (length of /n ) and we're more than 100 bytes from the start location
					#then we are at the first "blank" line after our start section (aka the end of the section)
					end = l
					break

				if (startStr in line) and (not start):
					start = l

				l += len(line) + len("\n") #increment length (bytes?) of current position

		return [start, end]

	def createDictionary (self, sectionTitle = defaultSection):

		"""
		Help info about this method.
		"""

		#preppedTempFilePath = self.readSectionAndCleanHeaders(sectionTitle) #pull relevant section and clean headers
		preppedTempFilePath = txt.extract_section_from_file(self.path, sectionTitle)
		if not preppedTempFilePath:
			return None #if nothing was found, do nothing

		passedHeaders = False

		with open(preppedTempFilePath) as file:
			the_dict = {}
			for line in file:

				if len(line) <=3 and not ";" in line: break
				if not passedHeaders:
					passedHeaders = True
					continue

				#check if line is commented out (having a semicolon before anything else) and skip accordingly
				if ";" == line.replace(" ", "")[0]:
					continue #omit this entire line

				line = line.split(";")[0] #don't look at anything to right of a semicolon (aka a comment)

				line = ' '.join(re.findall('\"[^\"]*\"|\S+', line))
				rowdata = line.replace("\n", "").split(" ")
				the_dict[str(rowdata[0])] = rowdata[1:] #create dictionary row with key and array of remaing stuff on line as the value

		os.remove(preppedTempFilePath)

		return the_dict

class rpt(SWMMIOFile):

	#creates an accessible SWMM .rpt object, inherits from SWMMIO object
	defaultImageDir = r"P:\Tools\Pipe Capacity Graphics\Scripts\image"
	def __init__(self, filePath):

		SWMMIOFile.__init__(self, filePath) #run the superclass init

		with open (filePath) as f:
			for line in f:
				if "Starting Date" in line:
					simulationStart = line.split(".. ")[1].replace("\n", "")
				if "Ending Date" in line:
					simulationEnd = line.split(".. ")[1].replace("\n", "")
				if "Report Time Step ........." in line:
					timeStepMin = int(line.split(":")[1].replace("\n", ""))
					break

		self.simulationStart = simulationStart
		self.simulationEnd = simulationEnd
		self.timeStepMin = timeStepMin

		#grab the date of analysis
		with open (filePath) as f:
			f.seek(self.file_size - 500) #jump to 500 bytes before the end of file
			for line in f:
				if "Analysis begun on" in line:
					date = line.split("Analysis begun on:  ")[1].replace("\n", "")

		self.dateOfAnalysis = date

		#assign the header list
		#self.headerList = swmm_headers.rptHeaderList
		self.byteLocDict = None #populated if necessary elsewhere (LEGACY, can prob remove)
		self.elementByteLocations = {"Link Results":{}, "Node Results":{}} #populated if necessary elsewhere

	def createByteLocDict (self, sectionTitle = "Link Results"):

		#method creates a dictionary with Key = to SWMM element ID and
		#Value as the starting byte location of its time series in the rpt file
		#for rapidly accessing large rpt files

		#create set of other headers that are not the desired one, use to find end of section
		possibleNextSections = set(['Link Results', 'Node Results', 'Subcatchment Results']) - set([sectionTitle])

		print possibleNextSections

		startByte = self.findByteRangeOfSection(sectionTitle)[0] #+ len('\n  ************') #move past the first asterisks

		id_byteDict = {}
		with open(self.path) as f:

			f.seek(startByte) #jump to general area of file if we know it
			l = startByte
			for line in f:

				#if "<<<" in line and ">>>" in line:
				if "<<<" and ">>>" in line:	#cr
					#found the begining of a link's section
					lineCleaned = ' '.join(re.findall('\"[^\"]*\"|\S+', line))
					rowdata = lineCleaned.replace("\n", "").split(" ")

					#add to the dict
					id_byteDict.update({rowdata[2]:l})

				if any(header in line for header in possibleNextSections):
					#checks if line includes any of the other headers,
					#if so, we found next section, stop building dict
					break

				l += len(line) + len("\n") #increment length (bytes) of current position

		self.byteLocDict = id_byteDict
		self.elementByteLocations.update({sectionTitle:id_byteDict})
		return id_byteDict

	def returnDataAtDTime(self, id, dtime, sectionTitle="Link Results", startByte=0):

		#this is a slow ass function, when the file is big - can we improve this?
		byteLocDict = self.elementByteLocations[sectionTitle]
		if byteLocDict:
			startByte = byteLocDict[id]

		elif startByte == 0:
			startByte = self.findByteRangeOfSection(sectionTitle)[0]
			print 'startByte ' + str(startByte)

		with open(self.path) as f:

			f.seek(startByte) #jump to general area of file if we know it
			subsectionFound = False

			for line in f:
				if id in line: subsectionFound = True

				if subsectionFound and dtime in line:
					line = ' '.join(re.findall('\"[^\"]*\"|\S+', line))
					rowdata = line.replace("\n", "").split(" ")
					return rowdata

class inp(SWMMIOFile):

	#creates an accessible SWMM .inp object
	#make sure INP has been saved in the GUI before using this

	def __init__(self, filePath):
		#is this class necessary anymore?
		SWMMIOFile.__init__(self, filePath) #run the superclass init



#end
