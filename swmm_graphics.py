#!/usr/bin/env python
#coding:utf-8

#graphical functions for SWMM files
#import swmmio
import swmm_utils as su
import draw_utils as du
import parcels
from time import gmtime, strftime
import re
import os
import numpy
#from PIL import Image, ImageDraw, ImageFont
import svgwrite
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF

from images2gif import writeGif
import glob
import shutil
import math
import datetime
from datetime import timedelta
import pickle



def saveImage(img, model, imgName=None, imgDir=None, antialias=True, open=True, fileExt=".png", verbose=False):

	#if imgName not specified, define as model name
	if not imgName:
		imgName = model.inp.name

	#create the saving location as necessary
	if not imgDir:
		inp = model.inp
		rpt = model.rpt
		standardDir = os.path.join(inp.dir, "img")
		if not os.path.exists(standardDir):
			#if no directory is specified and none exists at
			#standard location create a new directory
			os.makedirs(os.path.join(inp.dir, "img"))
		newFile = os.path.join(standardDir, imgName) + fileExt
	else:
		#imDir is specified by user
		if not os.path.exists(imgDir):
			#if directory doesn't exist, create new
			os.makedirs(imgDir)
		newFile = os.path.join(imgDir, imgName) + fileExt

	if verbose: print "saving image to: " + newFile

	#shrink for antialiasing niceness (though this blows up the file about 5x)
	if fileExt != '.png':

		img.saveas(newFile)
		drawing = svg2rlg(newFile)
		pdf = os.path.join(standardDir, imgName)+".pdf"
		renderPDF.drawToFile(drawing, pdf)
		if open:
			#os.startfile(newFile)
			os.startfile(pdf)
	else:
		if antialias:
			imgSize = (img.getbbox()[2], img.getbbox()[3])
			size = (int(imgSize[0]*0.5), int(imgSize[1]*0.5))
			img.thumbnail(size, Image.ANTIALIAS)

		img.save(newFile)
		if open:
			os.startfile(newFile)



def createFeaturesDict(options={},  bbox=None, shiftRatio=None, width=1024):

	#create dictionary containing draw coordinate data for the given basemap options

	gdb = options['gdb']
	features = options['features']

	import arcpy
	if arcpy.Exists(gdb):
		#if the gdb exists, then we can move forward
		#for feature, data in features.iteritems():
		for featureOps in features:
			feature = featureOps['feature']
			if arcpy.Exists(os.path.join(gdb, feature)):

				#this check so that we don't repeat this heavy op over an over
				featureDict = su.shape2Pixels(feature, bbox=bbox, targetImgW=width, where=None, cols=featureOps['cols'],
														shiftRatio=shiftRatio, gdb=gdb)

				featureOps.update({'featureDict':featureDict}) #retain this for later if necessary
				#featureOps['featureDict'].update{featureDict}) #retain this for later if necessary
			else:
				print '{} not found'.format(feature)
	else:
		print '{} not found'.format(gdb)
		return None

	return features


def drawParcels(draw, parcel_flooding_results,  options={}, bbox=None, width=1024, shiftRatio=None, svg=True):


	parcels_pixels = su.shape2Pixels("PWD_PARCELS_SPhila_ModelSpace", where=None,
									cols = ["PARCELID", "SHAPE@"], targetImgW=width,
									shiftRatio=shiftRatio, bbox=bbox)

	#if 'parcels' in parcel_flooding_results:
	#	parcel_flooding_results = parcel_flooding_results['parcels'] #HACKKK
	threshold = options['threshold']
	newflood = moreflood = floodEliminated= floodlowered= 0
	for PARCELID, parcel in parcel_flooding_results.iteritems():
		fill = du.lightgrey #default
		try:
			if parcel.is_delta:
				#we're dealing with "compare" dictionary
				if parcel.delta_type == 'increased_flooding':
					#parcel previously flooded, now floods more
					fill = du.red

				if parcel.delta_type == 'new_flooding':
					#parcel previously did not flood, now floods in proposed conditions
					fill = du.purple

				if parcel.delta_type == 'decreased_flooding':
					#parcel flooding problem decreased
					fill =du.lightblue #du.lightgrey
					floodlowered += 1

				if parcel.delta_type == 'eliminated_flooding':
					#parcel flooding problem eliminated
					fill =du.lightgreen


			elif parcel.flood_duration > threshold:
				fill = du.col2RedGradient(parcel.flood_duration + 0.5, 0, 3)

			parcel_pix = parcels_pixels['geometryDicts'][PARCELID]
			if not svg:
				draw.polygon(parcel_pix['draw_coordinates'], fill=fill, outline=options['outline'])
			else:
				color = 'rgb' + str(fill)

				if options['outline']:
					outline = 'rgb' + str(options['outline'])
				else:
					outline = None

				polygon = draw.polygon(parcel_pix['draw_coordinates'], fill=color)
				draw.add(polygon)

		except:
			pass

	#print "newflood = {}\nmoreflood={}\nfloodEliminated={}\nfloodlowered={}".format(newflood, moreflood, floodEliminated, floodlowered)
	#saveImage(img, None, imgName=model.inp.name + "_parcels", imgDir=r'C:\Users\Adam.Erispaha\Desktop\S Phila SWMM\img')

def drawBasemap(draw, img=None, featureDicts=None, options={}, bbox=None, shiftRatio=None, width=1024, xplier=1, svg=True):



	if not featureDicts:
		#generate the dict containing drawable data
		featureDicts = createFeaturesDict(options,  bbox=bbox, shiftRatio=shiftRatio, width=width) #will return None if probs finding data

	#if featureDicts:
	#	features = featureDicts['features']

	polyDrawCount = 0
	anno_streets = []
	#for feature, data in features.iteritems():
	for feature in featureDicts:
		featureDict = feature['featureDict']
		for poly, polyData in featureDict['geometryDicts'].iteritems():
			if polyData['geomType'] == 'polygon':
				if not svg:
					draw.polygon(polyData['draw_coordinates'], fill=feature['fill'], outline=feature['outline'])
				else:
					color = 'rgb' + str(feature['fill'])
					outline = 'rgb' + str(feature['outline'])
					polygon = draw.polygon(polyData['draw_coordinates'], fill=color)
					draw.add(polygon)
				polyDrawCount += 1

			elif polyData['geomType'] == 'polyline':

				if not svg:
					draw.line(polyData['draw_coordinates'], fill=feature['fill'], width=feature['width']*xplier)
					if 'ST_NAME' in polyData:
						su.annotateLine(img, polyData, annoKey='ST_NAME', labeled = anno_streets)
						polyDrawCount += 1
				else:
					color = 'rgb' + str(feature['fill'])

					polyline = draw.polyline(polyData['draw_coordinates'], stroke=color, opacity=0)
					draw.add(polyline)
					if 'ST_NAME' in polyData:
						du.annotateLine(draw, polyData, annoKey='ST_NAME', labeled = anno_streets, svg=True)
						polyDrawCount += 1


def drawModel (model, **kwargs):

	#unpack and update the options
	ops = du.default_draw_options()
	for key, value in kwargs.iteritems():
		ops.update({key:value})
	#return ops
	width = ops['width']
	xplier = ops['xplier']
	bbox = ops['bbox']
	imgName = ops['imgName'] # for some reason saveImage() won't take the dict reference
	imgDir = ops['imgDir'] # for some reason saveImage() won't take the dict reference
	svg = ops['svg']

	focusConduits = []
	for node in ops['traceUpNodes']:
		#return list of elements upstream of node
		focusConduits += su.traceFromNode(model, node, mode='up')['conduits']

	for node in ops['traceDnNodes']:
		#return list of elements downstream of node
		focusConduits += su.traceFromNode(model, node, mode='down')['conduits']

	#antialias X2
	xplier *= width/1024 #scale the symbology sizes
	width = width*2

	#parse out the main objects of this model
	inp = model.inp
	rpt = model.rpt

	#organize relavant data from SWMM files, #FIX TO NOT RUN ALL THIS IF NOT NECESSARY
	conduitData = model.organizeConduitData(bbox) #diction	ary of overall model data, dimension, and the conduit dicts
	conduits = conduitData['conduit_objects']
	pixelData = su.convertCoordinatesToPixels(conduits, targetImgW=width, bbox=bbox)
	shiftRatio = pixelData['shiftRatio']
	imgSize = pixelData['imgSize']

	if not bbox:
		bbox = pixelData['boundingBox'] #this is the box tightly wrapping all conduits, will be <= the user provided bbox

	#node data
	nodeData = model.organizeNodeData(bbox)
	nodes = nodeData['node_objects']
	su.convertCoordinatesToPixels(nodes, targetImgW=width, bbox=bbox)


	if svg:
		draw = svgwrite.Drawing(size=imgSize)
	else:
		img = Image.new('RGB', imgSize, ops['bg'])
		draw = ImageDraw.Draw(img)

	anno_results = {}

	#DRAW THE PARCELS
	if ops['parcelSymb']:

		#grab the parcel&node dict that associates node IDs to each parcel
		parcel_flooding_results = parcels.parcel_flood_duration(model, parcel_features=ops['parcelSymb']['feature'],
															threshold=ops['parcelSymb']['threshold'],
															gdb= ops['parcelSymb']['gdb'], bbox=None,
															anno_results=anno_results)

		drawParcels(draw, parcel_flooding_results['parcels'], options=ops['parcelSymb'], bbox=bbox, width=width, shiftRatio=shiftRatio, svg=svg)

	#DRAW THE BASEMAP
	if ops['basemap']:
		if not svg:
			drawBasemap(draw, img=img, options=ops['basemap'], width=width, bbox=bbox, shiftRatio=shiftRatio, xplier = xplier, svg=svg)
		else:
			drawBasemap(draw, options=ops['basemap'], width=width, bbox=bbox, shiftRatio=shiftRatio, xplier = xplier, svg=svg)

	drawCount = 0

	#DRAW THE CONDUITS
	if ops['conduitSymb']:
		for id, conduit in conduits.iteritems():

			if conduit.coordinates: #has coordinate? draw. protect from rdii junk
				su.drawConduit(conduit, draw, ops['conduitSymb'], rpt=rpt, xplier = xplier, highlighted=focusConduits, svg=svg)
				drawCount += 1

	#DRAW THE NODES
	if ops['nodeSymb']:
		for id, node in nodes.iteritems():
			if node.coordinates: #this prevents draws if no flow is supplied (RDII and such)
				#drawNode(node, draw, options, rpt=None, dTime=None, xplier=1):
				su.drawNode(node, draw, rpt=rpt, options=ops['nodeSymb'], dTime=None, xplier=xplier, svg=svg)
				drawCount += 1


	#if conduitDrawCount > 0 and conduitDrawCount % 2000 == 0: print str(conduitDrawCount) + " pipes processed "
	#su.annotateMap (draw, model, options=ops, results=anno_results)

	if svg:
		del ops
		#SAVE IMAGE TO DISK
		saveImage(draw, model, imgName, fileExt='.svg')
	else:
		del draw, ops
		saveImage(img, model, imgName, imgDir=imgDir)


def animateModel(model, startDtime=None, endDtime=None, **kwargs):


	#unpack and update the options
	ops = du.default_draw_options()
	for key, value in kwargs.iteritems():
		ops.update({key:value})
	#return ops
	width = ops['width']
	xplier = ops['xplier']
	bbox = ops['bbox']
	imgName = ops['imgName'] # for some reason saveImage() won't take the dict reference

	#for antialiasing, double size for now
	xplier *= width/1024 #scale the symbology sizes
	width = width*2

	#parse out the main object of this model
	inp = model.inp
	rpt = model.rpt

	#organize relavant data from SWMM files
	conduitData = model.organizeConduitData(bbox) #dictionary of overall model data, dimension, and the conduit dicts
	conduitDicts = conduitData['conduitDictionaries']
	pixelData = su.convertCoordinatesToPixels(conduitDicts, targetImgW=width, bbox=bbox)
	shiftRatio = pixelData['shiftRatio']
	imgSize = pixelData['imgSize']
	#node data
	nodeData = model.organizeNodeData(bbox)
	nodeDicts = nodeData['nodeDictionaries']
	su.convertCoordinatesToPixels(nodeDicts, targetImgW=width, bbox=bbox)

	#grab start and end simulation time if no range provided -> this will be huge if you're not careful!
	if not startDtime: startDtime = rpt.simulationStart
	if not endDtime: endDtime = rpt.simulationEnd

	#make sure dates are valid (within range)
	simStartDT = datetime.datetime.strptime(rpt.simulationStart, "%b-%d-%Y %H:%M:%S") #lower bound of time
	simEndDT = datetime.datetime.strptime(rpt.simulationEnd, "%b-%d-%Y %H:%M:%S") #upper bound of time
	if startDtime and endDtime:
		userStartDT = 	datetime.datetime.strptime(startDtime, "%b-%d-%Y %H:%M:%S")
		userEndDT = 	datetime.datetime.strptime(endDtime, "%b-%d-%Y %H:%M:%S")
		timeStepMod = userStartDT.minute % rpt.timeStepMin
		if userStartDT < simStartDT or userEndDT > simEndDT or timeStepMod != 0 or userEndDT < userStartDT:
			#user has entered fault date times either by not being within the
			#availble data in the rpt or by starting at something that doesn't fit the timestep
			print "PROBLEM WITH DATETIME ENTERED. Make sure it fits within data and start time rest on factor of timestep in minutes."
			print "userStartDT = ", userStartDT, "\nuserEndDT = ", userEndDT, "\nsimStartDT = ", simStartDT, "\nsimEndDT = ", simEndDT, "\nTIMESTEP = ", rpt.timeStepMin
			return None

	currentT = datetime.datetime.strptime(startDtime, "%b-%d-%Y %H:%M:%S") #SWMM dtime format needed
	endDtime = datetime.datetime.strptime(endDtime, "%b-%d-%Y %H:%M:%S") #SWMM dtime format needed
	delta = timedelta(minutes=rpt.timeStepMin) #SWMM reporting time step (or step to animate with)


	#use or create working dir
	if not os.path.exists( os.path.join(inp.dir, "img") ):
		os.makedirs(os.path.join(inp.dir, "img"))
	imgDir = os.path.join(inp.dir, "img")

	byteLocDictionaryFName = os.path.join(imgDir, inp.name) + "_rpt_key.txt"
	if not os.path.isfile(byteLocDictionaryFName):

		#this is a heavy operation, allow a few minutes
		print "generating byte dictionary..."
		#conduitByteLocationDict = rpt.createByteLocDict("Link Results")
		rpt.createByteLocDict("Link Results")
		rpt.createByteLocDict("Node Results")

		#save dict to disk
		dictSaveFile =  open (byteLocDictionaryFName, 'w')
		pickle.dump(rpt.elementByteLocations, dictSaveFile)
		dictSaveFile.close()

	else:
		rpt.elementByteLocations = pickle.load( open(byteLocDictionaryFName, 'r') )
		#rpt.byteLocDict = conduitByteLocationDict

	print "Started Drawing at " + strftime("%b-%d-%Y %H:%M:%S")
	log = "Started Drawing at " + strftime("%b-%d-%Y %H:%M:%S") + "\n\nErrors:\n\n"
	drawCount = 0
	conduitErrorCount = 0

	font = ImageFont.truetype(su.fontFile, 30)
	basemapFeatureDicts = createFeaturesDict(options=ops['basemap'],  bbox=bbox, shiftRatio=shiftRatio, width=width) #will be populated after first frame is produced
	while currentT <= endDtime:

		img = Image.new('RGB', imgSize, ops['bg'])
		draw = ImageDraw.Draw(img)
		currentTstr = currentT.strftime("%b-%d-%Y %H:%M:%S").upper()
		#print 'time =', currentTstr

		#DRAW THE BASEMAP
		if ops['basemap']:
			drawBasemap(draw, img=img, options=ops['basemap'], width=width, bbox=bbox, featureDicts=basemapFeatureDicts, xplier = xplier)

		#DRAW THE CONDUITS
		if ops['conduitSymb']:
			for conduit, coordPairDict in conduitDicts.iteritems():
				#coordPair = coordPairDict['coordinates']
				if 'maxflow' in coordPairDict: #this prevents draws if no flow is supplied (RDII and such)

					su.drawConduit(conduit, coordPairDict, draw, options=ops['conduitSymb'],  rpt=rpt, dTime = currentTstr,
									xplier = xplier)

					drawCount += 1

				if drawCount > 0 and drawCount % 2000 == 0: print str(drawCount) + " pipes drawn - simulation time = " + currentTstr

		#DRAW THE NODES
		if ops['nodeSymb']:
			for node, nodeDict in nodeDicts.iteritems():
				if 'floodDuration' in nodeDict: #this prevents draws if no flow is supplied (RDII and such)
					su.drawNode(node, nodeDict, draw, rpt=rpt, dTime=currentTstr, options=ops['nodeSymb'], xplier=xplier)
					drawCount += 1


		#DRAW THE ANNOTATION
		dtime = currentT.strftime("%b%d%Y_%H%M").upper()
		#su.drawAnnotation (draw, inp, rpt, imgWidth=width, title=None, currentTstr = currentTstr, fill=su.black)
		su.annotateMap (draw, model, options=ops, currentTstr = currentTstr)
		currentT += delta
		del draw


		#shrink for antialiasing niceness (though this blows up the file about 5x)
		tempImgDir = os.path.join(imgDir, "temp_frames")
		saveImage(img, model, dtime, imgDir=tempImgDir, open=False, verbose=False)

	#WRITE THE GIF
	frames = []
	for image in glob.glob1(tempImgDir, "*.png"):
		imgPath = os.path.join(tempImgDir, image)
		frames.append(Image.open(imgPath))

	print "building gif with " + str(len(glob.glob1(tempImgDir, "*.png"))) + " frames..."
	if not imgName: imgName = inp.name
	gifFile = os.path.join(imgDir, imgName) + ".gif"
	frameDuration = 1.0 / float(ops['fps'])
	writeGif(gifFile, frames, duration=frameDuration)

	shutil.rmtree(tempImgDir) #delete temporary frames directory after succesful GIF

	log += "Completed drawing at " + strftime("%Y%m%d %H:%M:%S")
	with open(os.path.join(imgDir, "log.txt"), 'w') as logFile:
		logFile.write(log)

	print "Draw Count =" + str(drawCount)
	print "Video saved to:\n\t" + gifFile

	os.startfile(gifFile)#this doesn't seem to work


def drawProfile (model, upNode=None, dnNode=None, imgName=None, imgDir=None, width = 1024, height=512, drawSizeMultiplier = .01,
				bgroundColor=(0,3,18), bbox = None, conduitSubset = None, extraData=None):


	w = width*2
	h = height*2
	if upNode:
		#trace from up node to downstream stop node
		traced = su.traceFromNode(model, startNode=upNode, mode='down', stopnode=dnNode)
		conduitSubset = traced['conduits']
		nodeSubset = traced['nodes']

	#match the coordinates for each conduit
	conduitData = model.organizeConduitData(subset = conduitSubset)
	nodeData = model.organizeNodeData(subset=nodeSubset)
	conduits = conduitData['conduit_objects']
	nodes = nodeData['node_objects']
	minEl = nodeData['minEl']
	maxEl = nodeData['maxEl']
	reachLength = su.length_of_conduits(conduits)

	#transform
	profileHeight = maxEl - minEl

	sX = w / reachLength # to scale down from coordinate to pixels
	sY = h / profileHeight # to scale down from coordinate to pixels

	print w, " x " , h
	print 'profileHeight = ', profileHeight
	print "minEl = ",  minEl
	print "maxEl = ",  maxEl
	print "shiftRatioY = ",  sY
	print "shiftRatioX = ",  sX
	print 'reachLength = ', reachLength

	#return nodes
	img = Image.new('RGB', (w, h))
	draw = ImageDraw.Draw(img)
	x=0 #starting x location
	print '{}, {}, {}, {}, {}, {}'.format('conduit.id', 'x', 'upinvert', 'crown', 'flood_el', 'dnnode.maxHGL')
	for c in conduitSubset:

		#collect data about the curren conduit
		conduit = conduits[c]#grab object from dict\
		l = conduit.length
		upnode = nodes[conduit.upNodeID]
		dnnode = nodes[conduit.downNodeID]

		upfloodel = upnode.invert + upnode.maxDepth
		dnfloodel = dnnode.invert + dnnode.maxDepth

		inletoffset = conduit.inletoffset
		outletoffset = conduit.outletoffset
		geom1 = conduit.geom1
		upcondinv = upnode.invert + inletoffset
		dncondinv = dnnode.invert + outletoffset

		#compute the draw coordinates for each line segment
		upcondinvxy = su.coordToDrawCoord((x, upcondinv),shiftRatio=sX, shiftRatioY=sY, width=w, height=h)
		dncondinvxy = su.coordToDrawCoord((x+l, dncondinv),shiftRatio=sX, shiftRatioY=sY, width=w, height=h)
		#s = '{}: x={} upinv={}, dninv={}, uppeak={}, downpk={}'.format(conduit.id,x, upnode.invert, dnnode.invert, upnode.maxHGL, dnnode.maxHGL)
		s = '{}, {}, {}, {}, {}, {}'.format(conduit.id,x, upcondinv, upcondinv + geom1, upfloodel, upnode.maxHGL)
		x += l
		#print s
		print (upcondinvxy, dncondinvxy)
		draw.line((upcondinvxy, dncondinvxy), fill = (200, 200, 200), width = 2)



		#coordToDrawCoord(coordinates, bbox=None, shiftRatio=None, width=None, height=None)

	print '{}, {}, {}, {}, {}, {}'.format(conduit.id,x, upcondinv, upcondinv + geom1, upfloodel, upnode.maxHGL)
	saveImage(img, model, imgName)
