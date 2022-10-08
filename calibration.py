import ij.IJ
import ij.gui
import ij.plugin
import ij.plugin.frame
import ij.process
import itertools
import ij.measure
import ij.plugin.filter
##HELPER FUNCTIONS
##Image Cropping Utility
	
def selectAndCrop(roi, image, duplicator = None, manager = None):
	duplicator = duplicator or ij.plugin.Duplicator()
	manager    = manager or ij.plugin.frame.RoiManager().getInstance()
	rect       = roi.getBounds()
	ij.IJ.makeRectangle(rect.x, rect.y, rect.width, rect.height)
	#return(duplicator.crop(image))
	return(image.crop())

def normalize (image, images):
	##Measure ALl Rois in Image
	def measure(image, roi):
			measurement = ij.measure.Measurements.MEDIAN
			image.setRoi(roi)
			return(image.getStatistics(measurement).median)	
	measurements = [measure(image,roi) for roi in ROIs][1:-1]
	factor       = max(measurements)
	norms        = [factor/value for value in measurements]
	[image.getChannelProcessor().multiply(norm) for image,norm in zip(images, norms)]

def resetAndAddROIsToManager(ROIs, manager = None):
	###Gets instance of ROI manager and adds ROIs of Slices to It.
	manager = manager or ij.plugin.frame.RoiManager().getInstance()
	manager.reset()
	for i,roi in enumerate(ROIs):
		manager.addRoi(roi)

def stackRegistration(stack):
	## Registration Using TemplateMatching
	ij.IJ.run(stack, "Align slices in stack...", "method=5 windowsizex=88 windowsizey=85 x0=93 y0=94 swindow=0 subpixel=false itpmethod=0 ref.slice=1")
	results = ij.IJ.getTextPanel().getResultsTable()
	dx      = list(results.getColumn(results.getColumnIndex('dX')))
	dy      = list(results.getColumn(results.getColumnIndex('dY')))
	##Registration end
	##CLEANING
	ij.IJ.selectWindow("Results"); 
	ij.IJ.run("Close");
	ij.IJ.selectWindow("Log"); 
	ij.IJ.run("Close");
	return(dx,dy)
	

image  = ij.IJ.getImage()
frames = list(itertools.product(range(-1,2),range(-1,2)))[1:-1]

###gets the current Selection and returns a bounding box
ij.IJ.makeRectangle(365, 369, 273, 270)
rect   = image.getRoi().getBounds()

ROIs   = [ij.gui.Roi(rect.x+rect.width*j, rect.y+rect.height*i, rect.width, rect.height) for i,j in frames]


##Creates Stack From cropped images using ROIs but also normalzes by intensity
images = [selectAndCrop(roi,image) for roi in ROIs]#Removes first and last frame
normalize(image, images)
stack     = ij.plugin.ImagesToStack().run(images)
ij.IJ.run(stack, "8-bit", "");
## Stack created ##

dx,dy = stackRegistration(stack)
dx.insert(0,0)
dy.insert(0,0)

adjustedROIs = [ij.gui.Roi(rect.x+rect.width*j-dx, rect.y+rect.height*i-dy, rect.width, rect.height) for dx,dy, (i,j) in zip(dx,dy,frames)]
resetAndAddROIsToManager(adjustedROIs)

images = [selectAndCrop(roi,image) for roi in adjustedROIs]
stack  = ij.plugin.ImagesToStack().run(images)
stack.show()
