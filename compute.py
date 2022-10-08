import ij.IJ
import ij.gui
import ij.plugin
import ij.plugin.filter
import ij.plugin.frame
import ij.process
import itertools
import ij.measure
import ij.io

import ij.WindowManager
import sys
sys.path.append("C:\\Users\\Robert Shelansky\\Desktop\\")
import slicer


###DIALOGUE##############################################
images = ij.WindowManager.getImageTitles()
dialog = ij.gui.GenericDialog('Title')
dialog.addChoice('Background Image',images,images[0])
dialog.addChoice('Stack'           ,images,images[0])
dialog.addStringField('output basename  :', 'Strain_Date_Experiment',72)
dialog.addStringField('Imaging Condition:', '30Power ND1 SSEMAMP GreenMFMonly 250ms Exposure', 72)
dialog.showDialog()
if dialog.wasCanceled():
    sys.exit()

backgroundTitle = dialog.getChoices()[0].getSelectedItem()
stackTitle      = dialog.getChoices()[1].getSelectedItem()
basename        = dialog.getStringFields()[0].getText()
description     = dialog.getStringFields()[1].getText()

backgroundImage =  ij.WindowManager.getImage(backgroundTitle)
stackImage      =  ij.WindowManager.getImage(stackTitle)
##END DIALOGUE######################################

##BACKGROUND SUBTRACT THE IMAGE
subtractedImage = ij.plugin.ImageCalculator().run("Subtract create stack",stackImage,backgroundImage)

##BLEACH CORRECTTION
ij.IJ.run(subtractedImage, "Bleach Correction", "correction=[Exponential Fit]")
bleachedImage    = ij.WindowManager.getCurrentImage()
ij.WindowManager.getImage("y = a*exp(-bx) + c").close()
logText = str(ij.IJ.getLog())
ij.IJ.beep()
ij.IJ.log(r"\\Clear")
ij.WindowManager.getWindow("Log").close()




manager = ij.plugin.frame.RoiManager().getInstance()
manager.deselect()
ROIs    = manager.getRoisAsArray()

##Slice Images
c,z,t = 1,7,stackImage.getStackSize()
order = "xyctz"
bleachedSlice   = slicer.sliceImage(bleachedImage  , ROIs)
Slice           = slicer.sliceImage(subtractedImage, ROIs)
bleachedImage  .close()
subtractedImage.close()

bleachedHyperstack = ij.plugin.HyperStackConverter().toHyperStack(bleachedSlice,c,z,t, order , "grayscale")
hyperstack         = ij.plugin.HyperStackConverter().toHyperStack(Slice,c,z,t, order , "grayscale")
bleachedSlice.close()
Slice        .close()
##Slicing Completed

#Z-projections
maxBleachedHyperstack= ij.plugin.ZProjector().run(bleachedHyperstack,'maxall')
maxHyperstack        = ij.plugin.ZProjector().run(hyperstack        ,'maxall')
#sumBleachedHyperstack= ij.plugin.ZProjector().run(bleachedHyperstack,'sumall')


##FILE HANDLING AT THE END
stackImage           .setTitle(str(basename))
bleachedHyperstack   .setTitle(str(basename)+"_BS_BC_SLICED")
hyperstack           .setTitle(str(basename)+"_BS_SLICED")
maxBleachedHyperstack.setTitle(str(basename)+"_BS_BC_SLICED_MAX")
maxHyperstack        .setTitle(str(basename)+"_BS_SLICED_MAX")
#sumBleachedHyperstack.setTitle(str(basename)+"_BS_BC_SLICED_SUM")

chooser = ij.io.DirectoryChooser("Select Output Directory.")
path    = chooser.getDirectory()
#path = "C:\\Users\\Robert\\Desktop\\" 
try:
	textFile = open(path + basename + ".txt", 'w')
	text = description + "\n\nMiddle Frame Adjusted 0.685626\n\n PhotoBleachCorrection:\n\n" + logText
	textFile.write(text)
except:
	textFile.close()
textFile.close()

ij.io.FileSaver(hyperstack)           .saveAsTiffStack(path + hyperstack.getTitle()            + ".tif")
ij.io.FileSaver(stackImage)           .saveAsTiffStack(path + stackImage.getTitle()            + ".tif")
ij.io.FileSaver(maxHyperstack)        .saveAsTiffStack(path + maxHyperstack.getTitle()         + ".tif")
ij.io.FileSaver(bleachedHyperstack)   .saveAsTiffStack(path + bleachedHyperstack.getTitle()    + ".tif")
ij.io.FileSaver(maxBleachedHyperstack).saveAsTiffStack(path + maxBleachedHyperstack.getTitle() + ".tif")
#ij.io.FileSaver(sumBleachedHyperstack).saveAsTiffStack(path + sumBleachedHyperstack.getTitle() + ".tif")
#END Handling


##SHOW IMAGES
bleachedHyperstack   .show()
hyperstack           .show()
maxBleachedHyperstack.show()
maxHyperstack        .show()
sumBleachedHyperstack.show()
