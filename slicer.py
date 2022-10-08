import ij.IJ
import ij.gui
import ij.plugin
import ij.process
import itertools
import ij.measure
import ij.plugin.filter

	
def selectAndCrop(roi, image, options = 'slice'):
	rect       = roi.getBounds()
	image.setRoi(rect.x, rect.y, rect.width, rect.height)
	return(image.crop(options))

def sliceImage(image, ROIs):
	images  = [selectAndCrop(roi, image, options = 'stack') for roi in ROIs]
	##NORMALIZATION
	ij.IJ.run(images[3],'Multiply...',"value=0.685626 stack")
	##ENDNORMALIZATION
	stack   = ij.plugin.Concatenator.run(images)
	return stack


def main():
	image   = ij.IJ.getImage()
	manager = ij.plugin.frame.RoiManager().getInstance()
	ROIs    = manager.getRoisAsArray()
	stack   = sliceImage(image,ROIs)
	c,z,t = 1,7,image.getStackSize()
	order = "xyctz"
	HyperStack         = ij.plugin.HyperStackConverter().toHyperStack(stack,c,z,t, order , "grayscale")
	#ij.IJ.run(stack, "Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
	#ij.plugin.Hyperstack_rearranger
	#stack     = fiji.stacks.Hyperstack_rearranger().reorderHyperstack(stack,"CTZ",False,False)
	HyperStack.show()

if __name__ == "__builtin__":
	main()