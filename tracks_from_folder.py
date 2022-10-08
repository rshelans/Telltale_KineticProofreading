#Modules Required by Program
import ij.plugin.frame.RoiManager as RoiManager
import ij
import ij.IJ
import ij.gui
import ij.io
import sys
import jarray.array
import os
import glob

def trackToRoi(file):
	"""
	Take a str formated: "x,y,intensity,frame,particle" get Rois
	"""
	track=[]
	for line in file:
		line    = line.split()
		if len(line) >= 5:#the number of expected columnsin a row
			frame   =int(float(line[3]))
			x       =int(float(line[0]))
			y       =int(float(line[1]))
			inensity=float(line[2])
 			roi=ij.gui.Roi(x-5,y-5,10,10)
 			roi.setPosition(frame+1)
			roi.setName("Track_"+str(frame+1))
 			track.append(roi)
 	return(track)


manager= RoiManager.getInstance() or RoiManager()
manager.reset()
directory=ij.IJ.getDirectory("Select Track Directory.")
tracks=glob.glob(directory+"*.trk")
for track in tracks:
		print(os.path.basename(track).split('.')[0])
		track=open(track)
		trk=   trackToRoi(track)
		track.close()
		for roi in trk:
			manager.addRoi(roi)
		