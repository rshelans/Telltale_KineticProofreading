#Modules Required by Program
import ij.plugin.frame.RoiManager as RoiManager
import ij
import ij.IJ
import ij.io
import ij.gui
import ij.gui.Plot
import java.awt.Color
import java.awt.event
import sys
import jarray.array
import os



def trackToRoi(file):
	"""
	Take a str formated: "x,y,intensity,frame,particle" get Rois
	"""
	track=[]
	for line in file:
		line    = [col for  col  in line.split() if col]
		if len(line) >= 4:#the number of expected columnsin a row
			frame   =int(float(line[3]))
			x       =int(float(line[0]))
			y       =int(float(line[1]))
			inensity=float(line[2])
 			roi=ij.gui.Roi(x-5, y-5, 11, 11)
 			roi.setPosition(frame+1)
			roi.setName("Track_" + str(frame+1))
 			track.append(roi)
 	return(track)

def trackToPlot(file):
	x=[]
	y=[]
	for line in file:
		line    = line.strip().split()
		if len(line) >= 4:#the number of expected columnsin a row
			x.append( int(float(line[3])) +1 )#frame
			y.append( float(line[2]) )#intensity		
	jx=jarray.array(x,'d')
	jy=jarray.array(y,'d')
	plot=ij.gui.Plot("Track","Frame","Intensity")
	plot.setLimits(min(x),max(x),min(y),max(y))
	plot.setColor(java.awt.Color.BLUE)
	plot.addPoints(jx,jy,ij.gui.Plot.LINE)
	return plot

class plotMouse(java.awt.event.MouseAdapter):
	def __init__(self,plot,stack):
		self.plot  =plot
		self.stack =stack
	def mouseClicked(self,e):
		if self.plot:
			x= int(self.plot.descaleX(e.getPoint().x))
			y= int(self.plot.descaleY(e.getPoint().y))
			self.stack.setPosition(x)
			#ij.IJ.log("("+str(x)+","+str(y)+") : "+str(self.stack ) )
		



stack   = ij.WindowManager.getCurrentImage()
manager = RoiManager.getInstance() or RoiManager()
#draws the plot of the track
#trackToPlot(track).show()
# sets the roi manager to have the track and cell calls

track = ij.io.OpenDialog("Choose A Track.")
track = open(track.getPath())

trk  = [line for line in track]
print(trk)
plot = trackToPlot(trk)
trk  = trackToRoi(trk)
track.close()
print(trk)

plot_window = plot.show()
manager.reset()
for roi in trk:
	print(roi.getPosition())
	manager.add(stack, roi, roi.getPosition())
	
plot_window.getCanvas().addMouseListener(plotMouse(plot,stack))
#print(track)
#ij.IJ.log("OPEN CONSOLE: "+str(stack))	