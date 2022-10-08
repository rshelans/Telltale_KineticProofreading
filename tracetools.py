import scipy
import pandas
class GaussianMask:
	"""
	PsfFit fits a gaussian mask to find the centroid of individual particles, algorithm from: 
	'Precise Nanometer Localization Analysis for Individual Fluorescent Probes' Russell E. Thompson, Daniel R. Larson, and Watt W. Webb
	"""
	
	def __init__(self, window_radius=25, psf_width=1.7):
		##FORCE  THE FIT  WINDOW TO  BE An ODD NUMBER
		self.window_radius = int(window_radius)
		self.window_size   = self.window_radius  * 2 + 1 
		self.psf_width     = float(psf_width)
		self.i,self.j      = scipy.meshgrid(scipy.r_[:self.window_size]-int(0.5*self.window_size) ,scipy.r_[:self.window_size]-int(0.5*self.window_size))
		self.N             = scipy.exp(-(self.i)**2/(2*psf_width**2)-(self.j)**2/(2*psf_width**2))/(2*scipy.pi*psf_width**2)
		self.N2            = scipy.sum(self.N**2)

	class GaussianMaskFit:
		"""
		GaussianMaskFit Represents a Fit that has been Completed.
		A fit that is incomplete has a centroid of 0,0 and intensity 0.
		"""
		def __init__(self,centroid=scipy.array([0.0,0.0]),intensity=0.0):
			self.centroid   = centroid
			self.intensity  = intensity

	def fit(self, image, coordinate, iterations=200, delta = 0.1, fit=None):
		"""
		Fit applys the GaussianMask to to a specific cooredinate.
		"""
		if scipy.amax(scipy.r_[coordinate-self.window_radius< 0, coordinate + self.window_radius > image.shape[::-1]]):
				raise ValueError("A {}px window centered at {}. Does not fit on image of shape {}.".format(self.window_size,  coordinate, image.shape))
		#Check to see if a Fit has already been attempted for this mask.
		fit = fit if fit else GaussianMask.GaussianMaskFit(coordinate)
		for i in range(iterations):
			#origin (topleft corner of box) of a window around the coordinate
			#The slice of the image by the window
			x , y = fit.centroid.astype(int) 
			S     = image[y-self.window_radius:y+self.window_radius+1,x-self.window_radius:x+self.window_radius+1]
			#calculate the centroid from the current fit
			# print(x,y,S.shape,image.shape)
			SN              = self.N * S
			centroid        = fit.centroid - scipy.array([scipy.sum(self.i*SN)/scipy.sum(SN),scipy.sum(self.j*SN)/scipy.sum(SN)])
			fit.intensity   = scipy.sum(SN) / self.N2
			fit._iterations_= i+1
			#return fit if converged
			if scipy.absolute(fit.centroid-centroid).max() <= delta:
				self.S        = S
				self.SN       = SN
				self.centroid = centroid
				return(fit)
			fit.centroid=centroid
		raise ValueError("A fit at coordinate {} did not converge (delta={}, iterations={}).".format(coordinate,delta,iterations))


##Pretty Sure these suck... But afraid to depricate...
class DefaultBackground:
		"""
		DefaultBackground is the base class for the Background Interface.
		All background corrections must have an apply function.
		Which takes two arguments the image to correct and the coordinate to
		which to apply the coorection.
		"""
		def apply(self,im,coordinate):
			"""
			This class simply returns the image. No Background Correction is Completed.
			"""
			return(im)

class FlatBackground:
	"""
	FlatBackground is the base class for the Background Interface.
	All background corrections must have an apply function.
	"""
	def __init__(self,plane_size,window_size=50):
		"""
		Flat Backgrounds use A GaussianMask reference to decide how to flatten the background.
		"""
		self.window_size = window_size
		self.plane_size  = plane_size

	def apply(self,im,coordinate):
		"""
		Subtracts the mean background of a square around a coordinate.
		"""
		if not scipy.prod(coordinate-self.plane_size/2.>=0)*scipy.prod(coordinate+self.plane_size/2.<=im.shape[::-1]):
			ValueError("A {}px window centered at {}. Does not fit on image of shape {}.".format(self.window_size,coordinate,im.shape))
			
		im = im.copy()
		x,y=(coordinate-int(self.plane_size)/2).astype(int)#origin of a box top left corner. 
		bg =scipy.median((im[y,x:x+self.plane_size],#top Row
			im[y+self.plane_size-1,x:x+self.plane_size],#bottom Row
			im[y:y+self.plane_size,x],#Left Column
			im[y:y+self.plane_size,x+self.plane_size-1]))#Right Column
		
		im[y:y+self.window_size,x:x+self.window_size] = (im[y:y+self.window_size,x:x+self.window_size]- bg).clip(min=0)
		return (im)

class PlaneBackground:
	"""
	PlaneBackground uses a polynomial fit to the square around a coordinate to subtract a planer background.
	"""
	def __init__(self,window_size=50):
		"""
		Plane Backgrounds use A GaussianMask reference to decide how to flatten the background.
		"""
		self.window_size=window_size

	def apply(self,im,coordinate):
		"""
		Subtracts the mean background of a square around a coordinate.
		"""
		if not scipy.prod(coordinate-self.window_size/2.>=0)*scipy.prod(coordinate+self.window_size/2.<=im.shape[::-1]):
			ValueError("A {}px window centered at {}. Does not fit on image of shape {}.".format(self.window_size,coordinate,im.shape))
		im        = im.copy()
		x,y       = (coordinate-int(self.window_size)/2).astype(int)
		xy        = scipy.r_[:2*self.window_size]%self.window_size-(self.window_size-1)/2
		mgx       = xy[:self.window_size] * scipy.polyfit(xy,scipy.r_[im[y:y+self.window_size,x],im[y:y+self.window_size,x+self.window_size-1]],1)[0]
		mgx.shape = (1,self.window_size)
		mgy       = xy[:self.window_size]*scipy.polyfit(xy,scipy.r_[im[y,x:x+self.window_size],im[y+self.window_size-1,x:x+self.window_size]],1)[0]
		mgy.shape = (self.window_size,1)

		bg = scipy.mean([
			im[y,x:x+self.window_size],#top Row
			im[y+self.window_size-1,x:x+self.window_size],#bottom Row
			im[y:y+self.window_size,x],#Left Column
			im[y:y+self.window_size,x+self.window_size-1]])#Right Column
		
		im[y:y+self.window_size,x:x+self.window_size]=(im[y:y+self.window_size,x:x+self.window_size]-mgx-mgy-bg).clip(min=0)
		return (im)

##Maybe deprecated to the bottom methods
# class Loci(pandas.DataFrame):
# 	def __init__(self,*args,**kw):
# 		super(Loci,self).__init__(*args,**kw)

# 	def track(self,path_start,start,end,max_distance=50):
# 		"""
# 		track finds the brightest loci at most max_distance from path_start in frames start:end.
# 		track will insert an puncta with intensity zero and x,y coordinate equal to the previous
# 		frame for all frames which contains no puncta within the desired distance.

# 		loc is a pandas data frame with the following structure:
# 			#loc.columns =["x","y","intensity","frame"]
# 		-path_start is a tuple (x,y), it is the anchor for all distance comparisons:
# 		-start  is the frame to begin the path.
# 		-end is the frame after the frame to end the path on	.
# 		-max_distance is the maximum distance a puncta can be from the path_start to be considered for a path.

# 		Returns a pandas DataFrame of the form track.columns =["x","y","intensity","frame"] for all frames start:end.
# 		"""
# 		#Distnce of everypoint from the path_start
# 		self["distance"] = ((path_start[0]-self.x.values)**2+(path_start[1]-self.y.values)**2)**0.5
# 		#Only those dots within distance from path_start To limit search space:limit to specific frames
# 		track=self[(self.distance < max_distance) & (self.frame >= start) & (self.frame < end)]
# 		#Grab the brightest dot in each frame
# 		track=track.sort_values('intensity',ascending=False).groupby("frame",as_index=False).first()
# 		return(Loci(track[["x","y","intensity","frame"]]))

# 	def fill(self,start,end):
# 		"""
# 		fill takes a track with at most 1 dot in each frame and fills in frames inbatween with the x,y coordinate
# 		of the last seen puncta. This will crash if multiple rows from the same frame are present.
# 		A valueError: cannot reindex from duplicate axis error from pandas proper.
# 		"""
# 		#add in empty frames
# 		track=self.set_index("frame").reindex(pandas.Index(range(start,end,1),name="frame")).reset_index()
# 		#fill the nas with the preceeding value
# 		track[['x','y']]=track[['x','y']].fillna(method='ffill')
# 		return(track)

# 	def track_from_loci(self,path,frames,distance):	
# 		assert(len(path)==len(frames) and len(path)==len(distance))
# 		tracks=(self.track(path[i],frames[i][0],frames[i][1],max_distance=distance[i]) 
# 			for i in range(len(path)))
# 		return(pandas.concat(tracks))

# 	def write_loc(self,file,sep="\t"):
# 		super(Loci, self[['x','y','intensity','frame']]).to_csv(file,sep=sep,header=False,na_rep=0,index=False)
		
# 	def write_trk(self,file,sep="\t"):
# 		track=self
# 		track["particles"]=0
# 		super(Loci, self[['x','y','intensity','frame','particles']]).to_csv(file,sep=sep,header=False,na_rep=0,index=False)
# 		del track["particles"]

# 	@property
# 	def _constructor(self):
# 		"""
# 		Required For Subclassing Pandas.DataFrame
# 		""" 
# 		return(Loci)

# 	@classmethod
# 	def read_csv(cls,file,sep="\t"):
# 		new_loci=pandas.read_csv(file,sep=sep,header=None)[[0,1,2,3]]
# 		new_loci.columns =["x","y","intensity","frame"]
# 		return(Loci(new_loci))

##Messing with ROIS From ImageJ
# import  read_roi
# rois = read_roi.read_roi_zip(r'E:\raw\Heta\20180513\YTL362D23_2\cell1.zip')
def getOvalRois(rois):
    for roi in rois.values():
        left,top,height,width,position = [roi[feature]  for feature in ['left','top','height','width','position']]
        yield(position-1, (left+width/2, top+width/2, width, height, top, left, position))

def forwardFill(rois, last_frame=None):
	first_frame = min(rois)
	last_frame  = last_frame  or max(rois) + 1
	current     = rois[min(rois)]
	for i in range(first_frame,last_frame):
		current = rois.get(i) or current
		yield i, current

def reverseFill(rois, first_frame=None):
	first_frame = first_frame if first_frame is not None else min(rois)
	last_frame  = max(rois) + 1
	current     = rois[max(rois)]
	for i in reversed(range(first_frame,last_frame)):
		current = rois.get(i) or current
		yield i, current

def localize(rois, loc):
    inOval = lambda x, y, cx, cy, w, h: (x-cx)**2/ (w/2)**2 + (y-cy)**2/ (h/2)**2 <= 1
    for frame in range(min(rois),max(rois)+1):
        thisRoi           = rois[frame]
        thisLoc           = loc[loc[:,3]==frame]
        x, y, cx, cy, w,h = thisLoc[:,0], thisLoc[:,1], *thisRoi[:4]
        potentialPuncta   = thisLoc[inOval(x,y,cx,cy,w,h)]
        if len(potentialPuncta):
            yield(max(potentialPuncta, key = lambda puncta: puncta[2]))

def trackFromCellKeyFrames(rois, loc, cell_number = 0):
	"""
	Takes Key Frames of Cell location and generates a track of the brightest puncta in each key frame.
	Feeding Forward frames with no puncta.
	"""
	rois         = dict(forwardFill(dict(getOvalRois(rois))))
	track        = scipy.array(list(localize(rois,loc)))
	##IN CASE THERE IS NO PUNCTA SEEN IN THE WHOLE TRACE
	track        = [[roi[2],roi[3],0,frame] for frame,roi in rois.items()] if len(track) == 0 else track
	#trad         = {puncta[3] : puncta for puncta in track}
	fdtrack      = dict(forwardFill( {int(puncta[3]) : list(puncta) for puncta in track}, last_frame=max(rois)+1))
	fdtrack      = dict(reverseFill(fdtrack,first_frame = min(rois)))
	fftrack      = scipy.array(list(map(fdtrack.__getitem__, sorted(fdtrack))))
	frames       = sorted(fdtrack)
	rough_track        = scipy.zeros((fftrack.shape[0],fftrack.shape[1]+1))
	rough_track[:,:-1] = fftrack
	rough_track[:,-1 ] = cell_number
	rough_track[:, 3 ] = frames
	return (rough_track)

# def _discover(im,psfPx=1.7,threshold=6.5,min_distance=1):
# 	imBpass   =imagetools.bpass(im,1.,psfPx).clip(0)
# 	maxima    =skimage.feature.peak_local_max(
# 		imBpass,
# 		min_distance=min_distance,
# 		threshold_abs=threshold*scipy.var(imBpass)**.5,
# 		indices=True)
# 	track=[]
# 	for coordinate in maxima[:,[1,0]]:
# 		track.append([coordinate[0],coordinate[1],0.0])
# 	track= pandas.DataFrame(track)
# 	track.columns=['x','y','intensity']
# 	return(track)


# def discover(movie,psfPx=1.7,threshold=6.5,min_distance=1):
# 	track=[]
# 	for frame,im in enumerate(movie):
# 		imBpass   =imagetools.bpass(im, 1., psfPx).clip(0)
# 		maxima    =skimage.feature.peak_local_max(imBpass,
# 						min_distance=min_distance,
# 						threshold_abs=threshold*scipy.var(imBpass)**.5,
# 						indices=True)
		
# 		for coordinate in maxima[:,[1,0]]:
# 			track.append([coordinate[0],coordinate[1],0.0,frame])
	
# 	track=Loci(track)
# 	track.columns=['x','y','intensity','frame']
# 	return(track)

# def punctaAnalyzer(psf = 1.7,  window_radius = 25):
# 	wr         = window_radius
# 	ws         = wr*2+1
# 	background = PlaneBackground(window_size   = ws)
# 	mask       = GaussianMask   (window_radius = wr, psf_width = psf)
# 	def analyzePuncta(image, x,  y):
# 		x,y        = int(x), int(y)
# 		coordinate =scipy.array([x,y]).astype(int)
# 		if all(coordinate - wr >  0) and all(coordinate + wr < image.shape[::-1]):
# 			centroid   = scipy.array([wr,wr]).astype(int)
# 			mini       = image[y-wr-1:y+wr,x-wr-1:x+wr]
# 			backed     = background.apply(mini, centroid)
# 			bpassed    = bpass(backed, r2 = psf).clip(0)
# 			fit        = mask.fit(bpassed, centroid, iterations = 1, delta = scipy.inf)
# 			return(fit.intensity)
# 		return(scipy.nan)
# 	return(analyzePuncta)

# def reAnalyzeTrack(stack, track, psfPx=1.5, window_radius=21):
# 	analyzer = punctaAnalyzer(psf=psfPx, window_radius=window_radius)
# 	for i in  range(min(track[:,3]).astype(int), max(track[:,3].astype(int))):
# 		x, y, _, frame,cell = track[i]
# 		yield(x, y, analyzer(stack[i], x, y), frame, cell )


##TRACK MANIPULATION
def mergeTracks(tracks):
	for cell, track in enumerate(tracks):
		track[:,4] = cell
	return(scipy.concatenate(tracks))

##For reanalyzing puncta
##Deperacated For Normalize and Analyze on 12/05/2019.
# def reAnalyzePuncta(stack, loc, psfPx = 1.5, window_radius = 40, smooth= 7):
# 	wr   = window_radius
# 	mask = GaussianMask   (window_radius = wr, psf_width = psfPx)
# 	for frame, image in enumerate(stack):
# 		image     = scipy.array(image).astype(int)
# 		smoothed  = scipy.ndimage.filters.median_filter(image, size=smooth)
# 		corrected = (image - smoothed).clip(0)
# 		bpassed   = bpass(corrected, 1., psfPx).clip(0)
# 		for puncta in loc[loc[:,3]==frame]:
# 			x,y,intensity,f,cell =  puncta
# 			x,y,f,cell = int(x),int(y),int(f),int(cell)
# 			coordinate =scipy.array([x,y])
# 			assert(f == frame)
# 			if all(coordinate - wr >  0) and all(coordinate + wr < image.shape[::-1]):
# 				centroid   = scipy.array([wr,wr])
# 				mini       = bpassed[y-wr-1:y+wr,x-wr-1:x+wr]
# 				fit        = mask.fit(mini, centroid, iterations = 1, delta = scipy.inf)
# 				yield(x, y, fit.intensity, frame,cell)
# 			else:
# 				yield(x, y, scipy.nan    , frame,cell)


##Functions to normalize and analyze movies: Acts on individual images at a time can be used with pool for fast 
##processing
##E.G.
# def analyze_helper(movie_url, nuc_url, turls, pool = None, normalize = None, analyze = None):
#     pool      = pool      or multiprocessing.Pool(processes=16)
#     normalize = normalize or test.Normalize(psfPx = 1.5, window_radius = 40, smooth = 7)
#     analyze   = analyze   or test.Analyze  (psfPx = 1.5, window_radius = 40, smooth = 7)

#     with pims.open(movie_url) as stack:
#         normalized = pool.map(normalize, stack)
        
#         tracks    = [scipy.loadtxt(track, delimiter='\t',converters={0:int,1:int,2:float,3:int,4:int}) for track in turls]
#         loc       = scipy.array(sorted(tracetools.mergeTracks(tracks), key = lambda x:x[3]))      
#         nucs      = scipy.loadtxt(nuc_url, delimiter='\t', converters = {0:int,1:float,2:int,3:int,4:int,5:float,6:int,7:float})
#         #[[image,loc] for image in normalized]
#         iterable  = [[image/scipy.median(nucs[nucs[:,3]==frame],axis=0)[0],loc[loc[:,3] == frame]] for frame,image in enumerate(normalized)]
#         loc       = pandas.concat(pool.map(analyze, iterable))
#     return loc
class Normalize:
    def __init__(self, psfPx = 1.5, window_radius = 40, smooth = 7, **kwargs):
        self.psfPx         = psfPx
        self.window_radius = window_radius
        self.smooth        = smooth

    def __call__(self, image):
        image     = scipy.array(image).astype(int)
        smoothed  = scipy.ndimage.filters.median_filter(image, size=self.smooth)
        corrected = (image - smoothed).clip(0)
        bpassed   = bpass(corrected, 1., self.psfPx).clip(0)##########
        return(bpassed)

class Analyze:
    def __init__(self, psfPx = 1.5, window_radius = 40, smooth = 7, **kwargs):
        self.psfPx         = psfPx
        self.window_radius = window_radius
        self.smooth        = smooth
        self.mask          = GaussianMask   (window_radius = window_radius, psf_width = psfPx)#########

    def __call__(self, item):
        frame = pandas.DataFrame(list(self.quantifyPuncta(*item)))
        frame.columns = ['x','y','intensity','frame','cell']
        return (frame)

    def quantifyPuncta(self, image, loc):
        wr = self.window_radius
        for puncta in loc:
            x,y,intensity,f,cell = puncta
            x,y,f,cell = int(x),int(y),int(f),int(cell)
            coordinate =scipy.array([x,y])
            if all(coordinate - wr >  0) and all(coordinate + wr < image.shape[::-1]):
                centroid   = scipy.array([wr,wr])
                mini       = image[y-wr-1:y+wr,x-wr-1:x+wr]
                fit        = self.mask.fit(mini, centroid, iterations = 1, delta = scipy.inf)
                yield(x, y, fit.intensity, f,cell)
            else:
                yield(x, y, scipy.nan    , f,cell)



##Finding potential puncta and scoreing them on intensity
import skimage.feature
def discover(image, psfPx=1.7, threshold = 6.5, min_distance=1):
	imBpass   = bpass(image, 1., psfPx).clip(0)
	maxima    = skimage.feature.peak_local_max(imBpass,
										min_distance=min_distance,
										threshold_abs=threshold*scipy.var(imBpass)**.5,
										indices=True)
	for x, y in maxima[:,[1,0]]:
		yield(x, y)

def getAllPunctaFromStack(stack, psfPx = 1.5, window_radius = 40, threshold = 5):
	wr   = window_radius
	mask = GaussianMask   (window_radius = wr, psf_width = psfPx)
	for frame, image in enumerate(stack):
		smoothed  = scipy.ndimage.filters.median_filter(image.astype(int),  size=7)
		corrected = (image.astype(int) - smoothed).clip(0)
		bpassed   = bpass(corrected, r2 = psfPx).clip(0)
		for x, y in discover(image, psfPx = psfPx, threshold = threshold, min_distance=1):
			coordinate = scipy.array([x,y]).astype(int)
			if all(coordinate - wr >  0) and all(coordinate + wr < image.shape[::-1]):
				centroid   = scipy.array([wr,wr])
				mini       = bpassed[y-wr-1:y+wr,x-wr-1:x+wr]
				fit        = mask.fit(mini, centroid, iterations = 1, delta = scipy.inf)
				yield(x, y, fit.intensity, frame)



import scipy.fftpack
###THIs comes from coulon et al.
####Computes the band pass of an image by calling as the following 
"""
####Example USE
##imBpass=imagetools.bpass(im,1.,psfPx)                       # Band-passed image
##tifffile.imsave("./out/band.tif", (imBpass).astype(int) )
##imBinary=(imBpass>thresholdSD*scipy.var(imBpass)**.5)*1.     # Binary image
##tifffile.imsave("./out/bindary.tif", (imBinary).astype(int) )
THIS CODE IS BURROWED FROMTINEKE WHO GOT IT FROM Dan  larsons student
"""

###FFT pack will shortly not be compatible with scipy... so everything is going to break eventually....
sHS=scipy.fftpack.fftshift # Swap half-spaces. sHS(matrix[, axes]). axes=all by default
def hS(m,axes=None):
		if axes==None: 
				axes=range(scipy.ndim(m))
		elif type(axes)==int: 
				axes=[axes]
		elif axes==[]: 
				return m
		return hS(m.swapaxes(0,axes[-1])[:m.shape[axes[-1]]/2].swapaxes(0,axes[-1]),axes[:-1])

def sHSM(m,axes=None):
		if axes==None: 
				axes=range(scipy.ndim(m))
		elif type(axes)==int: 
				axes=[axes]
		m=m.swapaxes(0,axes[0])
		max=m[1]+m[-1]
		m=(m+max/2)%max-max/2
		m=m.swapaxes(0,axes[0])
		return sHS(m,axes)

def bpass(im,r1=1.,r2=1.7):
	"""
	###FFT pack will shortly not be compatible with scipy... so everything is going to break eventually....
	"""
	ker1x=scipy.exp(-(sHS(sHSM(scipy.r_[:im.shape[1]]))/r1)**2/2)
	ker1x/=sum(ker1x)
	fker1x=scipy.fftpack.fft(ker1x)
	ker1y=scipy.exp(-(sHS(sHSM(scipy.r_[:im.shape[0]]))/r1)**2/2)
	ker1y/=sum(ker1y)
	fker1y=scipy.fftpack.fft(ker1y)
	ker2x=scipy.exp(-(sHS(sHSM(scipy.r_[:im.shape[1]]))/r2)**2/2)
	ker2x/=sum(ker2x)
	fker2x=scipy.fftpack.fft(ker2x)
	ker2y=scipy.exp(-(sHS(sHSM(scipy.r_[:im.shape[0]]))/r2)**2/2)
	ker2y/=sum(ker2y)
	fker2y=scipy.fftpack.fft(ker2y)
	fim=scipy.fftpack.fftn(im)
	return scipy.fftpack.ifftn((fim*fker1x).T*fker1y-(fim*fker2x).T*fker2y).real.T



import skimage.draw
def draw_box(image, x, y ,size=7,bit_depth=16):
	y_points=scipy.array([y+size,y+size,y-size,y-size])
	x_points=scipy.array([x+size,x-size,x-size,x+size])
	rr,cc   =skimage.draw.polygon_perimeter(y_points, x_points, shape=image.shape)
	image[rr,cc] = int(2**bit_depth)-1



###Toolkit for analyzing Nuclei Is ok not great.
class Nucleous:
	import itertools
	import scipy.spatial
	import scipy.ndimage

	@staticmethod
	def getPunctaCoords(loc, frame):
		coords = loc[ loc[:,3] == frame, 0:2 ].astype(int)
		return(zip(*coords))

	@staticmethod   
	def pick(labels, coords):
		centers      = scipy.around(scipy.ndimage.center_of_mass(labels, labels, index=scipy.arange(1,1+scipy.amax(labels)))).astype(int)
		tree         = scipy.spatial.cKDTree(centers)
		_,identity   = tree.query(coords)
		return(identity+1)

	@staticmethod
	def summary(image, labels):
		median = scipy.ma.median(scipy.ma.array(image, mask=labels == 0))
		var    = scipy.ma.var(scipy.ma.array(image, mask=labels == 0))
		bmedian= scipy.ma.median(scipy.ma.array(image, mask=labels > 0))
		bvar   = scipy.ma.var(scipy.ma.array(image, mask=labels > 0))
		return(median, var, bmedian, bvar)

	@staticmethod    
	def nucleiInfo(image, labels, identity):
		f    = lambda x: [scipy.ma.median(x), scipy.ma.var(x)]
		info = [f(scipy.ma.array(image, mask = labels!=i))   for i in identity]
		med, var = zip(*info)
		return(scipy.array(med),scipy.array(var))

	@staticmethod    
	def  label(image, sigma = 2):   
		mew         = scipy.mean(image)
		var         = scipy.var (image)
		binary = image > mew + sigma * scipy.sqrt(var)
		return binary 

	@staticmethod
	def identify(image,sigma=2, cell_width=15, smooth=4):
		fcells      = scipy.ndimage.filters.gaussian_filter(image,          cell_width)
		fnuclei     = scipy.ndimage.filters.gaussian_filter(fcells - image , smooth)
		nlabels     = Nucleous.label(fnuclei,sigma)
		labels,_   = scipy.ndimage.label(nlabels)
		return(labels)

	@staticmethod
	def buildFrame(info, desc, identity,frame):
		info          = pandas.DataFrame(info).T
		info.columns  = ['median','variance']
		info['cell' ] = scipy.arange(len(identity))+1
		info['frame'] = frame
		info['allMedian'],info['allVariance'],info['backgroundMedian'],info['backgroundVariance'] = desc
		return(info)

	@staticmethod
	def processImage(image, loc, frame):
		labels        = Nucleous.identify(image, sigma=1.7)
		x,y           = Nucleous.getPunctaCoords(loc, frame)
		identity      = Nucleous.pick(labels, list(zip(y,x)))
		info          = Nucleous.nucleiInfo(image,labels,identity)
		desc          = Nucleous.summary(image,labels) 
		info          = Nucleous.buildFrame(info, desc, identity, frame)
		return(info)

	@staticmethod
	def processImageMap(parms):
		image, loc, frame = parms
		labels        = Nucleous.identify(image, sigma=1.7)
		x,y           = Nucleous.getPunctaCoords(loc, frame)
		identity      = Nucleous.pick(labels, list(zip(y,x)))
		info          = Nucleous.nucleiInfo(image,labels,identity)
		desc          = Nucleous.summary(image,labels) 
		info          = Nucleous.buildFrame(info, desc, identity, frame)
		return(info)



##Deperacated for PIMS 
## with pims.open(movie_url) as stack:
##why remake the wheel.
#import tifffile
#import glob

# class VirtualStackAdapter:
# 	"""
# 	Stack Adapter Allows for the iteration through an image sequence that is typically too large to read fully into memmory.
# 	Corrently The Stack must have a strcuture of 13 z images and any number  of t positions. With an order of (zt)
# 	"""
# 	def __init__(self, path, z=13, reader= tifffile.imread):
# 		self.reader    = reader
# 		self.file_list = sorted(glob.glob(path+'/*.tif'))
# 		self.t         = -1
# 		self.z         = z

# 	def seek(self, frame):
# 		self.t = frame
# 		return(scipy.array( [self.reader(image_file) for image_file in self.file_list[self.t*self.z:self.t*self.z+self.z]]))
		
# 	def __iter__(self):
# 		return self._z_stack_gen_()

# 	def _z_stack_gen_(self):
# 		for i in range(0,len(self.file_list),self.z):
# 			yield scipy.array( [self.reader(image_file) for image_file in self.file_list[i:i+self.z]])
# class ImageJ:
# 	@staticmethod
# 	def stubb:
# 		pass