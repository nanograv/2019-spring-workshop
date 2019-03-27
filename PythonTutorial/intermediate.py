################### ANSWERS FOR THE INTERMEDIATE EXERCISE ######################
#import the packages we will need for this exercise 
import numpy as np
import matplotlib.pyplot as plt

#print the header for our table:
print 'SNR\tRMS(std)\tRMS(calc)' #the \t represents a tab in the string

#read in the file
f =open('intermediate.txt','r')

###################### LOOP OVER EACH LINE IN THE FILE #########################
#Note that each line contains one gaussian time series
for line in f.readlines():

################### CONVERT THE DATA TO A NUMPY ARRAY ##########################
	#strip the new line character (\n) off the end of the line, then split
	#the line into a list of entries
	data = line.strip().split(',')
	#convert the list into a numpy array (and convert all the entries from strings to floats)
	dataarray = np.array(data, dtype='float')
	
################################ MAKE A PLOT ###################################
	#determine the number of entries in your array, and then plot your data,
	#one data point per entry.
	Nentries = len(dataarray) 
	plt.plot(np.arange(0,Nentries,1),dataarray,'.') 
	#Note: arange generates a numpy array with syntax (first entry,length,step).
	#Further note: the '.' makes the it plot points - you can also use plt.scatter()
	#To display the plot to the screen we would use plt.show() here. However,
	#We're going to plot all the data on top of each other and show at the end,
	#so we don't use plt.show() right now 

####################### CALCULATE SIGNAL TO NOISE RATIO ########################	
	#find the peak flux
	peakflux = np.max(dataarray)

	#here are two different ways of calculating the RMS.
	#Note that we cheat here: we know the pulse is in the first half of the
	#time series, so we just calculate the rms noise on the second half
	rmsflux = np.std(dataarray[Nentries/2:])
	rmsflux2 = np.sqrt(np.mean(dataarray[Nentries/2:]**2))

	#print out your signal to noise calculation
	print '%(SN).2f\t%(peak).2f\t%(rms).2f\t%(rms2).2f' %{'SN':peakflux/rmsflux, 'peak':peakflux, 'rms':rmsflux, 'rms2':rmsflux2}

################################### END LOOP ###################################

############################# DISPLAY THE PLOT #################################
plt.xlabel('time (s)')
plt.ylabel('Flux (mJy)')
plt.title("Michael's Face When:")
plt.show()

