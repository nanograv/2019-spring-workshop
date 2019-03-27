##################### ANSWERS FOR THE BEGINNER EXERCISE ########################
#do you want to print intermediate steps?
verbose=True

############### READING IN THE DATA FILE INTO A LIST OF INTEGERS ###############
### read in the text file
f = open('beginner.txt','r')

### store the text file into a list, one entry per line
values = f.readlines()

#if you've set verbose=True (or verbose=1.), print out the list
if verbose:
	print "My list of lines:", values

### turn each value in the list from a string to an integer
	#Note the syntax of my for loop
	#len(values) gives the length of the list called values
	#range returns a list of integers from 0 to len(values)-1
for i in range(len(values)): 
	#values[i] gives the i-th element of the list values.
	#values[i] is then a string object which ends in a new line, '\n', so
	#we use the .strip() method to remove the new line symbol.
	#we then use the int() function to convert the remaining string to an
	#int, and reassign the ith element of values with this integer 
	values[i] = int(values[i].strip())

if verbose: 
	print "My list of values (as ints):", values

######################### SORTING THE LIST ##################################### 

### make a sorted list that we will fill with values
sortedlist = []

### loop through our list of integer values, using similar syntax to above 
for i in range(len(values)):
	currentmin = values[0] #store the 1st value in the list as our current minimum

	#loop through the remaining values in the list. Note the alternative syntax for the for loop.
	for value in values[1:]: #the [1:] takes all the elements except the 0th element of values
		if value < currentmin: #check which thing is lower, the value or our current minimum
			currentmin = value #if the value is lower, assign that to be the current min.

	#add the minimum value we found to the sorted list we created earlier
	sortedlist.append(currentmin)
	#remove that value from the values list
	values.remove(currentmin)

### Yay, we've created a sorted list. Let's have a look at it.
print "My sorted list:",sortedlist		
		
