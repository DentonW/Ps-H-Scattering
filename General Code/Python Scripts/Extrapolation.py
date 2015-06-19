#!/usr/bin/python

#TODO: Add checks for whether files are good
#TODO: Make relative difference function

import sys, scipy, pylab
import numpy as np
from math import * 
import matplotlib.pyplot as plt
from xml.dom.minidom import parse, parseString
from xml.dom import minidom


def NumTermsOmega(omega):    # Return the number of terms for a given omega
	"""Uses combination with repetition to determine the number of terms for a given omega"""
	f = factorial
	k = 6
	omega = omega + 1
	n = f(omega+k-1) / (f(k) * f(omega-1))
	return int(n)
	

def FindTerms(FileName, NumTerms):
	"""Reorders the first NumTerms of the output of Todd program to find omega breakpoints"""
	f = open(FileName, 'r')

	# Get the value of omega
	Omega = int(f.readline().split()[1])
	print "Omega =", Omega

	# Skip these lines
	for i in range(3):
		f.readline()

	Terms = []
	for line in f:
		s = line.split()
		if len(s) == 0:
			break
		if s[0].isdigit():
			Terms.append(int(s[0]))
	f.close()

	print NumTerms, len(Terms)
	if NumTerms > len(Terms):
		print("Requesting more terms than are available in file...exiting.")
		exit()

	print "Number of terms in file", FileName, ": ", len(Terms)
	print "Number of terms to use:", str(NumTerms)
	print

	TermsSub = Terms[0:NumTerms]
	TermsSub.sort()

	# Create a list of numbers of terms for the full set for omega = 1 through Omega
	FoundTerms = []
	OmegaTerms = []
	for i in range(Omega+1):
		OmegaTerms.append(NumTermsOmega(i))

	for i in range(Omega+1):
		for j in range(len(TermsSub)):
			if TermsSub[j] == OmegaTerms[i]:
				print "Found", OmegaTerms[i], "at position", j+1
				FoundTerms = FoundTerms + [j+1]
				break
			if TermsSub[j] > OmegaTerms[i]:
				#print "Found next term past", OmegaTerms[i], "at position", j+1
				#FoundTerms = FoundTerms + [j+1]
				print "Found term before", OmegaTerms[i], "at position", j
				FoundTerms = FoundTerms + [j]
				break

	if TermsSub[len(TermsSub)-1] != OmegaTerms[Omega]:
		print "Last term at", len(TermsSub), "is less than", OmegaTerms[Omega]
		FoundTerms = FoundTerms + [len(TermsSub)]
				
	# Just here to put some extra space after running			
	print
	return FoundTerms


def Extrapolate(Phases, Omega, OmegaPower, LowerOmega):
	"""Fits the data to a straight line use SciPy's polyfit"""
	xdata = range(LowerOmega, Omega+1)
	xdata[:] = [x**OmegaPower for x in xdata]
	ydata = []
	
	for i in range(LowerOmega, Omega+1):
		ydata = ydata + [tan(Phases[i])]

	fit = scipy.polyfit(xdata, ydata, 1, None, True)
	polycoeffs = fit[0]
	residuals = fit[1][0]

	ExtrapData = [polycoeffs, residuals, xdata, ydata]
	return ExtrapData


def ExtrapolatePlot(Phases, Omega, OmegaPower, LowerOmega):
	"""Plots the fitted line for the extrapolation"""
	ExtrapData = Extrapolate(Phases, Omega, OmegaPower, LowerOmega)
	yfit = scipy.polyval(ExtrapData[0], ExtrapData[2])
	print yfit
	yfit = np.append(yfit, ExtrapData[0][1])
	print yfit
	p1 = plt.plot(ExtrapData[2], ExtrapData[3], 'k.')
	ExtrapData[2].append(0.0)
	p2 = plt.plot(ExtrapData[2], yfit, 'r-')
	print ExtrapData[2]
	print ExtrapData[3]
	plt.show()
	return


def ToBool(s):
    if (s.lower() == 'true'):
        return True
    return False


def ReadXMLData(xmldoc, tag):
    """ Helper function for ReadPhaseShifts """
    itemlist = xmldoc.getElementsByTagName(tag)
    data = []
    for s in itemlist:
        data.append(str(s.childNodes[0].nodeValue))
    if len(data) > 1:
        print "More than one set found for ", tag
    if data == []:
        return None
    return data[0]


def ReadPhaseShifts(Filename, FoundTerms, NumTests):
    """ Reads the complete list of phase shifts from a given phase file and returns a 2D array. """
    xmldoc = minidom.parse(Filename)  # Read the XML file
    
    shortfile = ReadXMLData(xmldoc, 'shortfile')
    longfile = ReadXMLData(xmldoc, 'longfile')
    energyfile = ReadXMLData(xmldoc, 'energyfile')
    lvalue = int(ReadXMLData(xmldoc, 'lvalue'))
    numterms = int(ReadXMLData(xmldoc, 'numterms'))
    numsets = int(ReadXMLData(xmldoc, 'numsets'))
    shielding = ReadXMLData(xmldoc, 'shielding')
    #if shielding == None:  # Not found in the input file
    #    shielding = 2*lvalue + 1  #@TODO: Is this a valid assumption?
    #else:
    #    shielding = int(shielding)

    if shielding != None:
        shielding = int(shielding)
    explambda = ReadXMLData(xmldoc, 'lambda')

    # Read in nonlinear parameters
    #@TODO: Handle multiple sets
    alpha = float(ReadXMLData(xmldoc, 'alpha'))
    beta  = float(ReadXMLData(xmldoc, 'beta'))
    gamma = float(ReadXMLData(xmldoc, 'gamma'))
    kappa = float(ReadXMLData(xmldoc, 'kappa'))
    mu = float(ReadXMLData(xmldoc, 'mu'))
    ordering = ReadXMLData(xmldoc, 'ordering')
    # Boolean values
    paired = ReadXMLData(xmldoc, 'paired')
    reorder = ReadXMLData(xmldoc, 'reorder')
    paired = ToBool(paired)
    reorder = ToBool(reorder)
    
    # Read in the phase shift data
    data = str(ReadXMLData(xmldoc, 'data'))
    data = data.split('\n')
    data = data[1:len(data)-1]  # First and last entries are blanks from the newlines
    
    if len(data) != numterms+1:  # Include the +1 for the 0th entry
        return None

    phases = []
    for n,d in enumerate(data):
        if n not in FoundTerms:
            continue
        line = d.split()
        if n != int(line[0]):
            print "Phase shift file indices do not match!"
            return None
        if len(line) != NumTests+1:
            print "Missing phase shift data on line " + str(n)
            return None
        line = [float(i) for i in line[1:]]
        phases.append(line)
    
    return phases


# def GetPhaseShifts(f, FoundTerms, TotalTerms, NumTests):
# 	"""Reads phase shifts at specified terms"""
# 	Omega = len(FoundTerms)-1
#
# 	for i in range(3):
# 		f.readline()
#
# 	PhaseShifts = range(NumTests)
# 	for i in range(NumTests):
# 		PhaseShifts[i] = []
# 	j = 0  # Corresponds to Omega = 0
#
# 	for i in range(1,FoundTerms[Omega]+1):  # Assuming that the last term is the highest for Omega.
# 		#@TODO: Check for end of file somehow?
# 		line = f.readline()
# 		if line[0] == '0':
# 			line = f.readline()
# 		s = line.split()
# 		if (len(s) == 0):
# 			print " "
# 			print "Error reading phase shifts: line length of 0"
# 			exit()
# 		if (len(s) < NumTests):
# 			print " "
# 			print "Error reading phase shifts: line length of " + str(len(s)) + " < " + str(NumTests)
# 			exit()
#
# 		if i == FoundTerms[j]:
# 			j = j + 1
# 			if j > Omega+1:
# 				print "Internal error reading phase shifts"  # This shouldn't happen.
# 				return []
# 			for k in range(NumTests):
# 				#PhaseShifts[k+1] = PhaseShifts[k+1] + [float(s[k+1])]
# 				PhaseShifts[k].append(float(s[k+1]))
#
# 	# Skip rest of terms if we are not using them all
# 	print "Skipping " + str(TotalTerms-FoundTerms[Omega]+1) + " terms"
# 	for i in range(1,TotalTerms-FoundTerms[Omega]+1):
# 		f.readline()
#
# 	return PhaseShifts


#
# Main function follows
#

# These are hardcoded right now, but we could probably write something to read them in later.

# 109 of these!  #@TODO: Could also just read from file and match up, but that will probably be difficult.
Headings = [ "Kohn", "Inverse Kohn", "Complex Kohn (S)", "Complex Kohn (T)", "Gen Kohn tau = 0.0", "Gen Kohn tau = 0.1", "Gen Kohn tau = 0.2", "Gen Kohn tau = 0.3",
            "Gen Kohn tau = 0.4", "Gen Kohn tau = 0.5", "Gen Kohn tau = 0.6", "Gen Kohn tau = 0.7", "Gen Kohn tau = pi/4", "Gen Kohn tau = 0.8", "Gen Kohn tau = 0.9",
            "Gen Kohn tau = 1.0", "Gen Kohn tau = 1.1", "Gen Kohn tau = 1.2", "Gen Kohn tau = 1.3", "Gen Kohn tau = 1.4", "Gen Kohn tau = 1.5", "Gen Kohn tau = pi/2",
            "Gen Kohn tau = 1.6", "Gen Kohn tau = 1.7", "Gen Kohn tau = 1.8", "Gen Kohn tau = 1.9", "Gen Kohn tau = 2.0", "Gen Kohn tau = 2.1", "Gen Kohn tau = 2.2",
            "Gen Kohn tau = 2.3", "Gen Kohn tau = 3*pi/4", "Gen Kohn tau = 2.4", "Gen Kohn tau = 2.5", "Gen Kohn tau = 2.6", "Gen Kohn tau = 2.7", "Gen Kohn tau = 2.8",
            "Gen Kohn tau = 2.9", "Gen Kohn tau = 3.0", "Gen Kohn tau = pi", "Gen T Kohn tau = 0.0", "Gen T Kohn tau = 0.1", "Gen T Kohn tau = 0.2", "Gen T Kohn tau = 0.3",
            "Gen T Kohn tau = 0.4", "Gen T Kohn tau = 0.5", "Gen T Kohn tau = 0.6", "Gen T Kohn tau = 0.7", "Gen T Kohn tau = pi/4", "Gen T Kohn tau = 0.8",
            "Gen T Kohn tau = 0.9", "Gen T Kohn tau = 1.0", "Gen T Kohn tau = 1.1", "Gen T Kohn tau = 1.2", "Gen T Kohn tau = 1.3", "Gen T Kohn tau = 1.4",
            "Gen T Kohn tau = 1.5", "Gen T Kohn tau = pi/2", "Gen T Kohn tau = 1.6", "Gen T Kohn tau = 1.7", "Gen T Kohn tau = 1.8", "Gen T Kohn tau = 1.9",
            "Gen T Kohn tau = 2.0", "Gen T Kohn tau = 2.1", "Gen T Kohn tau = 2.2", "Gen T Kohn tau = 2.3", "Gen T Kohn tau = 3*pi/4", "Gen T Kohn tau = 2.4",
            "Gen T Kohn tau = 2.5", "Gen T Kohn tau = 2.6", "Gen T Kohn tau = 2.7", "Gen T Kohn tau = 2.8", "Gen T Kohn tau = 2.9", "Gen T Kohn tau = 3.0",
            "Gen T Kohn tau = pi", "Gen S Kohn tau = 0.0", "Gen S Kohn tau = 0.1", "Gen S Kohn tau = 0.2", "Gen S Kohn tau = 0.3", "Gen S Kohn tau = 0.4",
            "Gen S Kohn tau = 0.5", "Gen S Kohn tau = 0.6", "Gen S Kohn tau = 0.7", "Gen S Kohn tau = pi/4", "Gen S Kohn tau = 0.8", "Gen S Kohn tau = 0.9",
            "Gen S Kohn tau = 1.0", "Gen S Kohn tau = 1.1", "Gen S Kohn tau = 1.2", "Gen S Kohn tau = 1.3", "Gen S Kohn tau = 1.4", "Gen S Kohn tau = 1.5",
            "Gen S Kohn tau = pi/2", "Gen S Kohn tau = 1.6", "Gen S Kohn tau = 1.7", "Gen S Kohn tau = 1.8", "Gen S Kohn tau = 1.9", "Gen S Kohn tau = 2.0",
            "Gen S Kohn tau = 2.1", "Gen S Kohn tau = 2.2", "Gen S Kohn tau = 2.3", "Gen S Kohn tau = 3*pi/4", "Gen S Kohn tau = 2.4", "Gen S Kohn tau = 2.5",
            "Gen S Kohn tau = 2.6", "Gen S Kohn tau = 2.7", "Gen S Kohn tau = 2.8", "Gen S Kohn tau = 2.9", "Gen S Kohn tau = 3.0", "Gen S Kohn tau = pi" ]
NumTests = 109  #@TODO: Could just calculate the length of Headings
# Headings = [ "Kohn", "Inverse Kohn", "Complex Kohn (S)", "Complex Kohn (T)", "Gen Kohn tau = 0.0", "Gen Kohn tau = 0.1", "Gen Kohn tau = 0.2", "Gen Kohn tau = 0.3",
#              "Gen Kohn tau = 0.4", "Gen Kohn tau = 0.5", "Gen Kohn tau = 0.6", "Gen Kohn tau = 0.7", "Gen Kohn tau = pi/4", "Gen Kohn tau = 0.8", "Gen Kohn tau = 0.9",
#              "Gen Kohn tau = 1.0", "Gen Kohn tau = 1.1", "Gen Kohn tau = 1.2", "Gen Kohn tau = 1.3", "Gen Kohn tau = 1.4", "Gen Kohn tau = 1.5", "Gen Kohn tau = pi/2",
#              "Gen Kohn tau = 1.6", "Gen Kohn tau = 1.7", "Gen Kohn tau = 1.8", "Gen Kohn tau = 1.9", "Gen Kohn tau = 2.0", "Gen Kohn tau = 2.1", "Gen Kohn tau = 2.2",
#              "Gen Kohn tau = 2.3", "Gen Kohn tau = 3*pi/4", "Gen Kohn tau = 2.4", "Gen Kohn tau = 2.5", "Gen Kohn tau = 2.6", "Gen Kohn tau = 2.7", "Gen Kohn tau = 2.8",
#              "Gen Kohn tau = 2.9", "Gen Kohn tau = 3.0", "Gen Kohn tau = pi" ]
# NumTests = 39

if len(sys.argv) < 6:
	print """Usage: Extrapolation.py <energyfile> <phasefile> <outputfile> <# of terms in file> <# of terms to use> <lower omega> <optional: upper omega>
Example: Extrapolation.py energy.txt phase.txt output.txt 1216 1216 3"""
	exit()

if sys.argv[4].isdigit() == False:
	print "Error: The fourth argument must be a number."
	exit()
if sys.argv[5].isdigit() == False:
	print "Error: The fifth argument must be a number."
	exit()
if sys.argv[6].isdigit() == False:
	print "Error: The sixth argument must be a number."
	exit()

FoundTerms = FindTerms(sys.argv[1], int(sys.argv[5]))
Omega = len(FoundTerms)-1
UpperOmega = Omega
LowerOmega = int(sys.argv[6])

if len(sys.argv) > 7:
	if sys.argv[7].isdigit() == False:
		print "Error: The seventh argument must be a number."
		exit()
	UpperOmega = int(sys.argv[7])
	if UpperOmega < LowerOmega or UpperOmega < 0:
		print "Error: Upper omega must be in the range " + str(LowerOmega) + "-" + str(Omega)
		exit()

if LowerOmega > UpperOmega:
	print "Error: Lower omega must be in the range 0-" + str(UpperOmega)
	exit()



print
g = open(sys.argv[3], 'w')

g.write("Results from " + sys.argv[1] + " and " + sys.argv[2] + "\n")
g.write(" with " + str(sys.argv[5]) + " terms and starting at omega = " + str(sys.argv[6]) + "\n\n")
g.write("Extrapolated values\n")
g.write("-------------------\n")

PhaseShiftLists = range(NumTests)
ExtrapolationLists = range(NumTests)
DList = range(NumTests)
for i in range(NumTests):
	PhaseShiftLists[i] = []
	ExtrapolationLists[i] = []
	DList[i] = []

PhaseShifts = np.array(ReadPhaseShifts(sys.argv[2], FoundTerms, NumTests))
#print PhaseShifts
#print len(PhaseShifts[0])
#exit()

# Iterate over the sets of tests
for j in range(NumTests):
	RMin = 1.0e5  # Just some very high value
	MinVal = 0
	Phases = PhaseShifts[:,j]
	# This loop iterates from d = -6 to -0.1 in increments of 0.01, testing the extrapolation fit by
	#  comparing the residuals.  The d that gives the smallest residuals is used, and the extrapolation
	#  is saved.
	for i in range(0,690):
		Residuals = Extrapolate(Phases, UpperOmega, -7.0+i/100.0, LowerOmega)[1]
		if Residuals < RMin:
			RMin = Residuals
			MinVal = i
	print
	print "Results for " + Headings[j] + ":"
	print "Smallest residuals at", -7.0+MinVal/100.0, "of", RMin
		
	DList[j] = -7.0+MinVal/100.0
	PhaseShiftLists[j] = Phases
	Extrapolation = Extrapolate(Phases, UpperOmega, -7.0+MinVal/100.0, LowerOmega)
	ExtrapolationLists[j] = Extrapolation
	print "Extrapolated value =", atan(Extrapolation[0][1])
	print "Relative difference % =", abs((atan(Extrapolation[0][1]) - Phases[np.size(Phases)-1]) / (atan(Extrapolation[0][1]) + Phases[np.size(Phases)-1]) * 2) * 100
		
	print "Coefficients: ", Extrapolation[0]
		
	Line = Headings[j] + ": " + str(atan(Extrapolation[0][1])) + "\n"
	g.write(Line)
		
	print "w3 - w4: " + str(abs(Phases[3] - Phases[4]))
	if UpperOmega >= 5:
		print "w4 - w5: " + str(abs(Phases[4] - Phases[5]))
	if UpperOmega >= 6:
		print "w5 - w6: " + str(abs(Phases[5] - Phases[6]))
	if UpperOmega >= 7:
		print "w6 - w7: " + str(abs(Phases[6] - Phases[7]))


g.write("\n")
g.write("\n")
g.write("\n")
g.write("More detailed analysis\n")
g.write("----------------------\n")

g.write("\n")
g.write("Reordered terms:\n")
for i in range(len(FoundTerms)):
	g.write("Found " + str(NumTermsOmega(i)) + " at position " + str(FoundTerms[i]) + "\n")
g.write("\n")

for i in range(NumTests):
	g.write("\nResults for " + Headings[i] + ":\n")
	g.write("Phase shifts: ")
	for j in range(len(PhaseShiftLists[i])):
		g.write(str(PhaseShiftLists[i][j]) + " ")
	g.write("\n")

	g.write("Phase shift differences in omega: ")
	for j in range(len(PhaseShiftLists[i]) - 1):
		g.write(str(abs(PhaseShiftLists[i][j] - PhaseShiftLists[i][j+1])) + " ")
	g.write("\n")

	g.write("Phase shift difference ratios: ")
	for j in range(len(PhaseShiftLists[i]) - 2):
		#print PhaseShiftLists[i][j], PhaseShiftLists[i][j+1], PhaseShiftLists[i][j+2]
		g.write(str(abs( (PhaseShiftLists[i][j+1] - PhaseShiftLists[i][j+2]) / (PhaseShiftLists[i][j] - PhaseShiftLists[i][j+1]) )) + " ")
	g.write("\n")

	for j in range(LowerOmega+1,UpperOmega):
		if abs(PhaseShiftLists[i][j] - PhaseShiftLists[i][j+1]) > abs(PhaseShiftLists[i][j-1] - PhaseShiftLists[i][j]):
			g.write("No convergence pattern exists.\n")
	g.write("Smallest residuals at d = " + str(DList[i]) + " of " + str(ExtrapolationLists[i][1]) + "\n")
	g.write("Coefficients of " + str(ExtrapolationLists[i][0]) + "\n")
	reldiff = abs((atan(ExtrapolationLists[i][0][1]) - PhaseShiftLists[i][len(PhaseShiftLists[i])-1]) / (atan(ExtrapolationLists[i][0][1]) + PhaseShiftLists[i][len(PhaseShiftLists[i])-1]) * 2) * 100
	g.write("Relative difference % = " + str(reldiff) + "\n")
	g.write("Extrapolated value = " + str(atan(ExtrapolationLists[i][0][1])) + "\n")
	# This can be re-enabled to look at the fit lines with the phase shifts.
	#if i == 3:  # S-matrix
	#	ExtrapolatePlot(PhaseShiftLists[i], Omega, DList[i], LowerOmega)

g.close()


exit()
