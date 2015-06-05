#!/usr/bin/python

#TODO: Add checks for whether files are 

import sys, scipy, numpy, matplotlib, pylab
from math import * 


def NumTermsOmega(omega):    # Return the number of terms for a given omega
	"""Uses combination with repetition to determine the number of terms for a given omega"""
	f = factorial
	k = 6
	omega = omega + 1
	n = f(omega+k-1) / (f(k) * f(omega-1))
	return int(n)
	

def FindTerms(Filename, NumTerms):
	"""Reorders the first NumTerms of the output of Todd program to find omega breakpoints"""
	f = open(Filename, 'r')

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

	if NumTerms > len(Terms):
		print("Requesting more terms than are available in file...exiting.")
		exit()

	print "Number of terms in file", Filename, ": ", len(Terms)
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
				print i, ": Found", OmegaTerms[i], "at position", j+1
				FoundTerms = FoundTerms + [j+1]
				break
			if TermsSub[j] > OmegaTerms[i]:
				print i, ": Found next term past", OmegaTerms[i], "at position", j+1
				FoundTerms = FoundTerms + [j+1]
				break

	if TermsSub[len(TermsSub)-1] != OmegaTerms[Omega]:
		print Omega, ": Last term at", len(TermsSub), "is less than", OmegaTerms[Omega]
		FoundTerms = FoundTerms + [len(TermsSub)]
				
	# Just here to put some extra space after running			
	print
	return FoundTerms


#
# Main function follows
#

if len(sys.argv) < 3:
	print """Usage: ViewOmegaCutoffs.py <energyfile> <# of terms to use>
Example: ViewOmegaCutoffs.py energy.txt 1216"""
	exit()

if sys.argv[2].isdigit() == False:
	print "Error: The second argument must be a number."
	exit()

FoundTerms = FindTerms(sys.argv[1], int(sys.argv[2]))

exit()

