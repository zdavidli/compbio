#!/usr/local/bin/python
import sys
import operator
import math
from functools import reduce
#
# read data
#
data = open(sys.argv[1]).readlines()
#
#convert to numbers
#
x = []
for i in range(len(data)):
	x.append(float(data[i]))
#
# sum of absolute value of data
#
def abs_sum(numlist):
	abnum = list(map(abs, numlist))
	return reduce(operator.add, abnum, 0.0)

#
# sum of data
#
def sum(numlist):
	return reduce(operator.add, numlist, 0.0)
#
# mean of data
#
def mean(numlist):
	return sum(numlist)/len(numlist)
#
# standard deviation of data points
#
# first the sum of squares
#
def sum2(numlist):
	sum2 = 0
	for item in numlist:
		sum2 = sum2 + item*item
	return sum2
#
# the variance
#
def var(numlist):
	n = len(numlist)
	m = mean(numlist)
	dev = [0]*len(numlist)
	for i in range(len(numlist)):
		dev[i] = numlist[i] - m
	return sum2(dev)/float(n-1)
#
# the standard deviation
#
stdev = math.sqrt(var(x))
#
# Print output
#
print(('    N =%9i' %len(data)))
if ((sum(x) < 1.0) and (sum(x) > 0.0))| (sum(x) > 999999.9):
	print(('  Sum =%13.4e' %sum(x)))
else:
	print(('  Sum =%13.3f' %sum(x)))
 
if ((mean(x) < 1.0) and (mean(x) > 0.0)) | (mean(x) > 999999.9):
	print((' Mean =%13.4e' %mean(x)))
else:
	print((' Mean =%13.3f' %mean(x)))
 
if ((stdev < 1.0) and (stdev > 0.0)) | (stdev > 999999.9):
	print((' S.D. =%13.4e' %stdev))
else:
	print((' S.D. =%13.3f' %stdev))
 
if ((min(x) < 1.0) and (min(x) > 0.0))| (min(x) > 999999.9):
	print((' Min. =%13.4e' %min(x)))
else:
	print((' Min. =%13.3f' %min(x)))
 
if ((max(x) < 1.0) and (max(x) > 0.0))| (max(x) > 999999.9):
	print((' Max. =%13.4e' %max(x)))
else:
	print((' Max. =%13.3f' %max(x)))
 
if ((abs_sum(x) < 1.0) and (abs_sum(x) > 0.0)) | (abs_sum(x) > 999999.9):
	print(('|Sum| =%13.4e' %abs_sum(x)))
else:
	print(('|Sum| =%13.3f' %abs_sum(x)))

