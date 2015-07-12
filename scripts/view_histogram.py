import sys
import numpy as np
import pylab as plt

if len(sys.argv) < 2:
    print "Usage: " + sys.argv[0] + " [file] [bins=15]"
    sys.exit()
f_name = sys.argv[1]
num_bins = 15

def min_l(l):
	return reduce(lambda a,b: min(a,b), l)

def max_l(l):
	return reduce(lambda a,b: max(a,b), l)

if (len(sys.argv) == 3):
	num_bins = int(sys.argv[2])
print num_bins
f = open(f_name, 'r')

inst_name = f_name.split(".")[0]
methods = []
results = []
s = filter(lambda x : x != '' and x != '\r', f.read().split("\n"))
f.close()

count = 0
for line in s:
    l = filter(lambda x : x != '' and x != '\r', line.split(","))
    methods += [l[0]]
    results += [[]]
    for elt in l[1:]:
        results[count] += [int(elt)]
    count += 1

def samp_extension(s):
	if ('Column Generation + swapRound' in s):
		return '_cg_sr'
	if ('Max Entropy' in s):
		return '_me'
	if ('Column Generation' in s):
		return '_cg'
	if ('Splitting Off + swapRound' in s):
		return '_so_sr'
	if ('Christofides' in s):
		return '_cfds'
	if ('Splitting Off' in s):
		return '_so'
	raise Exception

sampled=[results[2],results[3],results[5]]
samp_names=[methods[2],methods[3],methods[5]]
print samp_names
overall_min = min_l(results[2]+results[3]+results[5])
overall_max = max_l(results[2]+results[3]+results[5])

for i in range(0,3):
	plt.hist(sampled[i], bins=num_bins,label=samp_names[i])
	plt.axis([overall_min,overall_max, 0, 250])
	plt.title(inst_name)
	plt.xlabel("Tour Cost")
	plt.ylabel("Frequency")
	plt.legend(fontsize=10)
	plt.savefig(inst_name+samp_extension(samp_names[i])+".png", dpi=200)
	plt.clf()
	#plt.show()