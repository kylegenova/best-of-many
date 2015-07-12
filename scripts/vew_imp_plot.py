import sys
import numpy as np
import pylab as plt

if len(sys.argv) < 2:
    print "Usage: " + sys.argv[0] + " [file] [cutoff=None]"
    sys.exit()
f_name = sys.argv[1]

cutoff = None
if len(sys.argv) == 3:
    cutoff = int(sys.argv[2]) - 1

def min_l(l):
	return reduce(lambda a,b: min(a,b), l)

def max_l(l):
	return reduce(lambda a,b: max(a,b), l)

f = open(f_name, 'r')

inst_name = f_name.split(".")[0]

s = map(lambda x: x.split(" "), f.read().split("\n"))
s = filter(lambda x: x != [''], s)
iterations = map(lambda x: x[0], s)
objs = map(lambda x: x[1], s)
plt.xlabel("# Samples")
plt.ylabel("Best Tour")
plt.plot(iterations,objs, '-', linewidth=3)
if (not(cutoff == None)):
    plt.plot(cutoff+1, objs[int(cutoff)], 'rs', markersize=7)
plt.savefig(inst_name+"_cgp.png", dpi=200)
plt.show()

plt.clf()
#for i in range(0,3):#
#	plt.hist(sampled[i], bins=num_bins,label=samp_names[i])
#	plt.axis([overall_min,overall_max, 0, 250])
#	plt.title(inst_name)
    #	plt.xlabel("Tour Cost")
    #	plt.ylabel("Frequency")
    #    plt.legend(fontsize=10)
    #	plt.savefig(inst_name+samp_extension(samp_names[i])+".png", dpi=200)
#	plt.clf()
	#plt.show()