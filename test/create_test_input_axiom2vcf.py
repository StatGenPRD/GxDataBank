#!/GWD/appbase/tools/bin/python2.7

import random
from subprocess import Popen, PIPE

def select_data(file_in,file_out,criteria):
	"""Output data from FILE_IN where sample name contains CRITERIA to FILE_OUT"""

	indices = list()
	f = Popen(['zcat', file_in], stdout=PIPE)
	with open(file_out,"w") as o:
		for line in f.stdout:
			if line.startswith("#"):
				o.write(line)
			else:
				cols = line.rstrip().split("\t")
				probeset_id = cols[0]
				data = cols[1:]

				if probeset_id == "probeset_id":
					indices = [i for i,header in enumerate(data) if header.find(criteria) > -1]
				
				selected = [data[i] for i in indices]
				o.write(probeset_id)
				o.write("\t")
				o.write("\t".join(selected))
				o.write("\n")

def sample_count(file_in):
	"""Return sample count from FILE_IN"""

	count = -1

	with open(file_in, "r") as f:
		for line in f:
			if line.startswith("probeset_id"):
				cols = line.rstrip().split("\t")
				count = len(cols) - 1
				break
	
	return count

def adjust_performance_values(file_in,file_out,n_old,n_new):
	"""Output performance values from FILE_IN to FILE_OUT, adjusted from N_OLD to N_NEW"""

	nMinorAllele = 7
	n_AA = 9
	n_AB = 10
	n_BB = 11
	n_NC = 12
	CR = 2

	f = Popen(['zcat', file_in], stdout=PIPE)
	with open(file_out,"w") as o:
		
		for line in f.stdout:
			if line.startswith("probeset_id"):
				o.write(line)
				continue

			cols = line.rstrip().split("\t")

			# Determine difference to account for
			difference = n_old - n_new
			accounted = difference

			indices = [n_NC,n_BB,n_AB,n_AA]

			# First make count columns numeric
			for i in indices:
				cols[i] = int(cols[i])

			# Adjust counts
			while accounted > 0:
				for index in indices:
					if cols[index] > 1:
						i = random.randint(0,min(accounted, cols[index]))
						i = min(i, cols[index] - 1)
						cols[index] = cols[index] - i
						accounted = accounted - i

			# Call Rate
			if cols[n_NC] == 0:
				cols[CR] = '100'
			else:
				n_Passing = (n_new - cols[n_NC]) * 1.0
				n_Total = n_new * 1.0
				cols[CR] = "%.13f" % (n_Passing / n_Total)

			# nMinorAllele
			cols[nMinorAllele] = cols[n_AB] + 2*cols[n_BB]

			cols = [str(x) for x in cols]
			o.write("\t".join(cols))
			o.write("\n")

CRITERIA = "_NA1852"
DIR = "/GWD/appbase/projects/RD-MDD-GX/QC0_BPS/BSA_2014-09-05_ClusterSet1/Axiom/BSA_2014-09-05_ClusterSet1/"

calls_in = DIR + "AxiomGT1.calls.txt.gz"
calls_out = "./AxiomGT1.calls.txt"

confs_in = DIR + "AxiomGT1.confidences.txt.gz"
confs_out = "./AxiomGT1.confidences.txt"

perf_in = DIR + "Ps.performance.txt.gz"
perf_out = "./AxiomGT1.performance.txt"

select_data(calls_in, calls_out, CRITERIA)
select_data(confs_in, confs_out, CRITERIA)

old_n = sample_count(calls_in)
new_n = sample_count(calls_out)

print "Old N = %d; New N = %d;" % (old_n, new_n)
adjust_performance_values(perf_in,perf_out,old_n,new_n)
