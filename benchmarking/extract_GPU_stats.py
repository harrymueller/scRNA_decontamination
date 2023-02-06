# NOTE: Mean max usage -> max

# libs
import pandas as pd
import numpy as np
import datetime

# paths
DIR = "/data/Perkins/benchmarking"
USAGE_PATH = DIR + "/cellbender_gpu/gpu_usage.csv"
OUTPUT_PATH = DIR + "/cellbender_gpu_gpu_usage.xlsx"
NUM_CELLS = {"BG1_BG20C": 2356, "BG3_BG21C": 4084, "BG5_BG22C": 5411, "BG52_BG20C_MeOH": 1886, "BG54_BG21C_MeOH": 4672, "BG56_BG22C_MeOH": 5069}



####################################
# USAGE
####################################
def get_usage(usage_path):
	baseline_usage = []
	usage = [["sample_id", "gpu_util", "mem_util", "mem_total", "mem_free", "mem_used"]]

	# reading in from file
	with open(usage_path, "r") as f:
		lines = f.readlines()

		for l in lines:
			if (l[0] == "u" or l == "\n"): continue # header row
			elif (l[0] == "c"): current_sample = l.split(" ")[-1][:-2]
			else:
				l_split = l.replace("\n", "").replace(" ", "").replace("MiB", "").replace("%", "").split(",")

				# gpu not begining used req. gpu at 0% usage and < 500MiB RAM
				if (l[0] == "0" and int(l_split[-1]) < 500): baseline_usage.append(int(l_split[-1])) # add stats for baseline
				else: usage.append([current_sample, *l_split])
		
	# converting to pandas df and selecting only important columns
	usage = pd.DataFrame(usage[1:], columns=usage[0])
	mem_used = round(sum(baseline_usage) / len(baseline_usage))

	# converting cols to correct type
	for col in ["gpu_util", "mem_total"]:
		usage[col] = usage[col].astype(int)	
	for col in ["mem_util", "mem_free", "mem_used"]:
		 usage[col] = usage[col].astype(np.double)	

	# calibrating mem_used based on background GPU usage
	usage["mem_used"] = usage["mem_used"].apply(lambda x: x - mem_used if (x - mem_used > 0) else 0)

	# recalculating mem_util based on calibrated value
	mem_total = usage["mem_total"].iloc[0]
	usage["mem_util"] = usage["mem_used"].apply(lambda x: x / mem_total * 100)

	# convert MiB to GiB
	for col in ["mem_free", "mem_used"]:
		 usage[col] = usage[col].apply(lambda x: x / 1000)	

	return usage


####################################
# Adding statistics to pid df
####################################
def summarise(usage):
	new_entries = []

	for (index, id) in enumerate(NUM_CELLS.keys()):
		new_entries.append(get_mean_for_sample(usage, id, index*2))
		
	new_entries.append(get_mean_for_sample(usage, "OVERALL", (index + 1)*2))

	summary = pd.concat([*new_entries])
	return summary


####################################
# Small mean function
####################################
def get_mean_for_sample(usage, id, i):
	rows = []
	for stat in (np.mean, np.std):
		new = {"sample_id": id, "statistic": "Mean" if (stat == np.mean) else "Std"}
		
		if id != "OVERALL":
			# selects on this PID
			mask = usage["sample_id"] == id[:-1] # last letter gets cutoff ?

		# averaging
		for col in ["gpu_util", "mem_util", "mem_used"]:
			if id != "OVERALL": #could be better
				new[col] = stat(usage[mask][col])
				new["max_" + col] = np.max(usage[mask][col]) if (stat == np.mean) else np.NaN
			else:
				new[col] = stat(usage[col])
				new["max_" + col] = np.max(usage[col]) if (stat == np.mean) else np.NaN

		rows.append(new)

	return pd.DataFrame(rows, index=[i, i+1])


####################################
# MAIN
####################################
def main():
	# input data
	usage = get_usage(USAGE_PATH)

	summary = summarise(usage)

	print(summary)

	with pd.ExcelWriter(OUTPUT_PATH) as writer:
		summary.to_excel(writer, sheet_name="GPU Summaries", index=False)  

if __name__ == "__main__":
	main()