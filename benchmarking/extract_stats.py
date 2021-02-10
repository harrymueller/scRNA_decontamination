# libs
import pandas as pd
import numpy as np
import datetime

# paths
pid_path = "/home/harry/gdrive/Education/Uni/GCRL2000/new_results/benchmarking/none_soupx/pid.txt"
usage_path = "/home/harry/gdrive/Education/Uni/GCRL2000/new_results/benchmarking/none_soupx/usage.txt"
output_path = "/home/harry/gdrive/Education/Uni/GCRL2000/new_results/benchmarking/usage.xlsx"

num_cells = {"BG1_BG20C": 2356, "BG3_BG21C": 4084, "BG5_BG22C": 5411, "BG52_BG20C_MeOH": 1886, "BG54_BG21C_MeOH": 4672, "BG56_BG22C_MeOH": 5069}



####################################
# FUNCTIONS: used for .apply(...)
####################################
# function to return ram in GB
def format_rss(x):
	if (x[-1] == "G"):
		return np.double(x[:-1])
	else:
		return np.double(x[:-1])/1000

# function to convert timedelta to nicely formatted string
def timedelta_to_str(x):
	seconds = int(x / np.timedelta64(1, 's'))
	return str(datetime.timedelta(seconds = round(seconds)))
	
	
	
####################################
# PIDS
####################################
pids = pd.read_table(pid_path, sep=" ")

# adding new columns to pid for statistics calculated later
empty = np.empty((len(pids), 6))
empty[:] = np.nan

pids[["num_cells", "time", "cpu_max", "cpu_avg", "ram_max", "ram_avg"]] = empty

pids["PID"] = pids["PID"].astype(np.int32)

for sample_id in num_cells:
	pids.loc[pids["SAMPLE_ID"] == sample_id, "num_cells"] = num_cells[sample_id]

pids["num_cells"] = pids["num_cells"].astype(np.int32)

	
####################################
# USAGE
####################################
# header
usage = [["Time", "UID", "PID", "%usr", "%system", "%guest", "%wait", "%CPU", "CPU", "minflt/s", "majflt/s", "VSZ", "RSS", "%MEM", "Command"]]

# reading in from file
with open(usage_path, "r") as f:
	f.readline()
	lines = f.readlines()

	for l in lines:
		if (l != "\n" and l[0] != "#"):
			usage.append(l.replace("\n", "").split())

# converting to pandas df and selecting only important columns
usage = pd.DataFrame(usage[1:], columns=usage[0])
usage = usage[["Time", "PID", "%CPU", "RSS"]]

# casting
usage["Time"] = [pd.Timedelta(x) for x in usage["Time"]]
usage["PID"] = usage["PID"].astype(np.int32)
usage["%CPU"] = usage["%CPU"].apply(lambda x: x[:-1]).astype(np.double)
usage["RSS"] = usage["RSS"].apply(lambda x: format_rss(x)).astype(np.double)



####################################
# Adding statistics to pid df
####################################
for i in pids["PID"]:
	# selects on this PID and removes the first 3 recordings (disrupted)
	subset = usage.loc[usage["PID"] == i,].iloc[3:]

	# removes any entries that are 3 std above the mean for either %CPU or RSS
	mask1 = ((subset["%CPU"] - np.mean(subset["%CPU"])) < 2 * np.std(subset["%CPU"])) 
	mask2 = ((subset["RSS"] - np.mean(subset["RSS"])) < 2 * np.std(subset["RSS"]))
	
	subset = subset[mask1 & mask2]

	# adding stats to pid df
	pids.loc[pids["PID"] == i, "time"] = subset.iloc[-1]["Time"] - subset.iloc[0]["Time"]
	
	pids.loc[pids["PID"] == i, "cpu_max"] = max(subset["%CPU"][mask1])
	pids.loc[pids["PID"] == i, "cpu_avg"] = np.mean(subset["%CPU"])
	
	pids.loc[pids["PID"] == i, "ram_max"] = max(subset["RSS"])
	pids.loc[pids["PID"] == i, "ram_avg"] = np.mean(subset["RSS"])
	
	#if (abs(float(pids.loc[pids["PID"] == i, "cpu_avg"]) - 100) < 3):
	#	pids.loc[pids["PID"] == i, "cpu_avg"] = 100
	
	
	
####################################
# adding tallies for each method
####################################
tallies = pd.DataFrame(columns = pids.columns[1:])

for (index, method) in enumerate(np.unique(pids["METHOD"])):
	mask = pids["METHOD"] == method
	
	# mean and std
	for stat in (np.mean, np.std):
		# template dictionary
		new = { i: None for i in tallies.columns }
		new["METHOD"] = method
		new["SAMPLE_ID"] = "Mean" if (stat == np.mean) else "Std"
	
		# looping through cols
		for i in pids.columns[3:]:
			if (i == "time"): # time
				seconds = (pids.loc[mask, i] / np.timedelta64(1, 's')).astype(int)
				new[i] = str(datetime.timedelta(seconds = round(stat(seconds))))
			elif (i == "cpu_max" or i == "ram_max"): #max
				new[i] = max(pids.loc[mask,i]) if (stat == np.mean) else np.std(pids.loc[mask,i])
			else: # otherwise
				new[i] = stat(pids.loc[mask, i])
	
		# turning in df
		new = pd.DataFrame(new, index=[2*index + (0 if (stat == np.mean) else 1)])
		
		# adding to tallies
		tallies = pd.concat([tallies, new])#.reset_index(drop=True)

		

####################################
# saving df to file
####################################
# converting timedelta to readable string
pids["time"] = pids["time"].apply(lambda x: timedelta_to_str(x)).astype(str)


with pd.ExcelWriter(output_path) as writer:
	tallies.to_excel(writer, sheet_name="Method Tallies", index=False)  
	pids.to_excel(writer, sheet_name="Methods and Sample IDs", index=False)  


print(tallies)
print(pids)
	
	
