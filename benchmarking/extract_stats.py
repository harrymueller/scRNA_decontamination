# NOTE: Mean max usage -> max

# libs
import pandas as pd
import numpy as np
import datetime

# paths
pid_path = "/home/harry/gdrive/Education/Uni/GCRL2000/new_results/benchmarking/all/pid.txt"
usage_path = "/home/harry/gdrive/Education/Uni/GCRL2000/new_results/benchmarking/all/usage.txt"
output_path = "/home/harry/gdrive/Education/Uni/GCRL2000/new_results/benchmarking/usage_all.xlsx"

num_cells = {"BG1_BG20C": 2356, "BG3_BG21C": 4084, "BG5_BG22C": 5411, "BG52_BG20C_MeOH": 1886, "BG54_BG21C_MeOH": 4672, "BG56_BG22C_MeOH": 5069}

order = ["no_decontamination", "soupx:autoEstCont", "soupx:background_genes", "soupx:top_background_genes", "decontx:no_cell_types", "decontx:with_cell_types"]


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

# removing rows w/ names
pids = pids.loc[pids["PID"] != "PID"].reset_index()

# adding new columns to pid for statistics calculated later
empty = np.empty((len(pids), 6))
empty[:] = np.nan

pids[["num_cells", "Time", "cpu_max", "cpu_avg", "ram_max", "ram_avg"]] = empty

pids["PID"] = pids["PID"].astype(np.int32)

for sample_id in num_cells:
	pids.loc[pids["SAMPLE_ID"] == sample_id, "num_cells"] = num_cells[sample_id]

pids["num_cells"] = pids["num_cells"].astype(np.int32)
pids = pids[pids.columns[1:]]

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
		if (l != "\n" and l[0] != "#" and l[0] != "L"):
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
	# removes first 3 entries <- too variable and inaccurate
	usage = usage.drop(index = usage.loc[usage["PID"] == i].iloc[0:3].index).reset_index(drop = True)
	
	# selects on this PID
	mask = usage["PID"] == i
	
	# removes any entries that are 3 std above the mean for either %CPU or RSS
	mask1 = (abs(usage[mask]["%CPU"] - np.mean(usage[mask]["%CPU"])) > 3 * np.std(usage[mask]["%CPU"])) 
	mask2 = (abs(usage[mask]["RSS"] - np.mean(usage[mask]["RSS"])) > 3 * np.std(usage[mask]["RSS"]))
	
	usage = usage.drop(index = usage[mask][mask1 | mask2].index).reset_index(drop=True) # dropping outliers from usage df
	mask = usage["PID"] == i # resets mask

	# adding stats to pid df
	pids.loc[pids["PID"] == i, "Time"] = usage[mask].iloc[-1]["Time"] - usage[mask].iloc[0]["Time"]
	
	pids.loc[pids["PID"] == i, "cpu_max"] = max(usage[mask]["%CPU"])
	pids.loc[pids["PID"] == i, "cpu_avg"] = np.mean(usage[mask]["%CPU"])
	
	pids.loc[pids["PID"] == i, "ram_max"] = max(usage[mask]["RSS"])
	pids.loc[pids["PID"] == i, "ram_avg"] = np.mean(usage[mask]["RSS"])
	
	#if (abs(float(pids.loc[pids["PID"] == i, "cpu_avg"]) - 100) < 3):
	#	pids.loc[pids["PID"] == i, "cpu_avg"] = 100
				
####################################
# adding tallies for each method
####################################
tallies = pd.DataFrame(columns = list(pids.columns[1:3]) + list(pids.columns[4:]))
print(pids)
for (index, method) in enumerate(order):
	mask = usage["PID"].isin(pids[pids["METHOD"] == method]["PID"])
	
	# mean and std
	for stat in (np.mean, np.std):
		# template dictionary
		new = { i: None for i in tallies.columns }
		new["METHOD"] = method
		new["SAMPLE_ID"] = "Mean" if (stat == np.mean) else "Std"
	
		# looping through cols
		for i in pids.columns[4:]:
			if (i == "Time"): # Time
				seconds = (pids.loc[pids["METHOD"] == method, "Time"] / np.timedelta64(1, 's')).astype(int)
				new[i] = str(datetime.timedelta(seconds = round(stat(seconds))))
			elif (i == "cpu_max" or i == "ram_max"): #max
				k = "%CPU" if (i == "cpu_max") else "RSS"
				new[i] = max(usage.loc[mask,k]) if (stat == np.mean) else np.std(usage.loc[mask,k])
			else: # otherwise
				k = "%CPU" if (i == "cpu_avg") else "RSS"
				new[i] = stat(usage.loc[mask, k])
	
		# turning in df
		new = pd.DataFrame(new, index=[2*index + (0 if (stat == np.mean) else 1)])
		
		# adding to tallies
		tallies = pd.concat([tallies, new])#.reset_index(drop=True)

tallies.columns = ["Method", "Statistic", "Time Taken", "Max CPU Usage", "Avg CPU Usage", "Max RAM Usage (GB)", "Avg RAM Usage (GB)"]
print(tallies)

####################################
# saving df to file
####################################
# converting timedelta to readable string
pids["Time"] = pids["Time"].apply(lambda x: timedelta_to_str(x)).astype(str)


with pd.ExcelWriter(output_path) as writer:
	tallies.to_excel(writer, sheet_name="Method Tallies", index=False)  
	pids.to_excel(writer, sheet_name="Methods and Sample IDs", index=False)  

	
