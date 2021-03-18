import pandas

usage_path = "usage.txt"
output_path = "out.txt"

# header
usage = [["Time", "UID", "PID", "minflt/s", "majflt/s", "VSZ", "RSS", "%MEM", "Command"]]

# reading in from file
with open(usage_path, "r") as f:
	f.readline()
	lines = f.readlines()

	for l in lines:
		if (l != "\n" and l[0] != "#" and l[0] != "L"):
			usage.append(l.replace("\n", "").split())

# converting to pandas df
usage = pd.DataFrame(usage[1:], columns=usage[0])
usage = usage[usage["Command"] == "R"]

#saving to csv
usage.to_csv(output_path)