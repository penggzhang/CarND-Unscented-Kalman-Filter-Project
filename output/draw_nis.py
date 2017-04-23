import matplotlib.pyplot as plt

nis_radar, nis_lidar = [], []

with open("output.txt", "r") as f:
	for line in f:
		items = line.split()
		if items[6] == "radar":
			nis_radar.append(float(items[7]))
		elif items[6] == "lidar":
			nis_lidar.append(float(items[7]))

nis = [nis_radar, nis_lidar]
title = ["NIS_radar", "NIS_lidar"]
chi_squared_value = [7.815, 5.991]

fig, axes = plt.subplots(2, 1)
for i, ax in enumerate(axes):
	ax.plot(nis[i], color="blue")
	ax.plot([0, len(nis[i])], [chi_squared_value[i], chi_squared_value[i]], color="green")
	ax.set_title(title[i])
#plt.show()
plt.savefig("NIS.png")
