import matplotlib.pyplot as plt
import numpy as np

species = (0,0.2,	0.4,	0.6,	0.8,	1,	1.5,	2)
values = {
    "Mirror": (82,	5,	5,	10,	4,	12,	2,	1),
    'preopt': (29,	0,	6,	10,	11,	41,	21,	3),
    'NT2': (32,	6,	10,	23,	18,	22,	5,	2),
    'Neb-TS': (65,	6,	8,	7,	7,	13,	11,	2),
}

x = np.arange(len(species))  # the label locations
width = 1/(len(values)+1)  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in values.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=2,rotation=90)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Number of ')
ax.set_title('Penguin attributes by species')
ax.set_xticks(x - width, species)
ax.legend(loc='upper right', ncols=4)
ax.set_ylim(0, 90)

plt.savefig(f"pictures/chart", dpi=300)