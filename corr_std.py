#!/usr/bin/env python3
import sys, os, glob, re
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 1:
    sys.exit('usage: corr_std.py <files>...')

names = ['Average']
all_mppkis = []

dirs = [dir for path in sys.argv[1:] for dir in glob.glob(path) or (path,)]
for i, dir in enumerate(dirs):
    mppkis = []
    for f in os.listdir(dir):
        if i == 0:
            names.append(os.path.splitext(f)[0])
        with open(os.path.join(dir, f)) as f:
            txt = f.read()
        m = re.search(r'Final Score Run1_Conditional_MPPKI:\s*(.*)', txt)
        if m:
            mppkis.append(float(m.group(1)))
    all_mppkis.append(mppkis)

all_mppkis = np.array(all_mppkis)
avgs = np.mean(all_mppkis, 1)
all_with_avg = np.column_stack((avgs, all_mppkis))

corrs = np.array([np.corrcoef(avgs, col)[0,1] for col in all_with_avg.T])
#corrs = np.maximum(corrs, 0)
stds = np.std(all_with_avg, 0)

#################

width = .3
n = len(corrs)//2
ind = np.arange(n)

ax1 = plt.subplot()
bars1 = ax1.bar(ind-width/2, corrs[:n], width, color='r')
ax1.set_xticks(ind)
ax1.set_xticklabels(names[:n])
ax1.set_xlabel('Traces')
ax1.set_ylabel('Correlation')

ax2 = ax1.twinx()
bars2 = ax2.bar(ind+width/2, stds[:n], width, color='g')
ax2.set_ylabel('Standard Deviation')

ax1.legend([bars1, bars2], ['Correlation', 'Standard Deviation'], loc='upper right')
plt.show()

#################

width = .3
ind = np.arange(len(corrs) - n)

ax1 = plt.subplot()
bars1 = ax1.bar(ind-width/2, corrs[n:], width, color='r')
ax1.set_xticks(ind)
ax1.set_xticklabels(names[n:])
ax1.set_xlabel('Traces')
ax1.set_ylabel('Correlation')

ax2 = ax1.twinx()
bars2 = ax2.bar(ind+width/2, stds[n:], width, color='g')
ax2.set_ylabel('Standard Deviation')

ax1.legend([bars1, bars2], ['Correlation', 'Standard Deviation'], loc='upper right')
plt.show()
