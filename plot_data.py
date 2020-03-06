from matplotlib import pyplot as plt
import csv
import os

diff = []

filelist = os.listdir()
csvlist =[f for f in filelist if '.csv' in f]

plt.ylim(0,10)

for i, filename in enumerate(csvlist):
    with open(filename) as csvfile:
        reader = csv.DictReader(csvfile)
        size = []
        gflops = []
        for row in reader:
            size.append(int(row['Size']))
            gflops.append(float(row['Gflops'].split('e')[0]))
            diff.append(row['Diff'])
        plt.plot(size, gflops,label=filename, color='rgbcmy'[i])

plt.xlabel('Matrix Size')
plt.ylabel('GFlops/sec')
plt.legend(loc='upper right')
plt.show()