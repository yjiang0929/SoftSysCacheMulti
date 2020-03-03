import matplotlib.pyplot as plt
import csv

size = []
gflops = []
diff = []

with open('output_MMult0.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        size.append(int(row['Size']))
        gflops.append(float(row['Gflops'].split('e')[0]))
        diff.append(row['Diff'])

plt.ylim(0,10)
plt.xlabel('Matrix Size')
plt.ylabel('GFlops/sec')
plt.plot(size, gflops)
plt.show()