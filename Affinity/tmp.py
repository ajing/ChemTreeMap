infile = "test.json"

fmax = 0
lmax = 0

flist = []
llist = []

for line in open(infile):
    if "position" in line:
        content = line.strip().split(":")[-1].strip()[:-1]
        first, last = map( float, content.replace("\"", "").split("-"))
        flist.append(first)
        llist.append(last)

print max(flist), min(flist), max(llist), min(llist)


# for link distance
linklist = []
for line in open(infile):
    if "dist" in line:
        content = line.strip().split(":")[-1].strip()[:-1]
        dist = float(content)
        linklist.append(dist)

print max(linklist),  min(linklist)

