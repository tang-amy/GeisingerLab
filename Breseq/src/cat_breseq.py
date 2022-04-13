## Yunfei Dai
## 2022/04/06

from os import listdir
from os.path import isfile, isdir, join
DIR = "/Users/geisingerlab/Yunfei/YFWGS-2022-04-02"
for item in listdir(DIR):
    file=join(DIR,item)
    if isdir(file)==True:
        print(join(file,"output","index.html"))
