import numpy as np
import hicstraw
# --------------------------------------------------------------------------------------------------
#HIC_FILE = "/home/data/Rao/GSE63525_GM12878_insitu_primary.hic"
#HIC_FILE = "/home/data/encode/HiC/ENCSR968KAY/ENCFF216QQM.hic"
HIC_FILE = "/home/data/4DN/HiC/4DNES4AABNEZ/4DNFIC3JD6O2.hic"
SUBJECT = "region_of_interest"
CHR = "1"
START = 178421513
END = 179491193
RES = 5000
NAME = "{0:s}_chr{1:s}_{2:d}-{3:d}_res{4:d}bp".format(SUBJECT, CHR, START, END, RES)
# --------------------------------------------------------------------------------------------------
hic = hicstraw.HiCFile(HIC_FILE)
mzd = hic.getMatrixZoomData(CHR, CHR, "observed", "KR", "BP", RES)
input_matrix = mzd.getRecordsAsMatrix(START, END - RES, START, END - RES)
records_list = mzd.getRecords(START, END - RES, START, END - RES)

np.savetxt("{0:s}.txt".format(NAME), input_matrix, fmt="%e")
print("Contact matrix size is {0:d}x{0:d}".format(input_matrix.shape[0]))

fp = open("{0:s}_list.txt".format(NAME), "w")
for i in range(len(records_list)):
    print("{0:d}\t{1:d}\t{2:f}".format(records_list[i].binX, records_list[i].binY, records_list[i].counts),
          file=fp)
fp.close()
