import tauFun as tf
import sys

#outname = "cat.txt"
try:
    cat = sys.argv[1]
except:
    print("Error: please specify category name.")
    sys.exit()
catn = tf.catToNumber(cat)
print(catn)

#outf = open(outname, "w")
#outf.write("%d"%catn)
#outf.close()
