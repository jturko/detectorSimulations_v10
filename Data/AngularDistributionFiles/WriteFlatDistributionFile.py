import sys
import math

# command line parameters
nbOfStates = int(sys.argv[1])
outputFile = sys.argv[2]

# open output file
outData = open(outputFile, "w")

for i in range(0, nbOfStates) :
    outData.write("179 \t")

outData.write("\n")

for th in range(1, 180) :
    for i in range(0, nbOfStates) :
        outData.write(str(th * 1.0) + "\t 1.00 \t\t")# + 
                      #str(math.sin(th*math.pi / 180.)) + "\t\t")
    outData.write("\n")

# close file
outData.close()

