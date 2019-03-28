
import ROOT as R
import sys
import math
import array


# ROOT plotting configurations
R.gROOT.SetStyle("Plain")
R.gStyle.SetPalette(1)


# file
fileIn = sys.argv[1]
nbOfStates = int(sys.argv[2])

# open ROOT files
data = R.TFile(fileIn)

# list with histos
graphsAllPoints = []

# get graphs
for i in range(0, nbOfStates) :
    graphsAllPoints.append(data.Get("Graph_" + str(i)))


print graphsAllPoints

graphs = []

for i in range(len(graphsAllPoints)) :
    graphs.append(R.TGraph(graphsAllPoints[i].GetN() / 4))
    for p in range(0, graphsAllPoints[i].GetN(), 4) :
        graphs[i].SetPoint(p/4, graphsAllPoints[i].GetX()[p], graphsAllPoints[i].GetY()[p])


print "NbOfPoints = ", graphs[0].GetN() 

###############################################################################
# statistic boxes, line colors, title offset
###############################################################################
colours = [1,2,R.kAzure+7, 8]
lineStyles = [1, 5, 7, 10, 5, 7]

for i in range(len(graphs)) :
    #graphs[i].SetStats(0)
    graphs[i].SetTitle("")
    graphs[i].SetLineColor(colours[i])
    graphs[i].SetMarkerStyle(20 + i)
    graphs[i].SetMarkerColor(colours[i])
    graphs[i].SetLineWidth(3)
    graphs[i].SetLineStyle(lineStyles[i])
    #graphs[i].Rebin(4)
    graphs[i].GetXaxis().SetTitle("#theta_{CM} [deg]")
    graphs[i].GetYaxis().SetTitle("d#sigma/d#Omega [mb/sr]")
    graphs[i].GetXaxis().SetTitleSize(0.05)
    graphs[i].GetYaxis().SetTitleSize(0.05)
    graphs[i].GetXaxis().SetTitleOffset(0.9)
    graphs[i].GetYaxis().SetTitleOffset(0.9)
    #graphs[i].GetXaxis().SetLimits(0, 180)
    graphs[i].GetYaxis().SetRangeUser(0.02, 20)
    graphs[i].GetXaxis().SetRangeUser(0, 180)
    #graphs[i].SetMinimum(0)

#graphs[0].SetFillColor(R.kGray)


###############################################################################
# legend
###############################################################################
legend = R.TLegend(0.656627 ,0.625 ,0.961847,0.924107) 

legend.SetHeader("^{133}Sn:")
names = ["g.s.", "854 keV", "1363 keV", "2005 keV"]

for i in range(len(graphs)) :
    legend.AddEntry(graphs[i], names[i], "L")
    
legend.SetFillColor(0)
legend.SetLineColor(0)
legend.SetShadowColor(0)
legend.SetTextSize(0.055) # 0.05


###############################################################################
# canvas
###############################################################################
can = R.TCanvas("can", "Miniball", 10, 10, 1000, 700)
can.SetTopMargin(0.05)
can.SetRightMargin(0.03)
can.SetLogy()
can.cd()


graphs[0].Draw("axis C")
graphs[0].Draw("C same")
graphs[1].Draw("L same")
graphs[2].Draw("L same")
graphs[3].Draw("L same")

legend.Draw()
can.Update()


raw_input("Press ENTER to continue.") 
    
