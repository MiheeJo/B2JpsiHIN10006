#!/bin/python
import sys,string,re,math,numpy
from os import listdir
from operator import itemgetter

foutput = open("fit_parameters","w")
line = ["rap","pT","cent","NSig","NSigErr","NBkg","NBkgErr","coefExp","coefExpErr","coeffGaus","coeffGausErr","meanSig1","meanSig1Err","sigmaSig1","sigmaSig1Err","sigmaSig2","sigmaSig2Err","alpha","alphaErr","enne","enneErr","enneW","enneWErr","fracRes","fracResErr","fracRes2","fracRes2Err","fracRes3","fracRes3Err","meanResSigW","meanResSigWErr","sigmaResSigW","sigmaResSigWErr","sigmaResSigN","sigmaResSigNErr","sigmaResSigM","sigmaResSigMErr","sigmaResSigO","sigmaResSigOErr","fLiving","fLivingErr","fpm","fpmErr","fbkgCtTot","fbkgCtTotErr","lambdam","lambdamErr","lambdap","lambdapErr","lambdasym","lambdasymErr","NLL","Prompt","PromptErr","NonPrompt","NonPromptErr","Bfrac","BfracErr","Resolution","ResolutionErr","nFullBinsResid","RSS"]
wline = ""
for idx in range(0,len(line)):
  wline = wline + str(line[idx]) + "\t"
wline = wline + "\n"
foutput.write(wline)

ftable = open("fit_table","w")
line2 = ["rap","pT","cent","NSig","NSigErr","NBkg","NBkgErr","NLL","PromptJ/psi","PromptJ/psiErr","Non-promptJ/psi","Non-promptJ/psiErr","Bfrac","BfracErr"]
wline = ""
for idx in range(0,len(line2)):
  wline = wline + str(line2[idx]) + "\t"
wline = wline + "\n"
ftable.write(wline)
'''
analysis_bins = ("0.0-2.4", "6.5-30.0", "0-100"), ("0.0-2.4", "6.5-10.0", "0-100"), ("0.0-2.4", "10.0-30.0", "0-100"),\
("0.0-1.2", "6.5-30.0", "0-100"), ("1.2-1.6", "5.5-30.0", "0-100"), ("1.2-1.6", "6.5-30.0", "0-100"), ("1.6-2.4", "3.0-30.0", "0-100"),("1.6-2.4", "6.5-30.0", "0-100"),\
("0.0-2.4", "6.5-30.0", "0-10"), ("0.0-2.4", "6.5-30.0", "10-20"),("0.0-2.4", "6.5-30.0", "20-30"),("0.0-2.4","6.5-30.0","30-40"),("0.0-2.4","6.5-30.0","40-50"),("0.0-2.4", "6.5-30.0", "50-100"),("0.0-2.4", "6.5-30.0", "0-20"),("0.0-2.4", "6.5-30.0", "20-100")


mean = numpy.zeros((len(analysis_bins),len(line2)));
rms =  numpy.zeros((len(analysis_bins),len(line2)));
rms_base =  numpy.zeros((len(analysis_bins),len(line2)));
'''

for f in range(1,len(sys.argv)):
  folder = sys.argv[f]
  print folder
  foutput.write(folder + "\n")
  ftable.write(folder + "\n")

  # Reading 1 fitting method, various bins
  filelist = listdir(folder)
  datapar = []
  datatab = []
  for files in filelist:  # 1 bin per a line
    data = []
    table = []
    if files == "fit_parameters":
      continue
    if files == "fit_table":
      continue
    fname = files.split("_")
    rap = re.search('\d\.\d-\d\.\d',fname[0]).group()
    cent = re.search('\d+-\d+',fname[1]).group()
    pt = re.search('\d+.\d-\d+.\d',fname[2]).group()
    data.append(rap)
    data.append(pt)
    data.append(cent)
    table.append(rap)
    table.append(pt)
    table.append(cent)

    try:
      finput = open("./"+folder+"/"+files+"/2D_GG.txt",'r')
    except:
      print "./"+folder+"/"+files+"/2D_GG.txt isn't able to open"
      continue;

    for line in finput: # Read one 2D_GG.txt file
      tmp = line.split(" ")
      for i in tmp:
        if i == "NLL":
          try:
            table.append (float(tmp[-1]))
          except:
            continue
        if i == "NSig" or i == "NBkg" or i == "PROMPT" or i == "NON-PROMPT" or i == "Bfraction":
          try:
            table.append (float(tmp[-2]))
            table.append (float(tmp[-1]))
          except:
            continue
        try:
          data.append( float(i) )
        except:
          continue

    # put 1 line for each bin (1 2D_GG.txt file)
    datapar.append(data)
    datatab.append(table)

  # sort ALL directory's results and put them into a file
  dataparfin = sorted(datapar, key=itemgetter(0,1,2))
  for i in dataparfin:
    line = ""
    for j in i:
      line = line+str(j)+"\t"
    foutput.write(line)
    foutput.write("\n")

  # sort results of only useful bins and put them into file
  datatabfin = sorted(datatab, key=itemgetter(0,1,2))
  for ridx in range(0,len(datatabfin)):
    row = datatabfin[ridx]
    rap = row[0]
    pt = row[1]
    cent = row[2]
    line = ""
    for i in row:   # from row[0] to the last column of row
      line = line+str(i)+"\t"
    '''
    for cidx in range(3,len(row)):   # from row[3] to the last column of row
      mean[ridx,cidx] = mean[ridx,cidx] + row[cidx];
    if f is 1:  #1st file is the default value
      for cidx in range(3,len(row)):   # from row[3] to the last column of row
        rms_base[ridx,cidx] = row[cidx];
    else:
      for cidx in range(3,len(row)):   # from row[3] to the last column of row
        rms[ridx,cidx] = rms[ridx,cidx] + math.pow(math.fabs(row[cidx]-rms_base[ridx,cidx]),2)
    '''            
    ftable.write(line)
    ftable.write("\n")
'''
# All files are processed
ftable.write("\n")
ftable.write("MEAN\n")
for ridx in range(0,len(datatabfin)):   # write mean values
  row = datatabfin[ridx]
  rap = row[0]
  pt = row[1]
  cent = row[2]
  for bins in analysis_bins:
    if (rap, pt, cent) == bins:
      line = str(rap) + "\t" + str(pt) + "\t" + str(cent) + "\t"
      for cidx in range(3,len(row)):   # from row[0] to the last column of row
        line = line+str(mean[ridx,cidx]/(len(sys.argv)-1))+"\t"
      ftable.write(line)
      ftable.write("\n")

ftable.write("\nRMS\n")
for ridx in range(0,len(datatabfin)):   # write rms values
  row = datatabfin[ridx]
  rap = row[0]
  pt = row[1]
  cent = row[2]
  for bins in analysis_bins:
    if (rap, pt, cent) == bins:
      line = str(rap) + "\t" + str(pt) + "\t" + str(cent) + "\t"
      for cidx in range(3,len(row)):   # from row[0] to the last column of row
        line = line+str(math.sqrt(rms[ridx,cidx]/(len(sys.argv)-2)))+"\t"
      ftable.write(line)
      ftable.write("\n")

ftable.write("\nRMS/Average\n")
for ridx in range(0,len(datatabfin)):   # write rms values
  row = datatabfin[ridx]
  rap = row[0]
  pt = row[1]
  cent = row[2]
  for bins in analysis_bins:
    if (rap, pt, cent) == bins:
      line = str(rap) + "\t" + str(pt) + "\t" + str(cent) + "\t"
      for cidx in range(3,len(row)):   # from row[0] to the last column of row
        line = line+str( (math.sqrt(rms[ridx,cidx]/(len(sys.argv)-2))) / (mean[ridx,cidx]/(len(sys.argv)-1)) )+"\t"
      ftable.write(line)
      ftable.write("\n")
'''
