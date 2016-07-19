#!/uboone/app/users/vgenty/bin/python

#by vic
#so this crunch bang is for vic...

import os,sys,re
import pprint as pp

########READ XML
XMLFILE=sys.argv[1]

contents = ""
with open(XMLFILE,"r") as file_:
    contents = file_.read() #why use etree if we don't have too

project = {}

for item in contents.split("!ENTITY")[1:-1]:
    s = re.search("(.+)\s+\"(.+)\"",item)
    key_ = s.group(1).strip()
    val_ = s.group(2).strip()
    project[key_] = val_

jobs    = float(re.search("<numjobs>([0-9]+)</numjobs>",contents).group(1))
project["jobs"] = jobs

logs    = str(re.search("<logdir>(.*)</logdir>",contents).group(1))
outs    = str(re.search("<outdir>(.*)</outdir>",contents).group(1))

for item in project:
    if re.search(item,logs) is not None:
        logs = logs.replace("&%s;"%item,project[item])

    if re.search(item,outs) is not None:
        outs = outs.replace("&%s;"%item,project[item])

project["logs"] = logs        
project["outs"] = outs
print "\n\t===> Project data for %s<===\n"%XMLFILE
pp.pprint(project)

########Check on the JOB
print "\n\t===> Check on the job <===\n"
LOGDIR=project["logs"]

print "\n\t===> How many files are missing? <===\n"
folders = [os.path.join(LOGDIR,f) for f in os.listdir(LOGDIR) if os.path.isdir(os.path.join(LOGDIR,f))]
folders = [f for f in folders if not f.endswith("_start")]
folders = [f for f in folders if not f.endswith("_stop")]

nmissing = project["jobs"] - len(folders)
nmissing = float(nmissing)
if nmissing < 0.0:
    print "Really? There is more output than number of jobs? Ok..."

badfraction = (nmissing / project["jobs"]) * 100.0

print "Missing jobs are: %f percent for %f out of %f"% (badfraction,nmissing,project["jobs"])
print ""

#Of the ones that are there, is the EXIT code OK?
print"\n\t===> Of the ones that are there, is the EXIT code OK? <===\n"

nbad = 0.0
bads = []

for folder in folders:

    #print "Checking folder ",folder
    larstat = os.path.join(folder,"lar.stat")

    if not os.path.exists(larstat):
        print "Missing larstat in ",folder
        nbad+=1.0
        continue

    data = None

    with open(larstat,'r') as f_:
        data = f_.read().strip()
    
    data = int(data)

    if data != 0:
        print "!!!Bad project larout in ",folder, " with data ",data,"!!!"
        nbad += 1.0
        bads.append(folder.split("/")[-1])
    
badfraction = (nbad / project["jobs"]) * 100.0
print "Bad larout jobs are: %f percent for %f out of %f"% (badfraction,nbad,project["jobs"])
badfraction = ((nbad + nmissing) / project["jobs"])*100.0
print "Now BAD jobs is: %f percent for %f out of %f"% (badfraction,nbad+nmissing,project["jobs"])



#Of the ones that exist, are the larlite files > 1k?
print "\n\t===> Of the ones that exist, are the larlite files > 1k byte? <===\n"
OUTDIR = project["outs"]
folders = [os.path.join(OUTDIR,f) for f in os.listdir(OUTDIR) if os.path.isdir(os.path.join(OUTDIR,f))]
folders = [f for f in folders if not f.endswith("_start")]
folders = [f for f in folders if not f.endswith("_stop")]

print "Already bad ones are\n",bads

lls = ["mcinfo","simch","wire","opdigit","opreco"]

if project["file_type"] == "data": 
    lls = lls[2:]

nsmall = 0.0

for folder in folders:
    
    if folder.split("/")[-1] in bads:
        continue

    files = [os.path.join(folder,f) for f in os.listdir(folder) if f.startswith("larlite")]
    
    for l in files:
        if l.split("/")[-1].split("_")[1] not in lls:
            nsmall+=1.0
            print "\t>larlite file %s doesn't exist "%l
            continue

        size = os.path.getsize(l)
        if size < 1024:
            print "\t>larlite file %s is %s too small!"%(l,size)
            nsmall += 1.0
            bads.append(folder)
            break

print "Now bad ones are\n",bads
badfraction = (nsmall / project["jobs"]) * 100.0
print "Small or missing larlite files are ",badfraction
nbad = float(len(bads))
badfraction = ( (nbad + nmissing) / project["jobs"] )* 100.0
print ""
print "\t====================================================================="
print "\tNow BAD jobs is: %f percent for %f out of %f"% (badfraction,nbad+nmissing,project["jobs"])
print "\t====================================================================="

