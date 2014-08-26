import os,commands,subprocess,sys
from subprocess import Popen, PIPE,STDOUT

TMP_WORK_DIR="tmp_mcshower"
if len(sys.argv) < 5:
    print 'Usage: %s $DATA_TOP_DIR $FILENAME $NUM_COMBINE $FCL'
    sys.exit(1)
INPUT_DIR=sys.argv[1]
FILENAME=sys.argv[2]
BATCH=int(sys.argv[3])
FCL_FILE=sys.argv[4]
if not os.path.isfile(FCL_FILE):
    print 'fcl file not found:',FCL_FILE
if not os.path.isdir(INPUT_DIR):
    print 'dir not exist:',INPUT_DIR

dirs=os.listdir(INPUT_DIR)
job_keys=[]
job_files={}
job_proc=-1
for d in dirs:
    if d.find("_")>0:
        proc = d[0:d.find("_")]
        id   = d[d.find("_")+1:len(d)]
        print proc,id
        if not proc.isdigit() or not id.isdigit():
            continue

#        if job_proc==-1: job_proc = proc
#        elif not proc == job_proc:
#            print "ERROR! un-matched proc: %d=>%d" % (proc,job_proc)
#            continue
        fname = "%s/%s/%s" % (INPUT_DIR,d,FILENAME)
        if os.path.isfile(fname):
            job_keys.append(int(id))
            job_files[int(id)]=fname
        else:
            print "ERROR! file not found for a key:", id

print "identified %d files..." % len(job_files)

if os.path.isdir(TMP_WORK_DIR):
    print "ERROR: process on-going..."
    sys.exit(1)

os.system("mkdir %s" % TMP_WORK_DIR)
job_keys.sort()
tmp_files={}
for key in job_keys:
    tmp_fname = "%s/tmp_%04d.root" % (TMP_WORK_DIR,key)
    os.system("ln -s %s %s" % (job_files[key],tmp_fname))
    tmp_files[key]=tmp_fname;

cmd=''
batch_ctr=0
combined_files=[]
for key in job_keys:
    
    cmd+='%s ' % tmp_files[key]

    if key and key%BATCH==0:

        combined_fname = "combined_mcshower_%02d.root" % batch_ctr
        combined_files.append(combined_fname)
        batch_ctr+=1
        if os.path.isfile(combined_fname):
            print "%s already exist! Skipping..." % combined_fname
        else:
            cmd = 'lar -c %s -s %s' % (FCL_FILE,cmd)
            print
            print cmd
            os.system(cmd)
            cmd = "mv larlight_mcshower.root %s" % combined_fname
            print cmd
            os.system(cmd)
        cmd=''
if cmd:

    combined_fname = "combined_mcshower_%02d.root" % batch_ctr
    combined_files.append(combined_fname)
    batch_ctr+=1
    if os.path.isfile(combined_fname):
        print "%s already exist! Skipping..." % combined_fname
    else:        
        cmd = 'lar -c %s -s %s' % (FCL_FILE,cmd)
        print
        print cmd
        os.system(cmd)
        cmd = "mv mcshower.root %s" % combined_fname
        print cmd
        os.system(cmd)
    cmd=''
    
os.system('rm -rf %s' % TMP_WORK_DIR)


        






