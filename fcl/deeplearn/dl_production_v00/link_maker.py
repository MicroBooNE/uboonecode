import os, sys,time

import os,sys,re,commands
import pprint as pp

########READ XML                                                                                                                                                                
XMLFILE=sys.argv[1]
print

if not XMLFILE.endswith('.xml'):
    print '\033[93mERROR\033[00m: input does not have an .xml extension!\n'
    sys.exit(1)
UNIQUE_NAME = XMLFILE.replace('.xml','')
UNIQUE_NAME = UNIQUE_NAME[UNIQUE_NAME.rfind('/')+1:len(UNIQUE_NAME)] + '/'
OUTPUT_DIR  = '/uboone/data/users/kterao/dl_production_symlink_v00/%s' % UNIQUE_NAME
print 'SymLink directory:\033[95m',UNIQUE_NAME,'\033[00m'
if os.path.isdir(OUTPUT_DIR) and len(os.listdir(OUTPUT_DIR)):
    print '\033[93mERROR\033[00m: SymLink directory is not empty... aborting!\n'
    sys.exit(1)
if not os.path.isdir(OUTPUT_DIR):
    os.system('mkdir %s; chmod -R 775 %s' % (OUTPUT_DIR,OUTPUT_DIR))

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

INPUT_DIR=project['outs']

jobdirs = [x for x in os.listdir(INPUT_DIR) if ( len(x.split('_'))==2 and 
                                                 x.split('_')[0].isdigit() and 
                                                 x.split('_')[1].isdigit() )
           ]

lite_files_map={}
ctr = 0
for d in jobdirs:

    jobid = int(d.split('_')[-1])
    jobdir = '%s/%s' % (INPUT_DIR,d)
    lite_files = [x for x in os.listdir(jobdir) if ( ( (x.startswith('larlite') and x.endswith('.root')) or
                                                       (x.startswith('ana_hist'))
                                                       ) and
                                                     os.path.getsize('%s/%s' % (jobdir,x)) > 1000000 )
                  ]
    for f in lite_files:

        flavor = f.split('_')[1]
        #if flavor == 'mc': flavor = 'mcinfo'
        if not flavor in lite_files_map: 
            lite_files_map[flavor] = {}
        if jobid in lite_files_map[flavor]:
            print 'ERROR: duplicate jobid found (%d)' % jobid
            raise Exception
        lite_files_map[flavor][jobid]='%s/%s' % (jobdir,f)

    if len(lite_files):
        ctr += 1
        sys.stdout.write('Job count: %-3d\r' % ctr)
        sys.stdout.flush()
print

print 'Listing found lite file flavors...'
success_list={}
for f in lite_files_map:
    jobids = lite_files_map[f].keys()
    sys.stdout.write('\033[95m%10s\033[00m ... \033[93m%-3d\033[00m files\n' % (f,len(jobids)))
    sys.stdout.flush()

    for j in jobids:
        if not j in success_list:
            success_list[j]=[]
        success_list[j].append(f)

success_jobs=[]
for j,lite_flavors in success_list.iteritems():

    if len(lite_flavors) != len(lite_files_map): continue
    success_jobs.append(j)

print
print 'Found \033[93m%d jobs\033[00m successful for ALL flavors...' % len(success_jobs)    
print
print 'Making sym-link under %s...' % OUTPUT_DIR
print

for flavor,fmap in lite_files_map.iteritems():

    sys.stdout.write('\033[95mProcessing: %-10s\033[00m\r' % flavor)
    sys.stdout.flush()

    ctr=len(success_jobs)
    for jobid in success_jobs:

        target_fname = fmap[jobid]
        if flavor == 'hist':
            target_link = '%s/anatree_%04d.root' % (OUTPUT_DIR,jobid)
        else:
            target_link = '%s/larlite_%s_%04d.root' % (OUTPUT_DIR,flavor,jobid)
        
        cmd = 'ln -s %s %s' % (target_fname,target_link)
        #print cmd
        #break
        os.system(cmd)
        ctr -=1
        sys.stdout.write('\033[95mProcessing: %-10s\033[00m ... Remaining: %-3d/%-3d\r' % (flavor,ctr,len(success_jobs)))
        sys.stdout.flush()
    print 'Processing: %-10s done (%d files)                       ' % (flavor,len(success_jobs))
print
print '\033[95mDone\033[00m'
print
