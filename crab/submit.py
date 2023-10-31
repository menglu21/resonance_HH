import os,sys

CONFIG_PATH=sys.argv[1]
skip_signal=sys.argv[2]
def get_configs():
  for root, dirs, files in os.walk(CONFIG_PATH,topdown = False):
        samples=files
  return samples

configs=get_configs()

for ifile in configs:
  if skip_signal==1 and 'Signal' in ifile:continue
  if '.pyc' in ifile:continue
  cmds='crab submit -c %s/%s'%(CONFIG_PATH,ifile)
  cmds2="grep requestName %s/%s | grep -o \"'[^']*'\""%(CONFIG_PATH,ifile)
  output = os.popen(cmds2).read().replace("\n", "").replace("'", "")
  cmds3='rm crab_%s/inputs/*.tgz'%(output)
  print('preparing submiting process ',ifile)
  if os.path.exists('crab_'+output):
    print('crab already exsits, skip this process')
    continue
    print('submit jobs:',cmds)
    os.system(cmds)
    print('rm tgzs:',cmds3)
    os.system(cmds3)
