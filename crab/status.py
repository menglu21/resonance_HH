import os,sys

CONFIG_PATH=sys.argv[1]
def get_configs():
  for root, dirs, files in os.walk(CONFIG_PATH,topdown = False):
        samples=files
  return dirs

configs=get_configs()

for idir in configs:
  if not idir.startswith('crab_'):continue
  cmds='crab status -d %s'%(idir)
  print('status jobs:',cmds)
  os.system(cmds)
