import os, sys, json, re
from shutil import copyfile

pathin = sys.argv[1]
print(pathin)

keys_=[]
values_=[]

with open('./samples2017.json', 'r') as fin:
  data=fin.read()
  lines=json.loads(data)
  keys=lines.keys()
  for key, value in lines.items():
    keys_.append(key)
    values_.append(value)

for ik in range(0,len(keys_)):
  tmp_cmd='%s%s_cfg.py'%(pathin,keys_[ik])
  os.system('echo "crab submit -c %s" >> aa.sh'%(tmp_cmd))
  os.system('echo "rm crab_%s/inputs/*.tgz" >> aa.sh'%(values_[ik][0]))
  
