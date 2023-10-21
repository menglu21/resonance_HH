import os, sys, json
from shutil import copyfile

year = sys.argv[1]

if year=='2016apv':
  if os.path.isdir('config_crab_2016apv') is False:
      os.system('mkdir config_crab_2016apv')
  workdir='config_crab_2016apv/'
  datajson='Cert_271036-284044_13TeV_Legacy2016_Collisions16_preVPF_JSON.txt'
  samplejson='samples2016apv.json'
  scriptpath='2016apv_script'

if year=='2016':
  if os.path.isdir('config_crab_2016') is False:
      os.system('mkdir config_crab_2016')
  workdir='config_crab_2016/'
  datajson='Cert_271036-284044_13TeV_Legacy2016_Collisions16_postVPF_JSON.txt'
  samplejson='samples2016.json'
  scriptpath='2016_script'

if year=='2017':
  if os.path.isdir('config_crab_2017') is False:
      os.system('mkdir config_crab_2017')
  workdir='config_crab_2017/'
  datajson='Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
  samplejson='samples2017.json'
  scriptpath='2017_script'

if year=='2018':
  if os.path.isdir('config_crab_2018') is False:
      os.system('mkdir config_crab_2018')
  workdir='config_crab_2018/'
  datajson='Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
  samplejson='samples2018.json'
  scriptpath='2018_script'

with open(samplejson, 'r') as fin:
  data=fin.read()
  lines=json.loads(data)
  for key, value in lines.items() :
    if len(value)==3:
      copyfile('data_cfg.py',workdir+key+'_cfg.py')
      value[1]=value[1].replace('/', 'sss')
      os.system(r'sed -i "6s/dummy/%s/g" %s' %(value[0],workdir+key+'_cfg.py'))
      if '_A' in value[0] or '_B' in value[0] or '_C' in value[0] or '_D' in value[0] or '_E' in value[0] or '_F' in value[0] or '_G' in value[0] or '_H' in value[0]:
        os.system(r'sed -i "13s/dummy/%s/g" %s' %(scriptpath+'sss'+value[2],workdir+key+'_cfg.py'))
      os.system(r'sed -i "13s/sss/\//g" %s' %(workdir+key+'_cfg.py'))
      os.system(r'sed -i "15s/dummy/%s/g" %s' %(datajson,workdir+key+'_cfg.py'))
      os.system(r'sed -i "19s/dummy/%s/g" %s' %(value[1],workdir+key+'_cfg.py'))
      os.system(r'sed -i "19s/sss/\//g" %s' %(workdir+key+'_cfg.py'))
      os.system(r'sed -i "23s/dummy/%s/g" %s' %(datajson,workdir+key+'_cfg.py'))
      os.system(r'sed -i "26s/dummy/%s/g" %s' %(value[0],workdir+key+'_cfg.py'))
    
    # for MC
    if len(value)==2:
      print(value)
      copyfile('mc_cfg.py',workdir+key+'_cfg.py')
      value[1]=value[1].replace('/', 'sss')
      os.system(r'sed -i "6s/dummy/%s/g" %s' %(value[0],workdir+key+'_cfg.py'))
      os.system(r'sed -i "13s/dummy/%s/g" %s' %(scriptpath+'sss'+'crab_script.sh',workdir+key+'_cfg.py'))
      os.system(r'sed -i "13s/sss/\//g" %s' %(workdir+key+'_cfg.py'))
      os.system(r'sed -i "19s/dummy/%s/g" %s' %(value[1],workdir+key+'_cfg.py'))
      os.system(r'sed -i "19s/sss/\//g" %s' %(workdir+key+'_cfg.py'))
      os.system(r'sed -i "26s/dummy/%s/g" %s' %(value[0],workdir+key+'_cfg.py'))
