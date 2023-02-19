#!/bin/env python3

import os
import sys
import re
from argparse import ArgumentParser

parser = ArgumentParser()

# https://martin-thoma.com/how-to-parse-command-line-arguments-in-python/
# Add more options if you like
parser.add_argument("--match", dest="match", default="",
                    help="When reading CRAB tasks, take only tasks whose names contain this string")

opts = parser.parse_args()

# https://stackoverflow.com/questions/287871/how-to-print-colored-text-to-the-terminal
class colors:
   colordict = {
                'RED'        : '\033[91m',
                'GREEN'      : '\033[92m',
                'BLUE'       : '\033[34m',
                'GRAY'       : '\033[90m',
                'WHITE'      : '\033[00m',
                'ORANGE'     : '\033[33m',
                'CYAN'       : '\033[36m',
                'PURPLE'     : '\033[35m',
                'LIGHTRED'   : '\033[91m',
                'PINK'       : '\033[95m',
                'YELLOW'     : '\033[93m',
                'BLINK'      : '\033[5m' ,
                'NORMAL'     : '\033[28m' ,
                "WARNING"    : '\033[33m',
                "CEND"       : '\033[0m',
                }

def getFinalCRABDir(opts, crabdir):
   taskNames  = crabdir
   
   if opts.match != "":
      if opts.match.lower() in crabdir.lower():
         return crabdir
      #else:
      #   raise Exception("Did NOT find proper match, try to remove some string") #FIXME
   else:
      return taskNames   
   
filesfolders = os.listdir(".")
folders = []
for filen in filesfolders:
    if os.path.isdir(os.path.join(os.path.abspath("."), filen)):
        if filen.startswith("crab_"):
           
           finalDir = getFinalCRABDir(opts, filen)
           if finalDir != None:
              print colors.colordict["GREEN"] + "Checking status for crab directory: " + colors.colordict['CEND'], colors.colordict['BLUE'] + finalDir + colors.colordict['CEND']
           
              print(os.popen("crab status -d "+filen).read())
              print 50*"="
           
            
