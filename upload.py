#!/usr/bin/python3
from secret import *
import sys
import requests

if(len(sys.argv) != 2):
    print("error")
    quit()

files = {'file': open(sys.argv[1], 'rb')}
param = {'token':token, 'channels':channel}
res = requests.post(url="https://slack.com/api/files.upload",params=param, files=files)
