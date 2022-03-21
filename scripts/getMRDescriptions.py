#!/usr/bin/env python3
# =================================================================
#
# Terms and Conditions of Use
#
# Unless otherwise noted, computer program source code of this
# distribution is covered under Crown Copyright, Government of
# Canada, and is distributed under the MIT License.
#
# The Canada wordmark and related graphics associated with this
# distribution are protected under trademark law and copyright law.
# No permission is granted to use them outside the parameters of
# the Government of Canada's corporate identity program. For
# more information, see
# http://www.tbs-sct.gc.ca/fip-pcim/index-eng.asp
#
# Copyright title to all 3rd party software distributed with this
# software is held by the respective copyright holders as noted in
# those files. Users are asked to read the 3rd Party Licenses
# referenced with those assets.
#
# Copyright (c) 2015 Government of Canada
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# =================================================================

## Ce script est inspire de ce script
##     http://gitlab.ssc.etg.gc.ca/woudc/woudc/blob/master/backup/backup_issue_tracker.py

## exemples d'appel:
###    ./gitlab-list-branches.py atmospheric-data-assimilation/gdps
###    ./gitlab-list-branches.py atmospheric-data-assimilation/derivate
###    ./gitlab-list-branches.py atmospheric-observations-quality-control/bgck.thinning1_radiosondes

import json
from urllib.request import urlopen
import sys
import os
import re

if len(sys.argv) < 1 or len(sys.argv) > 2:
    sys.stderr.write('Usage: {} [<private_token]\n'.format(sys.argv[0]))
    sys.exit(1)


PROJECT_ID_RAW="git@gitlab.science.gc.ca:atmospheric-data-assimilation/midas.git"
if len(sys.argv)==2:
    PRIVATE_TOKEN = sys.argv[1]
else:
    tokenfile='%s/.gitlab-private-token' % os.environ['HOME']
    if os.path.isfile(tokenfile):
        PRIVATE_TOKEN = open(tokenfile).read().strip()
    else:
        sys.stderr.write("%s: No token as been specified as the third argument and the file '%s' does not exist\n" % ('gitlab-list-branches.py',tokenfile))
        sys.stderr.write("%s: check this wiki page to know how to create the file '%s':\n" % ('gitlab-list-branches.py',tokenfile))
        sys.stderr.write("%s:         https://wiki.cmc.ec.gc.ca/wiki/Git/Doc#Interaction_avec_GitLab_en_utilisant_des_scripts\n\n" % ('gitlab-list-branches.py'))
        sys.exit(200)

if PROJECT_ID_RAW.startswith('http'):
    GITLAB_API = re.sub(r'(https?://[^/]*?)/.*',r'\1/api/v3',PROJECT_ID_RAW)
    PROJECT_PATH = re.sub(r'(https?://[^/]*?)/','',PROJECT_ID_RAW)
elif PROJECT_ID_RAW.startswith('git@'):
    GITLAB_API = re.sub(r'git@(.*?):.*',r'https://\1/api/v3',PROJECT_ID_RAW)
    PROJECT_PATH = re.sub(r'git@(.*?:)','',PROJECT_ID_RAW)
else:
    sys.stderr.write("The gitlab project URL must start with 'http://', 'https://' or 'git@'\n")
    sys.exit(1)

#print('GITLAB_API=%s' % GITLAB_API)
#print('PROJECT_PATH=%s' % PROJECT_PATH)

## This is the list of all issue and merge request that should have been mentionned in the CHANGELOG
##    See https://gitlab.science.gc.ca/atmospheric-data-assimilation/midas/issues/204
issues_mr = [[127, 131],
             [145, 129],
             [135, 135],
             [149, 136],
             [164, 143],
             [151, 145],
             [168, 146],
             [170, 153],
             [ 80, 147],
             [114, 163],
             [167, 166],
             [139, 172],
             [188, 173],
             [184, 177],
             [163, 175],
             [113, 174],
             [173, 176],
             [186, 179],
             [183, 187],
             [193, 190]]

### Il faut transformer certains caracteres avec la convention sous 'http://www.ascii.cl/htmlcodes.htm'
PROJECT_ID = PROJECT_PATH.replace('.git','').replace('/','%2F').replace('.','%2E')

all_mr_json=[]
number=1
while True:
    url='%s/projects/%s/merge_requests?private_token=%s&page=%d' % (GITLAB_API,PROJECT_ID,PRIVATE_TOKEN,number)
    print(url)
    content = urlopen(url).read()
    jsoncontent = json.loads(content)
    if len(jsoncontent)>0:
        all_mr_json.extend(jsoncontent)
    else:
        break
    number+=1

print(len(all_mr_json))

output = open("descriptions","w")
for mr_json in all_mr_json:
    mr_id=int(mr_json['iid'])
    for issue_mr in issues_mr:
        mr=issue_mr[1]
        if mr_id == mr:
            output.write("###########################################\n")
            output.write("(#%d and !%d)\n\n" % (issue_mr[0], issue_mr[1]))
            output.write(mr_json['title'])
            output.write("\n\n")
            output.write(mr_json['description'])
            output.write("\n###########################################\n")
            output.write("\n\n")
    #print(issue_json['title'], issue_json['iid'], issue_json['id'])
output.close()
