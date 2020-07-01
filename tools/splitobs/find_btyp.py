#!/usr/bin/python

btyps=[9312,15456,9322]

## obs, flag, O-P, obs, flag,O-P, relie a la surface
btyps=[9312,15456,9322,98,6242,106,3107]

btyps_obs = [9312,98]
btyps = [15456,9322,6242,106,3107]

## newbtyp_obs=ishft (btyp, -4)
## newbtyp_obs=iand(127,newbtyp_obs)
## 
## newbtyp=ishft (btyp_bloc, -4)
## newbtyp=iand(127,newbtyp)
## 
## si  newbtyp=newbtyp_obs
## ce bloc est a splitter

def is2split(btyp_obs, btyp):
    newbtyp_obs = btyp_obs >> 4 & 127;
    print 'newbtyp_obs = ', newbtyp_obs,
    newbtyp = btyp >> 4 & 127
    print ' newbtyp = ', newbtyp,
    if newbtyp_obs == newbtyp:
        print 'btyp_obs=%d btyp=%d et on splitte' % (btyp_obs,btyp)
    else:
        print 'btyp_obs=%d btyp=%d et on ne splitte pas' % (btyp_obs,btyp)

for btyp_obs in btyps_obs:
    for btyp in btyps:
        is2split(btyp_obs,btyp)
