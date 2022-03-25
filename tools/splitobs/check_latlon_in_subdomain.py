#!/usr/bin/env python3

import math
import sys

tolerance=0.01

if len(sys.argv)<6 or len(sys.argv)>7:
    sys.stderr.write("The script 'check_latlon_in_subdomain' takes only 5 or 6 arguments.  They are:\n")
    sys.stderr.write("filein: SQL output containing the information on each observation\n")
    sys.stderr.write("ilonband: number given by splitobs for the longitude band\n")
    sys.stderr.write("jlatband: number given by splitobs for the latitude band\n")
    sys.stderr.write("npex: number of longitude bands\n")
    sys.stderr.write("npey: number of latitude bands\n")
    sys.stderr.write(" verbose: 0 or 1 (1 means verbose and default to 0)\n")

filein=sys.argv[1]
ilonband=int(sys.argv[2])
jlatband=int(sys.argv[3])
npex=int(sys.argv[4])
npey=int(sys.argv[5])
if len(sys.argv)==7:
    verbose = True
else:
    verbose = False

def find_subdomain(lat,lon,npex,npey):
    if verbose:
        print "find_subdomain: lat=%f lon=%f npex=%d npey=%d" % (lat,lon,npex,npey)

    ilonband_r8 = ((lon+360)%360)/(360./npex)+1
    ilonband = int(math.floor(ilonband_r8))
    jlatband_r8 = (lat+90)/(180./npey)+1
    jlatband = int(math.floor(jlatband_r8))

    if verbose:
        message="find_subdomain: ((lon+360)%%360)/(360./npex)=%f (lat+90)/(180./npey)=%f returning ilonband=%d jlatband=%d" % (((lon+360)%360)/(360./npex),(lat+90)/(180./npey),ilonband,jlatband)
        print(message)

    return ilonband, jlatband, ilonband_r8, jlatband_r8

    
def check_latlon_in_subdomain(lat,lon,ilonband,jlatband,npex,npey):
    if verbose:
        message="check_latlon_in_subdomain: lat=%f lon=%f ilonband=%d jlatband=%d npex=%d npey=%d" % (lat,lon,ilonband,jlatband,npex,npey)
        print(message)

    if npex==1 and ilonband!=1:
        print("npex=1 but ilonband=%d" % ilonband)
        return False
    if npey==1 and jlatband!=1:
        print("npey=1 but jlatband=%d" % jlatband)
        return False

    thislonband, thislatband, lonband_r8, latband_r8 = find_subdomain(lat,lon,npex,npey)
    if npey!=1 and thislatband!=jlatband:
        print("jlatband=%d and thislatband=%d ecartlat=%f" % (jlatband,thislatband,jlatband-latband_r8))
        if abs(jlatband-latband_r8)<tolerance:
            return True
        else:
            return False
    if npey!=1 and thislatband!=jlatband:
        print("ilonband=%d and thislonband=%d ecartlon=%f" % (ilonband,thislonband,ilonband-lonband_r8))
        if abs(ilonband-lonband_r8)<tolerance:
            return True
        else:
            return False

    if verbose:
        print("check_latlon_in_subdomain: observation in the correct subdomain")

    return True
       
for line in open(filein).readlines():
    lat = float(line.split('|')[0])
    lon = float(line.split('|')[1])
    
    check = check_latlon_in_subdomain(lat,lon,ilonband,jlatband,npex,npey)
    if check and verbose:
        message="The function check_latlon_in_subdomain is returning True for args lat=%f lon=%f ilonband=%d jlatband=%d npex=%d npey=%d" % (lat,lon,ilonband,jlatband,npex,npey)
        print(message)

    if not check:
        message="The function check_latlon_in_subdomain is returning False for args lat=%f lon=%f ilonband=%d jlatband=%d npex=%d npey=%d" % (lat,lon,ilonband,jlatband,npex,npey)
        print(message)
        sys.exit(1)
