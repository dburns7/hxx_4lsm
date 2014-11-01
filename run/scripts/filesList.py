#!/usr/bin/env python
# $1 = source dir
# $2 = name

import sys
import os
import shutil

SRCDIR = str(sys.argv[1])
NAME   = str(sys.argv[2])
LIST   = os.listdir(SRCDIR)
FILE   = open(NAME,'w')

cnt = 0
for infile in LIST:
  cnt = cnt + 1
  #print infile
  if cnt%20 == 0:
    FILE.write('inputFiles=file:' + SRCDIR + '/' + infile + ' ' + '\\' + '\n')
  else:
    FILE.write('inputFiles=file:' + SRCDIR + '/' + infile + ' ')

FILE.close()

#os.rename(FILE, TARDIR)
#shutil.move(FILE, TARDIR)
