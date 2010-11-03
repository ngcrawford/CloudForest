#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Nick Crawford on 2010-10-28.
Copyright (c) 2010

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses

The author may be contacted at ngcrawford@gmail.com
"""

import os, sys, glob, shlex, subprocess
def mapper(args):
    """sends the individual """
    for line in sys.stdin:
        path = line.strip()
        # pathname = os.path.dirname(sys.argv[0])
        # pathname = os.path.abspath(pathname)
        # print pathname
        # this will probably need to be changed for aws
        command = './FastTree -nt -quiet %s' % (path)
        print command
        command = shlex.split(command)

            
def main():
    if sys.argv[1] == "mapper":
      mapper(sys.argv[2:])


if __name__ == '__main__':
    main()

