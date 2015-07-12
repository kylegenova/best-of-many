import shutil

shutil.copytree('./eigen', '../dependencies/Eigen')
shutil.copytree('./blossomV', '../dependencies/blossomV')
shutil.copytree('./delaunay', '../dependencies/blossomV/triangle')
shutil.copy('./concorde.exe/concorde.exe', '../dependencies/concorde.exe')
shutil.copy('./blossom-V.vcxproj', '../dependencies/blossomV/blossom-V.vcxproj')

def changeLine(lines, containing, to):
    return map(lambda line: line if containing not in line else to, lines)

def changeLineInFile(file, containing, to):
    f = open(file, 'r')
    lines = f.read().split('\n')
    lines = changeLine(lines,containing,to)
    f.close()
    f = open(file, 'w')
    f.write('\n'.join(lines))
    f.close()

# Stop blossomV from telling triangle to print.
changeLineInFile('../dependencies/blossomV/GEOM/GPMinit.cpp',
    'pczAevn', '\ttriangulate("pczAevnQ", &in, &out, &vorout);')

# Stop blossomV from throwing an assertion error when
# we use it to interact with triangle for other purposes than
# an MPM
changeLineInFile('../dependencies/blossomV/GEOM/GPMinterface.cpp',
    '# of points is odd', '')

# allow us access to hidden parts of blossomV... probably unwise but
# much easier than the alternative; coding up the same thing ourselves
changeLineInFile('../dependencies/blossomV/GEOM/GeomPerfectMatching.h',
    'private:', '')

# Tell blossomV that we have installed Triangle
changeLineInFile('../dependencies/blossomV/GEOM/GeomPerfectMatching.h',
    '//#define DELAUNAY_TRIANGLE', '#define DELAUNAY_TRIANGLE')

# Get rid of the warnings in Triangle
changeLineInFile('../dependencies/blossomV/triangle/triangle.c',
    'exponent = 2.0 * exponent + (length > SQUAREROOTTWO);',
    'exponent = (int)(2.0 * exponent + (length > SQUAREROOTTWO));')
changeLineInFile('../dependencies/blossomV/GEOM/GPMinit.cpp',
    '#include <stdio.h>', '#define _CRT_SECURE_NO_WARNINGS\n#include <stdio.h>')

# Tell Triangle we aren't on Unix
changeLineInFile('../dependencies/blossomV/triangle/triangle.c',
    '/* #define NO_TIMER */', '#define NO_TIMER')

# Some unused blossomV features generate an fopen warning; just remove it
changeLineInFile('../dependencies/blossomV/PMinterface.cpp', '#include <stdio.h>',
    '#define _CRT_SECURE_NO_WARNINGS\n#include <stdio.h>')
