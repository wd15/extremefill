#!/usr/bin/env python

import sys
import getopt
import time

from os import path

import sphinx

def main(argv=sys.argv):
    filebase = path.splitext(argv[-1])[0]
    master = filebase + ".rst"
    timestamp = time.asctime()
    
    f = open(master, "w")
    f.write("""\
.. %(filebase)s documentation master file, created by
   sphinx-oneoff on %(timestamp)s.

.. autosummary::
   :toctree: generated

   %(filebase)s
""" % locals())
    f.close()
    
    argv = [argv[0]] + ["-D", "source_suffix=.rst", 
                        "-D", "master_doc=" + path.splitext(master)[0],
                        "-D", "latex_documents=[('%(filebase)s', '%(filebase)s.tex', u'%(filebase)s Documentation',u'J. E. Guyer', 'howto'),]" % locals()] + argv[1:-1]
                       
    print >>sys.stderr, argv
                       
    return sphinx.main(argv)
    
if __name__ == '__main__':
    sys.exit(main(sys.argv))
