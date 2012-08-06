from fipy import *
import tables
import os


class DictTable:
    """
    Designed to save a dictionary of arrays at each timestep in a simulation.
    """
    IDprefix = 'ID'
    def __init__(self, h5filename, mode='a'):
        if mode == 'w' and os.path.exists(h5filename):
            print 'removing %s ' %  h5filename
            os.remove(h5filename)

        self.h5filename = h5filename
        
    def __setitem__(self, index, values):
        h5file = tables.openFile(self.h5filename, mode='a')
        h5file.root._v_attrs.latestIndex = index

        groupName = self.IDprefix + str(index)

        if hasattr(h5file.root, groupName):
            group = h5file.root._f_getChild(groupName)
            group._f_remove(recursive=True)
            
        group = h5file.createGroup(h5file.root, groupName)

        for k in values.keys():
            h5file.createArray(group, k, values[k])

        h5file.close()
        
    def __getitem__(self, index):
        h5file = tables.openFile(self.h5filename, mode='r')

        d = {}
        for array in h5file.listNodes('/' + self.IDprefix + str(index), classname='Array'):
            d[array.name] = array.read()

        h5file.close()
        return d

    def getLatestIndex(self):
        h5file = tables.openFile(self.h5filename, mode='r')
        latestIndex = h5file.root._v_attrs.latestIndex
        h5file.close()
        return latestIndex
    
