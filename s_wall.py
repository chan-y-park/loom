import logging

class SWall(object):
    def __init__(self, z_0=None, x1_0=None, x2_0=None, parent=None, 
                 label=None,):
        self.data = [[z_0, x1_0, x2_0]]
        self.parent = parent
        self.label = label

    def get_zs(self, ti=0, tf=None):
        if tf == None:
            tf = len(self.data)
        zs = []
        for t in range(ti, tf):
            zs.append(self.data[t][0])
        return zs 

#    def grow(self, t, step=None, relax=None, integrator='zvode',
#               method='adams', with_jacobian=False):
#        logging.debug('need to implement')

