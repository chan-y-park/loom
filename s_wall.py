import logging

class SWall(object):
    def __init__(self, z_0=None, x1_0=None, x2_0=None, parent=None, 
                 label=None,):
        self.z = [z_0,]
        self.x1 = [x1_0,]
        self.x2 = [x2_0,]
        self.parent = parent
        self.label = label

#    def grow(self, t, step=None, relax=None, integrator='zvode',
#               method='adams', with_jacobian=False):
#        logging.debug('need to implement')

