import logging

class SWall(object):
    def __init__(self, sw_curve=None, sw_diff=None, theta=0, z_0=None,
                 x1_0=None, x2_0=None):
        self.sw_curve = sw_curve
        self.sw_diff = sw_diff
        self.theta = theta
        self.z = [z_0,]
        self.x1 = [x1_0,]
        self.x2 = [x2_0,]

    def evolve(self, t, step=None, relax=None, integrator='zvode',
               method='adams', with_jacobian=False):
        logging.debug('need to implement')

