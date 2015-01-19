import logging


class SWall(object):
    def __init__(self, z_0=None, x1_0=None, x2_0=None, parent=None,
                 label=None,):
        self.data = [[z_0, x1_0, x2_0]]
        self.parent = parent
        self.label = label

    def get_zs(self, ti=0, tf=None):
        if tf is None:
            tf = len(self.data)
        zs = []
        for t in range(ti, tf):
            zs.append(self.data[t][0])
        return zs

    def grow(
        self,
        ode,
        ramification_point_zs,
        puncture_point_zs,
        z_range_limits=None,
        num_of_steps=None,
        size_of_small_step=None,
        size_of_large_step=None,
        size_of_neighborhood=None,
    ):
        steps = 0
        rpzs = ramification_point_zs
        ppzs = puncture_point_zs
        z_i, x1_i, x2_i = self.data[-1]
        ode.set_initial_value([z_i, x1_i, x2_i])

        while ode.successful() and steps < num_of_steps:
            if (len(rpzs) > 0 and 
                min([abs(z_i - rpz) for rpz in rpzs]) < size_of_neighborhood):
                dt = size_of_small_step
            else:
                dt = size_of_large_step
            z_i, x1_i, x2_i = ode.integrate(ode.t + dt)
            self.data.append([z_i, x1_i, x2_i])
            steps += 1
