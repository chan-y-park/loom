import logging
import pdb


class Joint:
    def __init__(self, z, x1, x2, parents=None, label=None,):
        self.z = z
        self.x1 = x1
        self.x2 = x2
        self.parents = parents 
        self.label = label

    def __eq__(self, other):
        return self.label == other.label


class SWall(object):
    def __init__(self, z_0=None, x1_0=None, x2_0=None, parents=None,
                 label=None,):
        self.data = [[z_0, x1_0, x2_0]]
        self.parents = parents
        self.label = label

    def get_zs(self, ti=0, tf=None):
        """
        return a list of (z.real, z.imag)
        """
        if tf is None:
            tf = len(self.data)
        zs = []
        for t in range(ti, tf):
            z = self.data[t][0]
            zs.append(z)
        return zs

    def get_zxzys(self, ti=0, tf=None):
        """
        return a list of (z.real, z.imag)
        """
        if tf is None:
            tf = len(self.data)
        zxzys = []
        for t in range(ti, tf):
            z = self.data[t][0]
            zxzys.append((z.real, z.imag))
        return zxzys

    def out_of_range(self, z_range_limits):
        z_f = self.data[-1][0]
        z_real_min, z_real_max, z_imag_min, z_imag_max = z_range_limits
        if (z_f.real < z_real_min or
            z_f.real > z_real_max or
            z_f.imag < z_imag_min or
            z_f.imag > z_imag_max):
            return True
        else:
            return False

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
        previous_length = len(self.data)
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

        return previous_length
