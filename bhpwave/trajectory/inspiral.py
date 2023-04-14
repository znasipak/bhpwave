from bhpwaveformcy import (TrajectoryDataPy,
                           InspiralGeneratorPy)
from bhpwave.constants import *

class InspiralGenerator:
    """
    A class for generating quasi-circular inspirals around a rotating massive black hole.
    Once instantiated, the class can be called to generate inspiral data.

    :param traj_data: A TrajectoryDataPy class which holds interpolants of the relevant trajectory data
    :type traj_data: TrajectoryDataPy or None, optional
    """
    def __init__(self, traj_data = None):
        if traj_data is None:
            self.traj_data = TrajectoryDataPy(dealloc_flag=False)
        else:
            self.traj_data = traj_data
        self.inspiral_generator = InspiralGeneratorPy(self.traj_data)

    def __call__(self, M, mu, a, r0, dt, T, num_threads=0):
        """
        Generates a quasi-circular inspiral of a point-particle around a rotating massive black hole.

        :param M: the mass of the rotating massive black hole (in solar masses)
        :type M: double
        :param mu: the mass of the smaller compact object (in solar masses)
        :type mu: double
        :param r0: the initial Boyer-Lindquist radius of the small body
        :type r0: double
        :param dt: time step in seconds for sampling the trajectory
        :type dt: double
        :param T: duration of the inspiral in years
        :type T: double

        :return: A class object containing trajectory data
        :rtype: InspiralContainerWrapper
        """
        massratio = mu/M
        dtM = dt/(M*Modot_GC1_to_S)
        TM = T*yr_MKS/(M*Modot_GC1_to_S)

        return self.inspiral_generator(massratio, a, r0, dtM, TM, num_threads)