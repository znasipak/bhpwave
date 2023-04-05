from bhpwaveformcy import (TrajectoryDataPy,
                           InspiralGeneratorPy)

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

    def __call__(self, massratio, a, r0, dt, T, num_threads=0):
        """
        Generates a quasi-circular inspiral of a point-particle around a rotating massive black hole.

        :param massratio: the small massratio :math:`\mu/M` of the binary
        :type massratio: double
        :param r0: the initial Boyer-Lindquist radius of the small body
        :type r0: double
        :param dt: time step in seconds for sampling the trajectory
        :type dt: double
        :param T: duration of the inspiral in years
        :type T: double

        :return: A class object containing trajectory data
        :rtype: InspiralContainerWrapper
        """
        return self.inspiral_generator(massratio, a, r0, dt, T, num_threads)