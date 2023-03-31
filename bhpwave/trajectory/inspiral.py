from bhpwaveformcy import (TrajectoryData,
                           InspiralGeneratorPy,
                           InspiralContainerWrapper)

class InspiralGenerator:
    def __init__(self, traj_data = None):
        if traj_data is None:
            self.traj_data = TrajectoryData(dealloc_flag=False)
        else:
            self.traj_data = traj_data
        self.inspiral_generator = InspiralGeneratorPy(self.traj_data)

    def __call__(self, *args):
        return self.inspiral_generator(*args)