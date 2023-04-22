import numpy as np



class Rocket:
    def __init__(self, name):
        self.name = name
        self.F = np.zeros((3,1))
        self.M = np.zeros((3,1))
        self.inputs = np.zeros(3)

    def load(self):
        #self.Iz = 1/2*self.mass*(self.radius**2)
        self.Iy = self.mass * (self.radius**2 / 4 + self.height**2 / 12)

        self.inertia = np.array([[self.Iy, 0, 0],
                                [0, self.Iy, 0],
                                [0, 0, self.Iy]])
    def input(self,thrust,alpha,beta):
        self.thrust = thrust
        self.alpha = alpha
        self.beta = beta
        self.inputs = [self.thrust,alpha,beta]

    def updateForces(self):
        Fx = self.max_thrust * self.thrust * np.sin(self.alpha) * np.cos(self.beta)
        Fy = self.max_thrust * self.thrust * np.sin(self.alpha) * np.sin(self.beta)
        Fz = self.max_thrust * self.thrust * np.cos(self.alpha)
        self.F = np.array([Fx,Fy,Fz])

    def updateMomenta(self):
        Mx = self.max_thrust * self.thrust * np.sin(self.alpha) * self.cg_thrust_length * np.sin(self.beta)
        My = self.max_thrust * self.thrust * np.sin(self.alpha) * self.cg_thrust_length * np.cos(self.beta)
        Mz = 0
        self.M = np.array([Mx,My,Mz])

    def update(self):
        self.updateForces()
        self.updateMomenta()
