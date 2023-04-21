import numpy as np



class Rocket:

    def __init__(self, name):
        self.name = name
        self.forces = [0,0,0]
        self.momenta = [0,0,0]
        self.state = [0,0,0,0,0,0]

    def load(self):
        #self.Iz = 1/2*self.mass*(self.radius**2)
        self.Iy = self.mass * (self.radius**2 / 4 + self.height**2 / 12)

    def input(self,thrust,thrust_angle):
        Fx = self.max_thrust * thrust * np.sin(thrust_angle)
        Fy = 0 #Fy = self.max_thrust * thrust * np.sin(np.abs(thrust_angle1)) * np.cos(thrust_angle2)
        Fz = self.max_thrust * thrust * np.cos(thrust_angle)
        
        Mx = 0
        My = self.max_thrust * thrust * np.sin(thrust_angle) * self.cg_thrust_length
        Mz = 0

        self.forces = [Fx,Fy,Fz]
        self.momenta = [Mx,My,Mz]
