import numpy as np

from Quat import Quat

rho = 1.293 #kg/m^3


class Rocket:
    def __init__(self, name):
        self.name = name
        self.F = np.zeros((3,1))
        self.M = np.zeros((3,1))
        self.inputs = np.zeros(3)

    def load(self):
        self.Aref = np.pi * self.radius**2
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
        Tx = self.max_thrust * self.thrust * np.sin(self.alpha) * np.cos(self.beta)
        Ty = self.max_thrust * self.thrust * np.sin(self.alpha) * np.sin(self.beta)
        Tz = self.max_thrust * self.thrust * np.cos(self.alpha)
        self.T = [Tx,Ty,Tz]

    def updateForces(self):
        self.F = self.T + self.Df

    def updateMomenta(self):
        self.M = self.cross(np.array([0,0,-self.cg_thrust_length]),self.T)
        print("M",self.M)
        self.M = self.M + self.cross(np.array([0,0,self.cg_cp_length]),self.Df)
        print("Cross",self.cross(np.array([0,0,self.cg_cp_length]),self.Df))

    def update(self,state):
        self.state = state
        self.aerodynamics()
        self.updateForces()
        self.updateMomenta()

    def aerodynamics(self):
        #low velocity
        if np.linalg.norm(self.state[3:6]) < 0.1:
            self.Df = np.array([0,0,0])
            return
        #print("V_inertial",self.state[3:6])
        V = Quat.intertialRotateInv(self.state[6:10],self.state[3:6])
        #print("V",V)
        D = 0.5*+rho*self.Aref*self.Cd*np.linalg.norm(V)**2
        self.Df = -D * V / np.linalg.norm(V)

    def cross(self,x,y):
        #reversed y-component of cross product (reason not known yet)
        return np.array([x[1]*y[2]-x[2]*y[1],-x[2]*y[0]+x[0]*y[2],x[0]*y[1]-x[1]*y[0]])