from scipy.integrate import solve_ivp
import numpy as np

from Quat import Quat


class Dynamics:
        def __init__(self, rocket):
                self.rocket = rocket
                self.state = np.zeros(13)
                self.state[9] = 1 #initial quaternion for vertical position
                self.states = np.array([self.state])

        #calculates the derivatives of the state vector
        def dState(self, t, state):
                x,y,z,vx,vy,vz,q0,q1,q2,q3,omega_x,omega_y,omega_z = state

                self.rocket.update(state)

                X = np.array([x,y,z])
                V = np.array([vx,vy,vz])
                Q = np.array([q0,q1,q2,q3])
                #normalize quaternion
                Q = Q / np.linalg.norm(Q)
                W = np.array([omega_x,omega_y,omega_z])

                dX = V
                dV = Quat.intertialRotate(Q,self.rocket.F) / self.rocket.mass + self.externalForces()/self.rocket.mass
                dW = np.linalg.inv(self.rocket.inertia) @ (self.rocket.M - np.cross(W, self.rocket.inertia @ W))
                print("dW",dW)
                dQ = 0.5 * np.array([[0,-W[0], -W[1], -W[2]],
                                        [W[0], 0, W[2], -W[1]],
                                        [W[1], -W[2], 0, W[0]],
                                        [W[2], W[1], -W[0], 0]]) @ Q
                
                newState = np.concatenate((dX,dV,dQ,dW))
                self.ground(state,newState)

                return newState
                
        def update(self, dt):
                newState = solve_ivp(self.dState, [0, dt], self.state)
                self.state = newState.y[:,-1]
                self.states = np.concatenate((self.states, [self.state]))
        
        def getQuat(self):
            return self.state[6:10]
    

        def getEuler(self):
            return Quat.eulerAngles2(self.getQuat())

        #Constraint functions
        def ground(self,state,newState):
                if state[2] < 0:
                        newState[0] = 0
                        newState[1] = 0
                        newState[2] = 0
                        newState[6] = state[6]
                        newState[7] = state[7]
                        newState[8] = state[8]
                        newState[9] = state[9]

        def externalForces(self):
                return np.array([0,0,-9.81])*self.rocket.mass