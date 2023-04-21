from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation as R
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

img = mpimg.imread('falcon9.png')
g = 9.81

class Dynamics:
        def __init__(self, rocket):
                self.rocket = rocket
                self.states = np.array([self.rocket.state])

        #converts quaternion to euler angles
        def eulerAngles(self, q):
                q0, q1, q2, q3 = q
                roll = np.arctan2(2*(q0*q1 + q2*q3), 1 - 2*(q1**2 + q2**2))
                pitch = np.arcsin(2*(q0*q2 - q3*q1))
                yaw = np.arctan2(2*(q0*q3 + q1*q2), 1 - 2*(q2**2 + q3**2))
                angles = np.array([roll, pitch, yaw])
                return angles
        
        #calculates the derivatives of the state vector
        def dState(self, t, state):
                x,y,z,vx,vy,vz,q0,q1,q2,q3,omega_x,omega_y,omega_z = state

                X = np.array([x,y,z])
                V = np.array([vx,vy,vz])
                Q = np.array([q0,q1,q2,q3])
                W = np.array([omega_x,omega_y,omega_z])

                dV = self.rocket.forces / self.rocket.mass
                dX = V
                dQ = 0.5 * np.array([[0,-W[0], -W[1], -W[2]],
                                        [W[0], 0, W[2], -W[1]],
                                        [W[1], -W[2], 0, W[0]],
                                        [W[2], W[1], -W[0], 0]]) @ Q
                dW = np.linalg.inv(self.rocket.inertia) @ (self.rocket.momenta - np.cross(W, self.rocket.inertia @ W))

                return [dX[0],dX[1],dX[2],dV[0],dV[1],dV[2],dQ[0],dQ[1],dQ[2],dQ[3],dW[0],dW[1],dW[2]]
                
        def intertialRotate(self,theta_y,F_):
                r = R.from_euler('y', theta_y)
                F = r.apply(F_)
                return F
        
        def increment(self, t, s):
                x,z,theta_y,dx,dz,dtheta_y = s
                F = self.intertialRotate(theta_y,self.rocket.forces)
                ddx = F[0] / self.rocket.mass
                ddz = F[2] / self.rocket.mass - g
                ddtheta_y = self.rocket.momenta[1] / self.rocket.mass
                return [dx,dz,dtheta_y,ddx,ddz,ddtheta_y]

        def update(self, dt):
                newState = solve_ivp(self.increment, [0, dt], self.rocket.state)

                self.rocket.state = newState.y[:,-1]
                self.states = np.append(self.states, [self.rocket.state], axis=0)
                return self.rocket.state
        
        def plot(self):
                
                dim1 = [4*self.rocket.radius, 0, 0]
                dim2 = [0, 0, self.rocket.height]

                dim1_ = self.intertialRotate(self.rocket.state[2],dim1)
                dim2_ = self.intertialRotate(self.rocket.state[2],dim2)

                dims = np.abs(dim1_) + np.abs(dim2_)

                img_rotated = ndimage.rotate(img, -self.rocket.state[2]*180/np.pi, reshape=True)

                plt.imshow(img_rotated, extent=[-dims[0] + self.rocket.state[0],
                                        dims[0] + self.rocket.state[0], 
                                        -dims[2] + self.rocket.state[1], 
                                        dims[2] + self.rocket.state[1]])
                
                plt.plot(self.states[:,0],self.states[:,1],'--')

                #equal axes
                plt.axis('equal')
                plt.show()
                return plt