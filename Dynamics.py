from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation as R
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

img = mpimg.imread('falcon9.png')
g = 9.81

class Dynamics:
        def __init__(self, rocket):
                self.rocket = rocket
                self.state = np.zeros(13)
                self.state[9] = 1 #initial quaternion for vertical position
                self.states = np.array([self.state])

                self.externalAccelerations = np.array([0,0,-9.81])
        
        #calculates the derivatives of the state vector
        def dState(self, t, state):
                x,y,z,vx,vy,vz,q0,q1,q2,q3,omega_x,omega_y,omega_z = state
                self.rocket.sate = state
                self.rocket.update()

                X = np.array([x,y,z])
                V = np.array([vx,vy,vz])
                Q = np.array([q0,q1,q2,q3])
                #normalize quaternion
                Q = Q / np.linalg.norm(Q)
                W = np.array([omega_x,omega_y,omega_z])

                dX = V
                dV = self.intertialRotate(Q,self.rocket.F) / self.rocket.mass + self.externalAccelerations
                dW = np.linalg.inv(self.rocket.inertia) @ (self.rocket.M - np.cross(W, self.rocket.inertia @ W))
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
        

        #Angle and quaternion functions
        def eulerToQuat(self, euler):
                r = R.from_euler('xyz', euler, degrees=False)
                return r.as_quat()
        
        def getQuat(self):
                return self.state[6:10]
        
        def getEuler(self):
                return self.eulerAngles2(self.getQuat())
        
                #converts quaternion to euler angles
        def eulerAngles(self, q):
                q = q / np.linalg.norm(q)
                q0, q1, q2, q3 = q
                roll = np.arctan2(2*(q0*q1 + q2*q1), 1 - 2*(q1**2 + q2**2))
                pitch = np.arcsin(2*(q0*q2 - q3*q1))
                yaw = np.arctan2(2*(q0*q3 + q1*q2), 1 - 2*(q2**2 + q3**2))
                angles = np.array([roll, pitch, yaw])
                return angles
        
        def intertialRotate(self, q, vector):
                qq = np.array([q[2],q[1],q[0],q[3]])
                r = R.from_quat(qq)
                return r.apply(vector)

        #converts quaternion to euler angles using scipy
        def eulerAngles2(self, q):
                qq = np.array([q[2],q[1],q[0],q[3]])
                r = R.from_quat(qq)
                return r.as_euler('xyz', degrees=False)
        
        #Constraint functions
        def ground(self,state,newState):
                if state[2] < 0:
                        newState[0] = 0
                        newState[1] = 0
                        newState[2] = 0


        #Plotting functions
        def plot(self):
                
                dim1 = [4*self.rocket.radius, 0, 0]
                dim2 = [0, 0, self.rocket.height]

                dim1_ = self.intertialRotate(self.state[6:10],dim1)
                dim2_ = self.intertialRotate(self.state[6:10],dim2)

                dims = np.abs(dim1_) + np.abs(dim2_)
                angle = self.getEuler()[1]
                img_rotated = ndimage.rotate(img, -angle, reshape=True)

                plt.imshow(img_rotated, extent=[-dims[0] + self.state[0],
                                        dims[0] + self.state[0], 
                                        -dims[2] + self.state[2], 
                                        dims[2] + self.state[2]])
                
                plt.plot(self.states[:,0],self.states[:,2],'--')

                #equal axes
                plt.axis('equal')
                plt.show()
                return plt
        
        
        def plot3d(self):
                ax = plt.figure(figsize=(20,10)).add_subplot(projection='3d')
                X, Y, Z = axes3d.get_test_data(0.05)


                ax.plot3D(self.states[:,0], self.states[:,1], self.states[:,2], 'red')
                ax.plot3D(self.states[:,0], self.states[:,1], 0*self.states[:,2], 'gray')
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                ax.set_zlabel('Z')
                ax.set_aspect('equal')

                """ for vec in self.rocket.inputs:
                        ax.quiver(vec[0], vec[1], vec[2], vec[3], vec[4], vec[5], color='b', arrow_length_ratio=0.1)
                ax.quiver(0, 0, 0, 1, 1, 1, color='b', arrow_length_ratio=0.1)
 """
                #change perspective
                #ax.view_init(view[0], view[1])

                #equal axis lengths

        