import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

from Quat import Quat

#img = mpimg.imread('falcon9.png')
class Plot3D:
    def __init__(self, rocket, states):
        self.rocket = rocket
        self.states = states

    def set(self,step,arrowlength):
        self.step = step
        self.arrowlength = arrowlength
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

    def origin(self):
        self.ax.plot(0,0,'o')

    def plot(self):
        self.ax = plt.figure(figsize=(20,10)).add_subplot(projection='3d')
        X, Y, Z = axes3d.get_test_data(0.05)


        self.ax.plot3D(self.states[:,0], self.states[:,1], self.states[:,2], 'red')
        self.ax.plot3D(self.states[:,0], self.states[:,1], 0*self.states[:,2], 'gray')
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')

        #make grid equal sized
        self.ax.set_aspect('equal')


    def attitude(self):
        for j in range(round(len(self.states)/self.step)):
                i = self.step*j
                state = self.states[i]
                quat = state[6:10]
                vec = Quat.intertialRotate(quat,[0,0,1])
                self.ax.quiver(state[0], state[1], state[2], vec[0], vec[1], vec[2], length=self.arrowlength, normalize=True, color='blue')

    def attitudeGraph(self):
        roll = []
        pitch = []
        yaw = []

        for j in range(round(len(self.states))):
                i = j
                state = self.states[i]
                quat = state[6:10]
                euler = Quat.eulerAngles2(quat)
                roll.append(euler[0])
                pitch.append(euler[1])
                yaw.append(euler[2])

        plt.plot(roll)
        plt.figure()
        plt.plot(pitch)
        plt.figure()
        plt.plot(yaw)
        plt.legend(['roll','pitch','yaw'])