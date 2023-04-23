import numpy as np
from scipy.spatial.transform import Rotation as R

class Quat:
    #Angle and quaternion functions
    @staticmethod
    def eulerToQuat(euler):
        r = R.from_euler('xyz', euler, degrees=False)
        return r.as_quat()
    
    #converts quaternion to euler angles
    @staticmethod
    def eulerAngles(q):
        q = q / np.linalg.norm(q)
        q0, q1, q2, q3 = q
        roll = np.arctan2(2*(q0*q1 + q2*q1), 1 - 2*(q1**2 + q2**2))
        pitch = np.arcsin(2*(q0*q2 - q3*q1))
        yaw = np.arctan2(2*(q0*q3 + q1*q2), 1 - 2*(q2**2 + q3**2))
        angles = np.array([roll, pitch, yaw])
        return angles
    
    @staticmethod
    def intertialRotate(q, vector):
        qq = np.array([q[2],q[1],q[0],q[3]])
        r = R.from_quat(qq)
        return r.apply(vector)
    
    @staticmethod
    def intertialRotateInv(q, vector):
        qq = np.array([q[2],q[1],q[0],q[3]])
        r = R.from_quat(qq)
        r = r.inv()
        return r.apply(vector)

    #converts quaternion to euler angles using scipy
    @staticmethod
    def eulerAngles2(q):
        qq = np.array([q[2],q[1],q[0],q[3]])
        r = R.from_quat(qq)
        return r.as_euler('xyz', degrees=False)