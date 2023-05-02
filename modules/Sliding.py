import numpy as np

class Sliding:
    def __init__(self,eta,lam,phi,alpha,A,B):
        self.eta = eta
        self.lam = lam
        self.A = A
        self.B_inv = np.linalg.inv(B)
        self.phi = phi
        self.alpha = alpha
        F = np.diag([1,1,1])
        D = np.diag([1,1,1])
        self.K = F + D + self.eta

        self.ref = np.array([[0],[0],[0]])
        self.dref = np.array([[0],[0],[0]])
        self.ddref = np.array([[0],[0],[0]])

        self.inputs = np.array([[0],[0],[0]])
        """ self.inputs1 = []
        self.inputs2 = []
        self.t = []
        self.dt = 0.01
        self.lam = 20
        self.a_hat = np.array([[0],[0],[0],[0]])
        self.t_prev = 0 """
    
    """ def update_Y(self,q1,q2,dq1,dq2,t):
        dqr1,dqr2 = self.ddesired(t)
        ddqr1,ddqr2 = self.dddesired(t)
        self.Y = np.array([[ddqr1, ddqr2,(2*ddqr1+ddqr2)*np.cos(q2)-(dq2*dqr1 + dq1*dqr2 + dq2*dqr2)*np.sin(q2),(2*ddqr1+ddqr2)*np.sin(q2)+(dq2*dqr1 + dq1*dqr2 + dq2*dqr2)*np.cos(q2)],
                         [0, ddqr1 + ddqr2,ddqr1*np.cos(q2)+dq1*dqr1*np.sin(q2),ddqr1*np.sin(q2)-dq1*dqr1*np.cos(q2)]])
    
    def update_a_hat(self,q,dq,t):
        q_tilde = q - self.desired(t)
        dq_tilde = dq - self.ddesired(t)
        s = dq_tilde + self.lam*q_tilde

        dt = t - self.t_prev
        self.t_prev = t
        return -self.R@(self.Y).T@s*dt """
    
    def updateB(self,B):
        self.B_inv = np.linalg.inv(B)
        print(self.B_inv)
    
    def updateA(self,A):
        self.A = A

    def sat(self,s):
        for i in range(len(s)):
            if abs(s[i,0]) <= self.phi:
                s[i,0] = s[i,0]/self.phi
            else:
                s[i,0] = np.sign(s[i,0])
        return s

    def input(self,q,dq):
        q_tilde = q - self.ref
        dq_tilde = dq - self.dref
        s = dq_tilde + self.lam*q_tilde
        Q = -self.ddref + self.lam*(dq - self.dref)
        input = -self.B_inv@(self.K@self.sat(s) + self.A(dq) + Q*10)

        self.inputs = np.hstack((self.inputs,input))
        return input
    
    def phi_dot(self):
        return -self.alpha*self.phi + self.eta[0,0]
