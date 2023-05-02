import numpy as np

class Adaptive:
    def __init__(self,eta,lam,phi,alpha,A,B):
        self.eta = eta
        self.lam = lam
        self.A = A
        self.B_inv = np.linalg.inv(B)
        self.phi = phi
        self.alpha = alpha
        self.K = self.eta

        self.ref = np.array([[0],[0],[0]])
        self.dref = np.array([[0],[0],[0]])
        self.ddref = np.array([[0],[0],[0]])

        self.inputs = np.array([[0],[0],[0]])
        self.states = np.array([[0],[0],[0]])
        self.dstates = np.array([[0],[0],[0]])

    def sigmoid(self,x):
        return 1/(1+np.exp(-x))
                  
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

        
        input = -self.B_inv@(self.K@self.sat(s) + self.A(q,dq) + Q)
        print(input)
        self.inputs = np.hstack((self.inputs,input))
        self.states = np.hstack((self.states,q))
        self.dstates = np.hstack((self.dstates,dq))
        return input
    
    def phi_dot(self):
        return -self.alpha*self.phi + self.eta[0,0] + 10