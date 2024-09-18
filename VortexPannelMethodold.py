import json
import numpy as np

class VortexPannelMethod:
    def __init__(self):
        pass
        # opens json file and sets airfoil text file as an x and y array
    def get_airfoil_array(self,filename):
        json_string= open(filename).read()
        #print json string
        input_dict=json.loads(json_string)
        # set up dictionary
        NACAfile= input_dict["airfoils"]
        # set up empty array
        filename = input_dict["airfoils"][0]
        with open(filename,"r") as f:
            points = [ [float(i) for i in line.split()] for line in f.readlines() ]
        # convert single aray into x y array
        xy = np.array(points)
        x=xy[:,0]
        y=xy[:,1]
        return x,y

    # get velocity from json file
    def get_velocity(self,filename):

        json_string=open(filename).read()
        #print json as string
        input_dict=json.loads(json_string)
        #open dict and retreave velocity
        velocity = input_dict["velocity"]
        #return velocity
        return velocity

    # get angle of atack in degrees
    def get_angle_of_atack(self,filename):
        json_string=open(filename).read()
        #print json as string
        input_dict=json.loads(json_string)
        #open dict and retreave alpha
        alpha = input_dict["alpha[deg]"]
        #return velocity
        return alpha

    # finds control points from an array of x and y values(4.21)
    def get_control_points(self,x,y):
        # equation (4.21)
        x_cp=(x[:-1]+x[1:])/2
        y_cp=(y[:-1]+y[1:])/2
        return x_cp, y_cp

    # calculates the length of the jth pannel(4.23)
    def get_length_of_jth_pannel(self,x,y,j):
        
        # sets up empty array
        #calculates length of the pannel saves as temp variable
        length_pannel=np.sqrt(np.power(x[j+1]-x[j],2)+np.power(y[j+1]-y[j],2))
        return length_pannel

    #calculates the xi and eta matrix(4.24)
    def get_xi_eta(self,x,y,x_cp,y_cp,l_j,i,j):
        xi_eta=[]
        #set up first matrix
        matrix1=[[(x[j+1]-x[j]), (y[j+1]-y[j])],
            [-(y[j+1]-y[j]), (x[j+1]-x[j])]]
        #set up second matrix
        matrix2=[[x_cp[i]-x[j]],
            [y_cp[i]-y[j]]]
        #multiply matrix 1 with matrix 2
        xi_eta=(1/l_j)*np.matmul(matrix1,matrix2)
        xi=xi_eta[0]
        eta=xi_eta[1]
        return xi[0],eta[0]

    # calculate phi
    def get_phi(self,eta,xi,l_j):
        phi=np.arctan2(eta*l_j,np.power(eta,2)+np.power(xi,2)-xi*l_j)
        return phi

    # calculate psi
    def get_psi(self,eta,xi,l_j):
        psi=(1/2)*np.log(np.abs((np.power(xi,2)+np.power(eta,2))/(np.power((xi-l_j),2)+np.power(eta,2))))
        return psi

    # get P matrix
    def get_P_matrix(self,x,y,x_cp,y_cp,i,j):
        # get length of jth pannel
        l_j=self.get_length_of_jth_pannel(x,y,j)
        # get xi and eta
        xi,eta=self.get_xi_eta(x,y,x_cp,y_cp,l_j,i,j)
        # gets phi
        phi=self.get_phi(eta,xi,l_j)
        # gets psi
        psi=self.get_psi(eta,xi,l_j)
        #sets up matrix 1
        matrix_1=[[(x[j+1]-x[j]), -1*(y[j+1]-y[j])],
        [(y[j+1]-y[j]), (x[j+1]-x[j])]]
        #setup matrix 2
        matrix_2=[[((l_j-xi)*phi+eta*psi), (xi*phi-eta*psi)],
        [(eta*phi-(l_j-xi)*psi-l_j), (-eta*phi-xi*psi+l_j)]]
        #calculate matrix mul for m1 and m2
        P_temp=np.matmul(matrix_1,matrix_2)
        #calculate random constant
        p_const=1/(2*np.pi*np.power(l_j,2))
        #calculate full p matrix
        P=p_const*P_temp
        #return P matrix
        return P

    #calculates A matrix
    def get_A_matrix(self,x,y,x_cp,y_cp):
        #determine number of points
        n=len(x)
        #initialize A matrix
        A=np.zeros((n,n))
        # calculate A Matrix
        for i in range(n-1):
            for j in range(n-1):
                #get P matrix
                p=self.get_P_matrix(x,y,x_cp,y_cp,i,j)
                #get length of pannel with i index
                l_i=self.get_length_of_jth_pannel(x,y,i)
                #calculate A (4.29)
                A[i,j]=A[i,j]+((x[i+1]-x[i])/l_i)*p[1,0]-((y[i+1]-y[i])/l_i)*p[0,0]
                A[i,j+1]=A[i,j+1]+((x[i+1]-x[i])/l_i)*p[1,1]-((y[i+1]-y[i])/l_i)*p[0,1]
        #set index to 1
        A[n-1,0]=1.0
        A[n-1,n-1]=1.0
        return A

    def get_b_matrix(self,x,y,filename):
        #determine number of points
        n=len(x)
        # initialize B matrix
        B=np.zeros((n,1))
        # get angle of atack
        alpha=self.get_angle_of_atack(filename)
        # calculates B matrix(4.32)
        for i in range(n-1):
            l_j=self.get_length_of_jth_pannel(x,y,i)
            B[i,0]=((y[i+1]-y[i])*np.cos(alpha*np.pi/180)-(x[i+1]-x[i])*np.sin(alpha*np.pi/180))/l_j
        B[n-1,0]=0.0
        #get velocity
        v=self.get_velocity(filename)
        B=v*B
        return B

    #get gamma
    def get_gamma_vector(self,A,B):
        #solve for gamma(4.32)
        gamma=np.linalg.solve(A,B)
        return gamma

    #get lift coefficient
    def get_CL(self,x,y,gamma,filename):
        n=len(x)
        c=x[0]
        #set cl to zero
        C_l=0
        #calculate sumation(4.36)
        for i in range(n-1):
            l_i=self.get_length_of_jth_pannel(x,y,i)
            v=self.get_velocity(filename)
            C_l+=(l_i/c)*((gamma[i]+gamma[i+1])/v)
        return C_l

    # get pitching moment at leading edge
    def get_cmle(self,x,y,gamma,filename):
        n=len(x)
        c=x[0]
        # set piching moment at leading edge to zero
        cmle_temp=0
        # Calculate sumation(4.37)
        for i in range(n-1):
            l_i=self.get_length_of_jth_pannel(x,y,i)
            alpha=self.get_angle_of_atack(filename)
            v=self.get_velocity(filename)
            cmle_temp+=(l_i/c)*(((2*x[i]*gamma[i]+x[i]*gamma[i+1]+x[i+1]*gamma[i]+2*x[i+1]*gamma[i+1])/(c*v))*np.cos(alpha*np.pi/180)+((2*y[i]*gamma[i]+y[i]*gamma[i+1]+y[i+1]*gamma[i]+2*y[i+1]*gamma[i+1])/(c*v))*np.sin(alpha*np.pi/180))
        # multiply total by (-1/3)
        Cmle=-(1/3)*cmle_temp
        return Cmle

    # get pitching moment at quarter chord
    def get_cmc4(self,C_l,C_mle,filename):
        alpha=self.get_angle_of_atack(filename)
        C_mc4=C_mle+C_l/4*np.cos(alpha*np.pi/180)
        return C_mc4

    # run function
    def run(self,filename):
    #get x y points
        x,y=self.get_airfoil_array(filename)

        #get control points
        x_cp,y_cp=self.get_control_points(x,y)

        #get A matrix
        A=self.get_A_matrix(x,y,x_cp,y_cp)

        #get B matrix
        B=self.get_b_matrix(x,y,filename)

        #get gamma
        gamma=self.get_gamma_vector(A,B)

        # get lift coefficient
        Cl=self.get_CL(x,y,gamma,filename)

        # get pitching moment at leadig edge
        CmLE=self.get_cmle(x,y,gamma,filename)

        # get pitching mometn at quarter chord
        Cmc4= self.get_cmc4(Cl,CmLE, filename)
        print(Cl[0])
        print(CmLE[0])
        print(Cmc4[0])

#main
if __name__=="__main__":
    q=VortexPannelMethod()
    q.run("input.json")
    