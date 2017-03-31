#!/usr/bin/env python3
# Name: Mark Gustincic
# Group Members: Thomas Richards (tarichar)


import math
 
class ProteinCoordinates :
    """
 
    Author: David Bernick
    Date: March 21, 2013
    This class calculates angles and distances among a triad of points.
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r
 
    """
 
    def __init__(self,p,q,r) :
        """ Construct a Triad.  
            p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,0,0) ).
        """
        self.p = p
        self.q = q
        self.r = r
# private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b)))
    
    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))
    
    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r))
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) /    math.sqrt(self.d2(self.q,self.p)*self.d2(self.r,self.p)))
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /  math.sqrt(self.d2(self.p,self.q)*self.d2(self.r,self.q)))
 
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /  math.sqrt(self.d2(self.p,self.r)*self.d2(self.q,self.r)))
   
    def main():

        """
        Prompts the user to enter three coordinates, and splits the inputted string at each ")" symbol.
        Assigns value to string in coords1[0] one position after the "(" symbol. This string is
        split again after each "," symbol, resulting in a list of three float values Cx, Cy, Cz.
        Object name Carbon is assigned to Cx, Cy, and Cz.
        """
        coordsX = input("Please enter three atomic coordinates: ")
        coords1 = coordsX.split(")")

        """
        Assigns value to string in coords1[0] one position after the "(" symbol. This string is
        split again after each "," symbol, resulting in a list of three float values Cx, Cy, Cz.
        Object name Carbon is assigned to Cx, Cy, and Cz.
        """
        start = coords1[0].find('(')+1
        C = coords1[0][start:]
        CNew = C.split(",")
        Cx = float(CNew[0])
        Cy = float(CNew[1])
        Cz = float(CNew[2])
        Carbon = (Cx,Cy,Cz)
        """
        Assigns value to string in coords1[1] one position after the "(" symbol. This string is
        split again after each "," symbol, resulting in a list of three float values Nx, Ny, Nz.
        Object name Nitrogen is assigned to Nx, Ny, and Nz.
        """
        start = coords1[1].find('(')+1
        N = coords1[1][start:]
        NNew = N.split(",")
        Nx = float(NNew[0])
        Ny = float(NNew[1])
        Nz = float(NNew[2])
        Nitrogen = (Nx,Ny,Nz)
        """
        Assigns value to string in coords1[2] one position after the "(" symbol. This string is
        split again after each "," symbol, resulting in a list of three float values Cax, Cay, Caz.
        Object name Carbon is assigned to Cax, Cay, and Caz.
        """  
        start = coords1[2].find('(')+1
        Ca = coords1[2][start:]
        CaNew = Ca.split(",")
        Cax = float(CaNew[0])
        Cay = float(CaNew[1])
        Caz = float(CaNew[2])
        Calcium = (Cax,Cay,Caz)
        """
        Instantiates name1 as instance of ProteinCoordinates with arguments Carbon, Nitrogen, and Calcium.
        Names are assigned to methods used with name1, and these names are formatted/printed to output. Main
        method is called.
        """
        name1 = ProteinCoordinates(Carbon, Nitrogen, Calcium)
        C_N = name1.dPQ()
        N_Ca = name1.dQR()
        C_N_Ca = (name1.angleQ()*180)/math.pi
        print("N-C bond length = %0.2f" % C_N)
        print("N-Ca bond length = %0.2f" % N_Ca)
        print("C-N-Ca bond angle = %0.1f" % C_N_Ca)
    
ProteinCoordinates.main()
        













        
        
