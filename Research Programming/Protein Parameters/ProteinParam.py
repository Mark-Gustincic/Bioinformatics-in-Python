#!/usr/bin/env python3
# Name: Mark Gustincic (mgustinc)
# Group Members: None

class ProteinParam:
 
    """
    The ProteinParam program is designed to take in a string of amino acids and
    return certain properties associated with the protein. Such properties
    include the number of amino acids, the protein's molecular weight, the molar
    extinction coefficient, the mass extinction coefficient, the theoretical
    isoelectric point, and the protein's amino acid composition expressed as
    percentages.
    
    The below section of code initializes dictionaries and value names as class
    attributes for use in methods within the ProteinParam class.
    """
    
    aa2mw = {
       'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
       'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
       'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
       'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
       }
    
    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

# the __init__ method requires a protein string to be provided, either as a
# string or list of strings that will be concatenated

    def __init__ (self, protein):
         """
         The init method cleans the input string and stores the aa2mw
         dictionary. Object l removes whitespace in inputted protein string,
         then splits the string into individual amino acids.

         The if statement within the for loop looks in self.protString for objects
         that are not included in the aa2mw dictionary. If such objects exist,
         they are removed with the 'replace' parsing function.

         The aa2mw dictionary is instantiated here for use in the aaComposition
         method.

         The aaDictionary is initialized here with amino acids as keys and their
         counts in the protein sequence as the associated values. The for loop
         counts each valid amino acid in sequence, and updates the value of each
         amino acid key with this count.
         """
         l = ''.join(protein).split()
         self.protString = ''.join(l).upper()
         k = self.protString
         
         for x in k:
            if x not in self.aa2mw:
                 k = k.replace(x,'')
         self.protString = k

         self.aa2mw = {

         'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
         'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
         'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
         'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
         }

         self.aaDictionary = {}
         for x in self.aa2mw.keys():
            count = self.protString.count(x)
            self.aaDictionary.update({x:count})


    def aaCount (self):
        """
        The aaCount method returns the sum of all valid amino acids included
        in self.protString.
        """
        
        totalAA = sum(x != '' for x in self.protString)
        return totalAA

    def pI (self):
        """
        The pI method returns the theoretical isoelectric point of the inputted
        protein sequence. This value is estimated by finding the pH that yields
        a neutral net charge. This pH is found through iteration of all pH values
        in order to find the one closest to zero.
        """
        a = 0
        b = 14
        middle = 0
        if self.protString != '':
            while (b-a) > .01:
                middle = (a + b)/2
                if self._charge_(middle) > 0:
                      a = middle
                if self._charge_(middle) < 0:
                      b = middle
        return middle

    def aaComposition (self):
        """
        The method aaComposition returns a dictionary of amino acid counts in the
        sequence. This dictionary is initialized in the init method.
        """
        count = self.aaDictionary
        return count
 

    def _charge_ (self, pH):
        """
        The _charge_ method returns the net charge on the protein at a specified
        pH. The charge on the N terminus is calculated by divding the N-terminus
        acidity constant by the sum of the the acidity constant and the hydrogen
        ion concentration. The charge on the C terminus is calculated by dividing
        the hydrogen ion concentration by the sum of the C terminus acidity
        constant and the hydrogen ion concentration. The final charge is calculated
        by subtracting the C terminus charge from the N terminus charge.
        """
        aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
        aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
        aaNterm = 9.69
        aaCterm = 2.34
        chargeN = 10**aaNterm/(10**aaNterm + 10**pH)
        chargeC = 10**pH/(10**aaCterm + 10**pH)
        charge = chargeN - chargeC
        """
        The for loop below is used to calculate the protein charge when charged
        amino acids are included in sequence. Charge is calculated by dividing
        the amino acid's acidity constant by the sum of the specific acidity
        constant and the hydrogen ion concentration.
        
        The charges of (+) amino acids are
        added to the original charge value, while the charges of (-) amino acids
        are subtracted from the original charge value.
        """       
        for a in self.protString:
            if a == 'K' and 'R' and 'H':
                 b = aa2chargePos.get(a)
                 charge += 10**b/(10**b + 10**pH)
            if a == 'D' and 'E' and 'C' and 'Y':
                 b = aa2chargeNeg.get(a)
                 charge -= 10**pH/(10**b + 10**pH)
        return charge
        

    def molarExtinction (self):
         """
         The method molarExtinction returns an estimation of the protein's
         extinction coefficient. This value indicates the amount of light the
         protein absorbs at a wavelength of 280 nm. Tyrosine's, tryptophan's,
         and cystine's molar extinction coefficient are multiplied by their
         relative counts in the amino acid sequence, then summed. 
         """
         aa2abs280 = {'Y':1490, 'W': 5500, 'C': 125}
         num_Y = self.protString.count('Y')
         num_W = self.protString.count('W')
         num_C = self.protString.count('C')
         molar = num_Y*aa2abs280['Y'] + num_W*aa2abs280['W'] + num_C*aa2abs280['C']
         return molar
    
    def massExtinction (self):
        """
        The method massExtinction returns the mass extinction coefficient. This
        value is calculated by dividing the molar extinction coefficient (found
        in the molarExtinction method) by the molecular weight value (obtained
        in the molecularWeight method).
        """
         
        MolWeight =  self.molecularWeight()
        unknown = 0.0
        if MolWeight <= 0.0:
            massExt = unknown
        else:
            massExt = self.molarExtinction() / MolWeight
        return massExt 
     

    def molecularWeight (self):
        """
        The method molecularWeight returns the calculated molecular weight value
        of the inputted protein. The number of all H2Os created by peptide bond
        formation in sequence are summed and multiplied by H2O's molecular weight.
        This value is subtracted from 0 and initialized as weight.
        
        The for loop iterates through all amino acids in sequence, and adds each
        amino acid's molecular weight to object weight's value. The if statement
        at the very end of the method handles molecular weight exceptions
        so that any non-protein sequence inputted will return a value of 0 instead
        of 18.
        """

        mwH2O = 18.015
        totalH2O = (len(self.protString)-1)
        tw = totalH2O*mwH2O
        w = 0 - tw
        
        for x in self.protString:
             num_aa = self.aa2mw[x]
             w += num_aa
        if w <= 18.015:
            w = 0
        return w
            
        
# Please do not modify any of the following.  This will produce a standard output that can be parsed

import sys
print("Input an amino acid sequence:")
for inString in sys.stdin :

    myParamMaker = ProteinParam(inString)
    c = 0
    
    boolean = False
    while not boolean:
        c += 1

        myAAnumber = myParamMaker.aaCount()

        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))

        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))

        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))

        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))

        print ("Amino acid composition:")

        myAAcomposition = myParamMaker.aaComposition()

        keys = list(myAAcomposition.keys())

        keys.sort()
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 

        for key in keys :

            print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))
        if c > 0:
            break
    prompt = input('Analyze another amino acid sequence?(yes/no)\n').lower()
    if prompt == 'yes':
        boolean = False
        print('New sequence:')
    else:
        sys.exit()
