#!usr/bin/env python3
# Name: Mark Gustincic (mgustinc) 
# Group Members: Thomas Richards (tarichar)

'''
input: fasta file or standard in containing template sequence. Primer sequences entered to standard in.

output: Primer statistics
        Amplicon sequence & statistics
        Restriction enzyme locator
        Specific restriction enzyme fragment calculator

        Protocol for running PCR
        Protocol for running gel electrophoresis
'''

import sys
class FastAreader :
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''       
        self.fname = fname
            
    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
 
    def readFasta (self):	
        header = ''
        sequence = ''       
        with self.doOpen() as fileH:			
            header = ''
            sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()
            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
						
        yield header,sequence

class PCR:
    """
    This class models what happens to a DNA sequence when it undergoes PCR using two primers to amplify a subsequence of the template DNA.
    Once the amplicon is determined, restriction enzyme cutsites are located within the sequence. The user is given the option to see fragment
    lengths created by specific enzymes within the sequence.
    """
    def __init__(self, seq):
        """
        The init method initializes lists and dictionaries to be used in other methods of the PCR class. The bases list contains a list of the four nucleotide bases
        found in DNA and is used for cleaning the input sequence. The self.compDict dictionary contains the complementary base pairs for each nucleotide, and is used
        to produce complentary sequences of DNA. The self.dnaMass dictionary contains the molecular weight of each nucleotide to be used to calculate the molecular
        weight of a DNA sequence. self.REdict dictionary contains enzyme name:cut sequence as a key value pair to be used to find cutsites within the amplicon. The code
        below the self.REdict dictionary is responsible for cleaning the inputted sequence so that it only contains the four bases present in DNA.

        Pseudo Code:
        bases = a list of the four nucleotide bases found in DNA
        self.compDict = a dictionary containing nucleotide basepairs as key:value pairs
        self.primerDict = a dictionary with forward and reverse as keys with values initialized at zero
        self.dnaMass = a dictionary containing the molecular weight of each nucleotide (the weight of its base pair is included to account for complementary strand
        self.REdict = a dictionary containg restriction enzymes and their cutsites as key:value pairs

        sequence with unwanted characters (rawSeq) is joined at all whitespaces, and then split into individual characters
        self.seq = join rawSeq at all whitespaces and capitalized
        cleanSeq = self.seq

        for character in cleanSeq:
            if character is not a nucleotide:
                character is replaced by whitespace
        self.seq is saved after removing non-nucleotide characters
        
        """
        bases = ['A', 'T', 'C', 'G']
        self.compDict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
        self.primerDict = {'f':0, 'r':0}
        self.dnaMass={'G': 652.41, 'C': 652.41, 'A': 651.43, 'T': 651.43}
        self.REdict = {'AatII':'GACGT|C', 'Acc65I':'G|GTACC', 'AccI': 'GT|AGAC', 'AccI': 'GT|ATAC', 'AccI': 'GT|CGAC', 'AccI': 'GT|CTAC',
                       'AclI': 'AA|CGTT', 'Afel': 'AGC|GCT', 'AflII': 'C|TTAAG', 'Agel': 'C|TTAAG', 'ApaI': 'GGGCC|C', 'ApaLI': 'GTGCAC',
                       'Apol': 'A|AATTC', 'Apol': 'G|AATTT', 'Apol': 'A|AATTT', 'Apol': 'G|AATTC', 'Ascl': 'GG|CGCGCC', 'AseI': 'AT|TAAT', 'AsiSI': 'GCGAT|CGC', 'AvrII': 'C|CTAGG',
                       'BamHI': 'G|GATCC', 'BclI': 'T|GATCA', 'BglII': 'A|GATCT', 'Bme1580I': 'GGGCA|C', 'Bme1580I': 'GTGCC|C', 'Bme1580I': 'GTGCA|C', 'Bme1580I': 'GGGCC|C', 'BmtI': 'GCTAG|C',
                       'BsaHI': 'GA|CGCC', 'BsaHI': 'GG|CGTC', 'BsaHI': 'GG|CGCC', 'BsaHI': 'GA|CGTC', 'BsiEI': 'CGAC|CG', 'BsiEI': 'CGGT|CG', 'BsiEI': 'CGGC|CG', 'BsiEI': 'CGAT|CG',
                       'BsiWI': 'C|GTACG', 'BspEI': 'T|CCGGA', 'BspHI': 'T|CATGA','BsrGI': 'T|GTACA', 'BssHII': 'G|CGCGC', 'BstBI': 'TT|CGAA', 'BstZ17I': 'GTA|TAC',
                       'BtgI': 'C|CACGG', 'BtgI': 'C|CGTGG', 'BtgI': 'C|CATGG', 'BtgI': 'C|CGCGG',
                       'ClaI': 'AT|CGAT', 'DraI': 'TTT|AAA', 'EaeI': 'C|GGCCA', 'EaeI': 'T|GGCCG', 'EaeI': 'T|GGCCA', 'EaeI': 'C|GGCCG',
                       'EagI': 'C|GGCCG', 'EcoRI': 'G|AATTC', 'EcoRV': 'GAT|ATC', 'FseI': 'GGCCGG|CC', 'FspI': 'TGC|GCA',
                       'HaeII': 'AGCGC|T', 'HaeII': 'GGCGC|C', 'HaeII': 'GGCGC|T', 'HaeII': 'AGCGC|C',
                       'HincII': 'GTC|AAC', 'HincII': 'GTT|GAC', 'HincII': 'GTC|GAC', 'HincII': 'GTT|AAC',
                       'HindIII': 'A|AGCTT', 'HpaI': 'GTT|AAC', 'KasI': 'G|GCGCC', 'KpnI': 'GGTAC|C',
                       'MfeI': 'C|AATTG', 'MluI': 'A|CGCGT', 'MscI':'TGG|CCA', 'MspA1I': 'CAG|CGG', 'MspA1I': 'CCG|CTG','MspA1I': 'CAG|CTG','MspA1I': 'CCG|CGG',
                       'NaeI': 'GCC|GGC', 'NarI': 'GG|CGCC', 'NcoI' : 'C|CATGG', 'NdeI': 'CA|TATG', 'NcoI': 'C|CATGG', 'NgoMIV': 'G|CCGGC', 'NheI': 'G|CTAGC',
                       'NotI': 'GC|GGCCGC', 'NruI': 'TCG|CGA', 'NsiI': 'ATGCA|T', 'NspI': 'ACATG|C', 'NspI': 'GCATG|T', 'NspI': 'ACATG|T', 'NspI': 'GCATG|C',
                       'PacI': 'TTAAT|TAA', 'PciI': 'A|CATGT', 'PmeI': 'GTTT|AAAC', 'PmlI': 'CAC|GTG', 'PsiI':	'TTA|TAA', 'PspOMI': 'G|GGCCC', 'PstI':	'CTGCA|G',
                       'PvuI': 'CGAT|CG', 'PvuII': 'CAG|CTG', 'SacI': 'GAGCT|C', 'SacII': 'CCGC|GG', 'SalI': 'G|TCGAC', 'SbfI': 'CCTGCA|GG', 'ScaI': 'AGT|ACT',
                       'SfcI': 'C|TACAG', 'SfcI': 'C|TGTAG', 'SfcI': 'C|TATAG', 'SfcI': 'C|TGCAG', 'SfoI': 'GGC|GCC', 'SgrAI': 'CA|CCGGCG', 'SgrAI': 'CG|CCGGTG',
                       'SgrAI': 'CA|CCGGTG', 'SgrAI': 'CG|CCGGCG', 'SmaI': 'CCC|GGG', 'SmlI': 'C|TCAAG', 'SmlI': 'C|TTGAG', 'SmlI': 'C|TCGAG', 'SmlI': 'C|TTAAG',
                       'SnaBI':	'TAC|GTA', 'SpeI': 'A|CTAGT', 'SphI': 'GCATG|C', 'SspI': 'AAT|ATT', 'StuI': 'AGG|CCT', 'SwaI': 'ATTT|AAAT', 'XbaI': 'T|CTAGA',
                       'XhoI': 'C|TCGAG', 'XmaI': 'C|CCGGG'
                       }
            

        rawSeq = ''.join(seq).split()
        self.seq = ''.join(rawSeq).upper()
        cleanSeq = self.seq

        for char in cleanSeq:
            if char not in bases:
                cleanSeq = cleanSeq.replace(char,'')
        self.seq = cleanSeq

    def amplicon(self):
        """
        The amplicon method is responsible for creating the virtual PCR product given the two inputted primers. If there is no PCR product, a string is
        outputted stating that there is no match for the inputted primers. If a PCR product does exist, the user is prompted if they want the amplicon
        printed. The molecular weight and GC content of the amplicon are calculated and printed. Exceptions are raised for calculating GC content to
        prevent a DivisionByZero error if no amplicon exists.

        Pseudo code:
        set name forward to value of the key f in self.primerDict
        set name reverse to value of the key r in self.primerDict
        product = pcr product found in pcrProd method

        if there is no product:
            print statement notifying user
            return -1
        prompt user if both primers are located on same strand
        if yes:
            print product
        count A, G, C, and T in product
        try:
            find GC content of product
        except:
            GC product = 0  
        calculate product molecular weight in kDa
        try:
            print GC content to output
        except:
            GC product = 0
        """
        forward = self.primerDict['f']
        reverse = self.primerDict['r']
        product = self.pcrProd()
        if len(product) == 0:
            print('There are no matches to {0} and {1} in the template sequence.'.format(forward, reverse))
            return -1
        print('')
        print('Amplicon sequence found.')
        prompt = input('Show sequence? (yes/no) \n')
        print('')
        if prompt == 'yes':
            print('Expected product sequence: \n5 {0} 3\n'.format(product))
        ac = product.count('A')
        gc = product.count('G')
        tc = product.count('T')
        cc = product.count('C')
        try:
            gcContent = ((gc+cc)/(gc+cc+ac+tc))*100
            tryb = True
        except ZeroDivisionError:
            tryb = False
        ampWeight = (ac*651.43)+(gc*652.41)+(cc*652.41)+(tc*651.43)
        ampWeight = ampWeight/1000
        print('Amplicon length: {0} base pairs'.format(len(product)))
        print('Molecular weight: {0} kDa'.format(round(ampWeight,2)))
        if tryb:
            print('GC content: {0}%\n'.format(round(gcContent,2)))
        else:
            return 0
            

    def protocol(self):
        """
        The protocol method contains a list of dynamic instructions for running PCR in a lab setting. The user is prompted to run PCR, and if yes is inputted, the user is prompted
        again for the desired reaction volume. This volume is saved and implemented in following steps of the protocol to produce dynamic instructions. An exception is raised when
        prompting the user for the volume to account for non-numbers being inputted. Once volume is obtained, it is converted to a float to be used in the print statements below.
        The PCR protocol is printed out line by line, requiring the user to press enter for the next step to be printed. This helps the user keep track of where they are in the
        protocol.

        Pseudo Code:
        prompt user to run PCR
        if user returns yes:
            prompt user for desired reaction volume
            while true:
                try:
                    make volume a float
                    break
                except ValueError:
                    ask user to input a number
            if volume is equal to '' or is less than 0:
                volume = input(prompt user for a positive number
            make volume into a float
            print instructions for running PCR line by line so that user must press enter for next step to be printed.
            some print statements are dynamic, and output is dependent on inputted value for volume
        """
        prompt = input('\nWould you like to run PCR now? (yes/no) \n')
        if prompt == 'yes':
            volume = input('What is the desired reaction volume in uL?\n')
            while True:
                try:
                    volumef = float(volume)
                    break
                except ValueError:
                    volume = input('Please enter a number\n')
                    continue

            if volume == '' or volumef <= 0:
                volume = input('Please enter a positive number\n')
                
            volumef = float(volume)
            print('\nThe following instructions demonstrate how to amplify DNA sequences using polymerase chain reaction(PCR).\nPress enter to proceed to the next step.')
            a = input('')
            print('Ingredients needed: \n    1.deionized Water\n    2.Taq polymerase buffer\n    3.forward primer\n    4.reverse primer\n    5.dNTPs\n    6.Taq polymerase\n    7.template DNA\n')
            print('If running multiple reactions, multiply given volumes by number of reactions desired.\nSet up PCR master mix in a snap-cap tube. Keep all reaction components and master mix on ice during preparation.'.format(volume))
            c = input('')
            print('This guide is for a {0} uL PCR reaction.\nTo master mix:\n1. Add {1} uL of deionized Water'.format(volume,round(volumef*.565, 2)))
            d = input('')
            print('2. Add {0} uL of 10X Taq polymerase buffer.'.format(round(volumef*.1, 2)))
            e = input('')
            print('3. Add {0} uL of forward primer concentrated at 10 uM.'.format(round(volumef*.02, 2)))
            f = input('')
            print('4. Add {0} uL of reverse primer concentratated at 10 uM.'.format(round(volumef*.02, 2)))
            g = input('')
            print('5. Add {0} uL of dNTPs concentrated at 10 mM.'.format(round(volumef*.02, 2)))
            h = input('')
            print('6. Add {0} uL of Taq polymerase \nIMPORTANT: enzyme is tempterature sensitive. Keep on ice to prevent denaturation.'.format(round(volumef*.025, 2)))
            i = input('')
            print('Aliquot {0} uL of master mix into thin-walled PCR tubes.'.format(round(volumef*.75, 1)))
            j = input('')
            print('Add {0} uL of template DNA to each PCR tube, bringing total reaction volume to {1} uL.'.format(round(volumef*.25, 1), volume))
            k = input('')
            print('Place all PCR tubes in thermal cycler.')
            l = input('')
            print('A typical cycle contains 3 steps:\n1. denaturation: 30 seconds at 95 C\n2. annealing: 30 seconds at 50-65* C \n3. extension: 30-60** seconds at 72*** C\n\n*All primers have different optimal annealing temperatures \n**Longer amplicons need longer extension times\n***Temperature is polymerase-dependent')
            m = input('')
            print('These three steps will be repeated 25-35 times, depending on amount of amplicon desired.\nIn addition, an ititial denaturation of 4 minutes & a final extension of 5 minutes are included.')
            n = input('')
            print('Reaction will take between 2 and 4 hours to complete, depending on extension time and number of cycles chosen.')
            o = input('')
            

    def gelProtocol(self):
        '''
        The gelProtocol method is responsible for printing dynamic instructions for how to load and run a gel to estimate DNA length. The numbers list contains integers 0-9 and is used to clean the
        input so that it only contains integers. After prompting the user, steps are printed out such that the user must press enter to proceed to the next step. After the third print statement, the user
        is prompted to enter the % agarose gel desired. This value is stored as p, and is cleaned so that it only contains integers. The rest of the protocol is then printed.

        Pseudo Code:
        numbers = [0-9]
        ask user if they want protocol for running gel electrophoresis
        if user returns yes:
            print introduction to gel electrophoresis
            p = input(ask user for percent agarose gel desired)
            join p by whitespace, then split into individual characters
            join p by whitespace and capatilize
            for character in cleaned up p:
                if character is not a number 0 - 9:
                    replace the character with whitespace
                    make p an integer
            set d equal to float of p/10
            print protocol for running gel step by step so that user must press enter to proceed to next step
        '''
        numbers = ['0','1','2','4','5','6','7','8','9']
        prompt = input('Print protocol for gel electrophoresis? (yes/no)\n')
        if prompt == 'yes':
            print('Gel electrophoresis is an useful tool for estimating DNA molecular weight/length.')
            print('This procedure demonstrates how to run an agarose gel. Press enter to procede to the next step.')
            a = input()
            print('Ingredients needed: \n    1.DNA\n    2.Molecular weight marker (ladder) \n    3.10X DNA stain (makes DNA visible under UV light)\n    4.6X loading dye')
            p = input('What percent agarose is desired? Note: most gels are between .7% and 2.0%\n')
            rawSeq = ''.join(p).split()
            cleanP= ''.join(rawSeq).upper()
            for char in cleanP:
                if char not in numbers:
                    p = cleanP.replace(char,'')
                    p = int(p)
            d = float(p/10)
            print('This guide is for a {0}% gel.\nAdd {0} g aragose powder to 100 mL deionized water. Mix and microwave until all agarose is dissolved (1-2 minutes).'.format(round(d,2)))
            a = input()
            print('While waiting for agarose to cool, set up casting tray containing a gel tray with well combs and prepare samples/ladder for loading.')
            a = input()
            print('For 4.4 uL of each sample DNA in a snap-cap tube, add:\n    1 uL of 6X loading dye\n    .6 uL DNA stain\nMix by flicking tube.')
            a = input()
            print('For 5 uL of ladder in a snap-cap tube, add:\n    1.35 uL of 6X loading dye\n    .5 uL DNA stain\nMix by flicking tube.')
            print('\nKeep samples and molecular weight marker on ice until ready to load gel.')
            a = input()
            print('Once agarose is cool but still liquid, pour into gel tray.')
            a = input()
            print('Fill casting tray with 1X running buffer until gel is covered.\nGently remove well combs.')
            a = input()
            print('Load 2.5 uL of DNA ladder into the first well to be used as a reference.')
            a = input()
            print('Load 5 uL of each sample mixture into each following well.\nBe sure to document what each well contains.')
            a = input()
            print('Once all samples have been loaded, run the gel at 80-150V.\nIMPORTANT: DNA is negativily charged, always run gel towards red (positive) wire!')
            a = input()
            print('Once dye line is approximately 75-80% of the way down the gel, turn off power source.\nVisualize your DNA using any device with UV light capabilities.\n')


    def REfinder(self):
        """
        The REfinder method is responsible for finder restriction enzyme cutsites within the amplicon. After scanning the amplicon sequence, a new sequence is created that shows
        the locations of the cutsites and the restriction enzyme that recognizes the cutsite. A count of all enzymes found in sequence is also listed. If no cutsites are found
        within the amplicon, an exception is raised and a message is printed to the user notifying them that no cutsites were found. The user is then prompted if they want the
        fragment lengths created by treating the amplicon with specific enzyme(s). If the inputted enzyme does not have a cutsite within the amplicon, an exception is raised
        stating that the inputted enzyme does not create any fragments of the amplicon.

        Pseudo code:
        set amp equal to pcr product
        ititilize enzyme as an empty dictionary to store specific enzymes found in amplicon
        itilitze enzymeList as an empty list
        for enzyme name, cut sequence in self.REdict.items():
            set cutseq equal to cutsite without | (denotes cleavage site)

            if cutSeq is found in amp:
                add name:cutsite to enzyme
                add name to enzymeList
                start position of product = finding position of forward primer in sequence
                end position of product = finding position of reverse primer in sequence + lenght of reverse primer
                cut sequence = start to stop position in amp
                sequence with cutsites = replace cutsequence with __cutsite__
                set amp to sequence with cutsites so that each enzyme gets added to sequence

        ask user if they want cutsites shown in amplified sequence
        if yes:
            try:
                print sequence with cutsites
                print list of enzymes with cutsites in sequence
                for key in self.REdict.keys():
                   if key is in sequence with cutsites:
                       count each enzyme
                       print enzyme name and count
            except UnbounDLocalError:
                print that no cutsites exist in sequence

            prompt user to show specific cutsites in amplicon
            if yes:
                prompt user for desired restriction enzymes
                add input to list
                make a dictionary copying key:value from REdict with specified restriction enzymes
                make copy of pcr product
                make another copy of pcr product
                for name, cutsite in newly created dictionary:
                    remove | from cutsite
                    if cutsite is in amp:
                        find start position
                        find stop postion
                        specificRE = replace start:stop from sequence with cutsite
                        noName = replace start:stop from sequence with cutsite(no enzyme name)
                        set amp to specificRE
                        set amp2 to noName
                        make list by splitting sequence on |
                    print expected fragment lengths
                    try:
                        loop through fragments
                        print length of each fragment
                    except:
                        print that no fragments are made with specified restriction enzyme
        """
        amp = self.pcrProd()
        enzyme = {}
        thisEnzyme = {}
        enzymeList = []
        for key, cutsite in self.REdict.items():
            cutSeq = cutsite.replace('|', '')
            
            if cutSeq in amp:
                enzyme[key] = cutsite
                enzymeList.append(key)
                start = amp.find(cutSeq)
                end = start + len(cutSeq)
                cs = amp[start:end]
                insertRE = amp.replace(cs, '__'+ cutsite +'__('+key+')')
                amp = insertRE

        prompt = input('Find restriction enzyme cutsites in amplicon? (yes/no) \n')
        if prompt == 'yes':
            try:
                print('5 ' + insertRE +' 3\n')
                print('Cutsites present in sequence: ')
                for key in self.REdict.keys():
                    if key in insertRE:
                        REcount = insertRE.count(key)
                        print('{0}: {1}'.format(key, REcount))
            except UnboundLocalError:
                print('No cutsites found in sequence.\n')
            thisPrompt = input('Digest amplicon with specific primer(s)? (yes/no)\n')
            if thisPrompt == 'yes':
                RE = input('List desired restriction enzyme names separated by a space \n')
                RElist = RE.split(' ')
                thisEnzyme = {key:value for key, value in self.REdict.items() if key in RElist}
                amp = self.pcrProd()
                amp2 = amp
                for key, cutsite in thisEnzyme.items():
                    cutSeq = cutsite.replace('|', '')
                    if cutSeq in amp:
                        start = amp.find(cutSeq)
                        end = start + len(cutSeq)
                        cs = amp[start:end]
                        specificRE = amp.replace(cs, '__'+ cutsite +'('+key+')__')
                        noName = amp.replace(cs, cutsite)
                        amp = specificRE
                        amp2 = noName
                        fragments = noName.split('|')
                print('Expected fragment lengths:')
                try:
                    for frag in fragments:
                        fragLen = len(frag)
                        print ('{0} base pairs'.format(fragLen))
                except:
                    print('\nNo fragments made with this restriction enzyme.')


        else:
            pass
        

    def pcrProd(self):
        """
        The pcrProd method sets forward and reverse to the forward and reverse primers, respectively. The start position of the PCR product is determined by finding the
        forward primer in the template sequence. The end position is determined by finding the start position of the reverse primer, and adding the length of the reverse
        primer to this value. The amplicon is saved as a fraction of the template sequence from start value to end value.

        Pseudo code:
        set forward equal to the value of f in primerdict
        set reverse equal to the value or r in primerdict
        find start in sequence
        find stop in sequence
        set product equal to start:stop in sequence
        return product
        """
        forward = self.primerDict['f']
        reverse = self.primerDict['r']
        start = self.seq.find(forward)
        end = self.seq.find(reverse) + len(reverse)
        product = self.seq[start:end]
        return product

    
    def getInfo(self):
        """
        The getInfo method is responsible for getting the two primer sequences from the user, and saving the sequences in self.primerDict to be referenced in other methods.
        The GC content and molecular weight is found and stored for each primer. The user is then prompted to answer if the primers are both located on the same strand of DNA.
        If the user inputs no, then the reverse complement of the reverse primer is found so that the program can find the primer in the single-stranded template sequence.
        The GC content and molecular weight for both primers are then printed. If the GC content for a primer is below 50%, a warning is attached advising the user to select
        another primer due to specificity issues.

        Pseudo code:
        ask for forward primer
        add forward primer to dictionary
        ask for reverse primer
        add reverse primer to dictionary

        for primer in primerdict.values:
            make primer a string
            count A, G, T, C
            calculate GC content
            calcualte molecular weight

        prompt user if primers are from same strand of DNA
        if no:
            find reverse complement of primer
            replace this sequence in primer dictionary

        for primer in primerdict.values:
            make primer a string
            count A, G, T, C
            calculate GC content
            calcualte molecular weight
            print forward primer with statistics
            print reverse primer with statistics
            print warning if GC content < 50%

        """
        forward = input('Please enter the forward primer sequence: ')
        self.primerDict['f'] = forward.upper()
        reverse = input('Please enter the reverse primer sequence: ')
        self.primerDict['r'] = reverse.upper()
        
        answer = input('Reverse compliment of reverse primer sequence? (yes/no) \n')
        if answer == 'no':
            temp = list(reverse.upper())
            rb = [self.compDict[nt] for nt in reversed(temp)]
            reverse = ''.join(rb)
            self.primerDict['r'] = reverse
            
        for primer in self.primerDict.values():
            primert = str(primer)
            ac = primert.count('A')
            gc = primert.count('G')
            tc = primert.count('T')
            cc = primert.count('C')
            gcContent = ((gc+cc)/(gc+cc+ac+tc))*100
            ampWeight = (ac*651.43)+(gc*652.41)+(cc*652.41)+(tc*651.43)
            ampWeight = ampWeight/1000
            print('')
            if self.primerDict['f'] is primer:
                print('Forward primer')
            if self.primerDict['r'] is primer:
                print('Reverse primer')
            print('Primer length: {0} nucleotides'.format(len(primert)))
            print('Molecular weight: {0} kDa'.format(round(ampWeight, 2)))
            print('GC content: {0}%'.format(round(gcContent,2)))
            if gcContent < 50:
                print('WARNING: low GC content. Consider choosing a different primer.')

    def main():
        """
        The main method is responsible for getting the template DNA sequence, either from a file or from input. Once the template DNA is obtained, an instance of
        the class PCR is made (bp) using this input. A boolean function is included and set to false in case the user wants to rerun the program after getting to the end of the program.
        Methods of this instance are then individually called on. Once all methods have been cycled though, the user is prompted if they would like to rerun the program
        with a differnt pair of primers. If the user inputs yes, then the boolean is returned false and the program starts over.

        Take input from standard in or from a file
        if input from standard in:
            if input > 0:
                make instance of PCR class with sequence
                make boolean set to false
                while no boolean:
                    run through methods of PCR class
                    if amplicon length is 0:
                        continue
                    continue running though methods
                    set p equal to input(run again?)
                        if yes:
                            boolean = false
                        else:
                            boolean = true
            else:
                end program
        else:
            input from file
            loop through head and sequence of input file
            if input > 0:
                make instance of PCR class with sequence
                make boolean set to false
                while no boolean:
                    run through methods of PCR class
                    if amplicon length is 0:
                        continue
                    continue running though methods
                    set p equal to input(run again?)
                        if yes:
                            boolean = false
                        else:
                            boolean = true
        """
        refFile = input('Please enter the name of the file containing the template sequence: \n(type "no file" and/or press enter to input your own sequence)\n')
        if refFile == 'no file' or refFile == '':
            terminalSeq = input('')
            if len(terminalSeq) > 0:
                bp = PCR(terminalSeq)
                booleani = False
                while not booleani:
                    bp.getInfo()
                    bp.pcrProd()
                    amplicon = bp.amplicon()
                    if amplicon == -1:
                        continue
                    bp.REfinder()
                    bp.protocol()
                    bp.gelProtocol()
                    p = input('Would you like to run another PCR reaction with a different primer set? (yes/no)\n')
                    if p == 'yes':
                        booleani = False
                    else:
                        booleani = True
            else:
                print('No sequence entered, program terminating.')
        else:
            myReader = FastAreader(refFile)#refFile
            for head, seq in myReader.readFasta():
                thisSeq = seq
                bp = PCR(thisSeq)
                booleani = False
                while not booleani:
                    bp.getInfo()
                    bp.pcrProd()
                    amplicon = bp.amplicon()
                    if amplicon == -1:
                        continue
                    bp.REfinder()
                    bp.protocol()
                    bp.gelProtocol()
                    p = input('Would you like to run another PCR reaction with a different primer set? (yes/no)\n')
                    if p == 'yes':
                        booleani = False
                    else:
                        booleani = True

PCR.main()
            
        

        
    
