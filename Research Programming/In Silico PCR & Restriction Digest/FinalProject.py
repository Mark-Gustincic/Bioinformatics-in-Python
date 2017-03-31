#!usr/bin/env python3
# Name: Mark Gustincic (mgustinc)
# Group Members: Thomas Richards (tarichar)

import sequenceAnalysis

class PCR:

    def __init__(self, seq):

        bases = ['A', 'T', 'C', 'G']
        self.compDict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
        self.primerDict = {'f':0, 'r':0}
        self.dnaMass={'G': 652.41, 'C': 652.41, 'A': 651.43, 'T': 651.43}
        self.guide = {12:'Place all PCR tubes in thermal cycler.',
                      13:'A typical cycle contains 3 steps:\n1. denaturation: 30 seconds at 95 C\n2. annealing: 30 seconds at 50-65* C \n3. extension: 30-60** seconds at 72*** C\n\n*All primers have different optimal annealing temperatures \n**Longer amplicons need longer extension times\n***Temperature is polymerase-dependent',
                      14:'These three steps will be repeated 25-35 times, depending on amount of amplicon desired.\nIn addition, an ititial denaturation of 4 minutes & a final extension of 5 minutes are included.',
                      15:'Reaction will take between 2 and 4 hours to complete, depending on extension time and number of cycles chosen.'
                      }
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
        forward = self.primerDict['f']
        reverse = self.primerDict['r']
        product = self.pcrProd()
        if len(product) == 0:
            print('There are no matches to {0} and {1} in the template sequence.'.format(forward, reverse))
            self.getInfo()
            self.pcrProd()
        print('')
        prompt = input('Show amplicon sequence? (yes/no) ')
        print('')
        if prompt == 'yes':
            print('Expected product sequence: \n5 {0} 3\n'.format(product))
        for x in product:
            ac = product.count('A')
            gc = product.count('G')
            tc = product.count('T')
            cc = product.count('C')
        gcContent = ((gc+cc)/(gc+cc+ac+tc))*100
        ampWeight = (ac*651.43)+(gc*652.41)+(cc*652.41)+(tc*651.43)
        ampWeight = ampWeight/1000
        print('Amplicon length: {0} base pairs'.format(len(product)))
        print('Molecular weight: {0} kDa'.format(round(ampWeight,2)))
        print('GC content: {0}%\n'.format(round(gcContent,2)))
        return ampWeight
            

    def protocol(self):
        x = 12
        prompt2 = input('\nPrint PCR protocol? (yes/no) ')
        if prompt2 == 'yes':
            volume = input('What is the desired reaction volume in uL?\n')
            while True:
                try:
                    volumef = float(volume)
                    break
                except ValueError:
                    volume = input('Please enter a number\n')
                    continue
            if volume == '':
                volume = input('Please input a number\n')
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
            print('Add {0} uL of template DNA to each PCR tube, bringing total reaction volume to {1} uL.\n'.format(round(volumef*.25, 1), volume))
            for step in self.guide.keys():
                if step == x:
                    print(self.guide[step])
                    x += 1
                    z = input()
            return volume

    def gelProtocol(self):
        numbers = [0,1,2,4,5,6,7,8,9]
        prompt = input('Print protocol for gel electrophoresis? (yes/no)\n')
        if prompt == 'yes':
            print('Gel electrophoresis is an useful tool for estimating DNA molecular weight/length.')
            print('This procedure demonstrates how to run an agarose gel. Press enter to procede to the next step.')
            a = input()
            print('Ingredients needed: \n    1.DNA\n    2.Molecular weight marker (ladder) \n    3.10X SyBr green (makes DNA visible under UV light)\n    4.6X loading dye')
            a = input()
            p = input('What percent agarose is desired? Note: most gels are between .7% and 2.0%\n')
            rawSeq = ''.join(p).split()
            cleanP= ''.join(rawSeq).upper()

            for char in cleanP:
                if char not in numbers:
                    p = cleanP.replace(char,'')
            p = float(p)
            print('This guide is for a {0}% gel.\nAdd {0} g aragose powder to 100 mL deionized water. Mix and microwave until all agarose is dissolved (1-2 minutes).'.format(p))
            a = input()
            print('While waiting for agarose to cool, set up casting tray containing a gel tray with well combs and prepare samples/ladder for loading.')
            a = input()
            print('For 4.4 uL of sample DNA in a snap-cap tube add:\n    1 uL of 6X loading dye\n    .6 uL SyBr green\nMix by flicking tube.')
            a = input()
            print('For 5 uL of ladder in a snap-cap tube Add:\n    1.35 uL of 6X loading dye\n    .5 uL Sybr green\nMix by flicking tube.')
            print('Keep samples and molecular weight marker on ice until ready to load gel.')
            a = input()
            print('Once agarose is cool but still liquid, pour into gel tray.')
            a = input()
            print('Fill casting tray with 1X TAE or TBE (buffer) until gel is covered.\nGently remove well combs.')
            a = input()
            print('Load 2.5 uL of DNA ladder into the first well to be used as a reference.')
            a = input()
            print('Load 5 uL of sample mixture into each following well.\nBe sure to document what each well contains.')
            a = input()
            print('Once all samples have been loaded, run the gel at 80-150V.\nIMPORTANT: DNA is negativily charged, always run gel towards red(positive) wire!')
            a = input()
            print('Once dye line is approximately 75-80% of the way down the gel, turn off power source.\nVisualize your DNA using any device with UV light capabilities.')

                       
    def REfinder(self):
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
                insertRE = amp.replace(cs, '__'+ cutsite +'('+key+')__')
                amp = insertRE
            
        x = 1
        prompt = input('Show restriction enzyme cutsites in amplicon? (yes/no) \n')
        if prompt == 'yes':
            if insertRE:
                print('5 ' + insertRE +' 3\n')
            else:
                print('No cutsites present in sequence')

        else:
            pass
        print('Cutsites present in sequence: ')
        for key in self.REdict.keys():
            if key in insertRE:
                REcount = insertRE.count(key)
                print('{0}: {1}'.format(key, REcount))

        thisPrompt = input('Show specific cutsites within amplicon? ')
        if thisPrompt == 'yes':
            RE = input('List desired restriction enzyme names separated by a space \n')
            RElist = RE.split(' ')
            thisEnzyme = {key:value for key, value in self.REdict.items() if key in RElist}
            amp = self.pcrProd()    
            for key, cutsite in thisEnzyme.items():
                cutSeq = cutsite.replace('|', '')
                if cutSeq in amp:
                    start = amp.find(cutSeq)
                    end = start + len(cutSeq)
                    cs = amp[start:end]
                    specificRE = amp.replace(cs, '__'+ cutsite +'('+key+')__')
                    amp = specificRE
                    #print(specificRE)
                    fragments = specificRE.split('|')
                    #fragList = [fragments]
            print('Expected fragment lengths:')
            #for fragments in fragList:
            for frag in fragments:
                fragLen = len(frag)
                print ('{0} base pairs'.format(fragLen))
            if not insertRE:
                print('No cutsites present in amplicon')
        return insertRE
        

    def pcrProd(self):
        forward = self.primerDict['f']
        #print (forward)
        reverse = self.primerDict['r']
        #print (reverse)
        start = self.seq.find(forward)
        end = self.seq.find(reverse) + len(reverse)
        product = self.seq[start:end]
        return product

    
    def getInfo(self):
        forward = input('Please enter the forward primer sequence: ')
        #forward = 'CCCTATCTCATCTGTCTCCCTTA'
        self.primerDict['f'] = forward.upper()
        reverse = input('Please enter the reverse primer sequence: ')
        #reverse = 'TTGGTACTCGTCAACGAGGATATATTT'
        self.primerDict['r'] = reverse.upper()

        for primer in self.primerDict.values():
            primer = str(primer)
            ac = primer.count('A')
            gc = primer.count('G')
            tc = primer.count('T')
            cc = primer.count('C')
            gcContent = ((gc+cc)/(gc+cc+ac+tc))*100
            ampWeight = (ac*651.43)+(gc*652.41)+(cc*652.41)+(tc*651.43)
        
        answer = input('Are primers on the same strand? (yes/no) ')
        if answer == 'no':
            temp = list(reverse.upper())
            rb = [self.compDict[nt] for nt in reversed(temp)]
            reverse = ''.join(rb)
            self.primerDict['r'] = reverse
            
        for primer in self.primerDict.values():
            primer = str(primer)
            ac = primer.count('A')
            gc = primer.count('G')
            tc = primer.count('T')
            cc = primer.count('C')
            gcContent = ((gc+cc)/(gc+cc+ac+tc))*100
            ampWeight = (ac*651.43)+(gc*652.41)+(cc*652.41)+(tc*651.43)
            ampWeight = ampWeight/1000
            print('')
            print(primer)
            print('Primer length: {0} nucleotides'.format(len(primer)))
            print('Molecular weight: {0} kDa'.format(round(ampWeight, 2)))
            print('GC content: {0}%'.format(round(gcContent,2)))
            if gcContent < 50:
                print('WARNING: low GC content. Consider choosing a different primer.')

    def main():
        refFile = input('Please enter the name of the file containing the template sequence: ')
        myReader = sequenceAnalysis.FastAreader(refFile)#refFile
        for head, seq in myReader.readFasta():
            thisSeq = seq
            bp = PCR(thisSeq)
            bp.getInfo()
            bp.pcrProd()
            bp.amplicon()
            bp.REfinder()
            bp.protocol()
            bp.gelProtocol()
            break

PCR.main()
            
        

        
    
