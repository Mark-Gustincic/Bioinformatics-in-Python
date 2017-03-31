#!/usr/bin/env python3
# Name: Mark Gustincic
# Group Members: Thomas Richards(tarichar)

class FastQParse (str):
    def __init__(self,seq):
        self.mySeqNew = seq
    def parse(self):
        seq = self[1:]
        seq2 = seq.split(":")
        return(seq2)
    """
    Prompts user for input, and assigns input to name 'seqNameLine'.
    seqNameLine is run through parse method, and then assigned the
    name pn. Names are assigned to objects 0-6 of seqNew.
    """
    def main():
        seqNameLine = input("Please input the seqname line of a FASTQ file: ")
        seqNew = FastQParse(seqNameLine)
        pn = seqNew.parse()
        Instrument = pn[0]
        Run_ID = pn[1]
        Flow_Cell_ID = pn[2]
        Flow_Cell_Lane = pn[3]
        Tile_Number = pn[4]
        x_coord = pn[5]
        y_coord = pn[6]
    
        print('Instrument = {0}'.format(Instrument))
        print('Run ID = {0}'.format(Run_ID))
        print('Flow Cell ID = {0}'.format(Flow_Cell_ID))
        print('Flow Cell Lane = {0}'.format(Flow_Cell_Lane))
        print('Tile Number = {0}'.format(Tile_Number))
        print('x-coord = {0}'.format(x_coord))
        print('y-coord = {0}'.format(y_coord))
        
        """
        Formats the seven parts of the parsed input and prints them
        to output on separate lines. Method main is called to invoke
        function.
        """
        

FastQParse.main()
        


