from . import algorithms as alg
from .pre_algorithms import k_mer_indexing


class InputChecker():

    @staticmethod
    def dnacheck(seqs):
        resultdict = {}

        nucleotides = 'ACGT'
        check = False

        for char in seqs[0]:
            if nucleotides.find(char) == -1:
                check = True
                break

        if not check:
            for char in seqs[1]:
                if nucleotides.find(char) == -1:
                    check = True

        return check


    @staticmethod
    def patternlengthchecker(seqs):

        if len(seqs[0]) > len(seqs[1]):
            return True

        return False
    

    @staticmethod
    def boyerMooreApproximatematchingpatternlengthchecker(seqs):

        if len(seqs[0]) <= 5 or len(seqs[0]) > len(seqs[1]):
            return True

        return False


    @staticmethod
    def hammingdistanceseqlenghchecker (seqs):

        if len(seqs[0]) != len(seqs[1]):
            return True

        return False
         


    @staticmethod
    def indexAndDistanceInputTypeChecker(textinputs):

        try:
            return not textinputs['index'].isdigit()
        except KeyError:
            try:
                return not textinputs['max distance'].isdigit()
            except KeyError:
                pass       

        return False
        
    

    @staticmethod
    def localAlignmentScoreTypeChecker(textinputs):

        if  textinputs['match'].isdigit() and  textinputs['mismatch'].isdigit() and  textinputs['gap'].isdigit():
            return True

        return False


    @staticmethod
    def globalAlignmentScoreTypeChecker(textinputs):

        if  textinputs['match'].isdigit() and  textinputs['transition'].isdigit() and textinputs['transversion'].isdigit() and textinputs['gap'].isdigit():
            return False

        return True



    @staticmethod
    def indexAndDistanceSignChecker(textinputs, algotype):

        if algotype == 1:
            if int(textinputs['index']) > len(textinputs['seq1']) or int(textinputs['index']) < 1:
                return True

        elif algotype == 2:
            if int(textinputs['max distance']) == 0 :
                return True

        return False



    @staticmethod
    def localAlignmentScoresSignChecker(textinputs):

        if int(textinputs['match']) <= 0 or int(textinputs['mismatch']) >= 0 or int(textinputs['gap']) >= int(textinputs['mismatch']):
            return True

        return False 


    @staticmethod
    def globalAlignmentScoresSignChecker(textinputs):

        if int(textinputs['match']) < 0 or int(textinputs['transition']) <= int(textinputs['match']) or int(textinputs['transversion']) <= int(textinputs['transition']) or int(textinputs['gap']) <= int(textinputs['transversion']) :
            return True

        return False 

