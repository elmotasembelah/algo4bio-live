from django.utils.datastructures import MultiValueDictKeyError
from . import Filereader as fr


class InputHandler:

    @staticmethod
    def fileinputhandler(request):
        fileinputs = []

        # getting the first uploaded file from the form 
        try:
            fileinput = request.FILES['document1']
        except MultiValueDictKeyError:
            fileinput = ''
            fileinputs.append(fileinput)
        else:
             
            fileinputs.append(fr.ReadingFiles.fastafilereader(request, 'document1').upper())
 
        # second uploaded file
        try:
            fileinput = request.FILES['document2']
        except MultiValueDictKeyError:
            fileinput = ''
            fileinputs.append(fileinput)
        else:
            fileinputs.append(fr.ReadingFiles.fastafilereader(request, 'document2').upper())

        return fileinputs

    @staticmethod
    def seqEmptyInputHandler(textinputs, fileinputs):


        if (textinputs['seq1'] == '' and fileinputs[0] == '') or (textinputs['seq2'] == '' and fileinputs[1] == ''):
            return True
            
        return False


    @staticmethod
    def multipleinputhandler(textinputs, fileinputs):
        if (textinputs['seq1'] != '' and fileinputs[0] != '') or (textinputs['seq2'] != '' and fileinputs[1] != ''):
            return True
            
        return False

    
    @staticmethod
    def indexAndMaxDistanceEmptyInputHandler(textinputs):

        try:
            if textinputs['index'] == '':
                return True 
            
        except KeyError:
            try:
                if textinputs['max distance'] == '':
                    return True 
            except KeyError:
                pass        

        return False

    @staticmethod
    def localAlignmentScoreEmptyInputHandler(textinputs):

        if textinputs['match'] == '' or textinputs['mismatch'] == '' or textinputs['gap'] == '':
            return True

        return False

    
    @staticmethod
    def globalAlignmentScoreEmptyInputHandler(textinputs):

        if textinputs['match'] == '' or textinputs['transition'] == '' or textinputs['transversion'] == '' or textinputs['gap'] == '':
            return True

        return False
        

    


    @staticmethod
    def inputtypehandler(textinputs, fileinputs):
        seqs = []

        if textinputs['seq1'] != '':
            seqs.append(textinputs['seq1'])
        else:
            seqs.append(fileinputs[0])

        
        if textinputs['seq2'] != '':
            seqs.append(textinputs['seq2'])   
        else:
            seqs.append(fileinputs[1])


        return seqs