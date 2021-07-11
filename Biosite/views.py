from django.http import QueryDict
from django.shortcuts import render, redirect
import numpy as np
from .classes import InputChecker as ic 
from .classes import Inputhandler as ih 
from .classes import algorithms as alg
from .classes.pre_algorithms import k_mer_indexing


def homepage(request):

        return render(request, template_name='main/homepage.html')
 


def globalAlignmentPage(request):

    if 'back' in request.GET or 'back' in request.POST:
        return redirect ('/home/')

    if request.method == 'GET':
        if'sister' in request.GET:
            return redirect('/localAlignment/')
        else:
            return render(request, template_name='main/algorithms pages/globalAlignment.html')
 
    elif request.method == 'POST':
        result = {}
          

        textinputs = QueryDict.copy(request.POST)
        fileinputs = ih.InputHandler.fileinputhandler(request) 

        textinputs['seq1'] = str(textinputs['seq1'].upper())
        textinputs['seq2'] = str(textinputs['seq2'].upper())


        result['error'] = checkEverythingBeforeAlgorithm(textinputs, fileinputs, 5) 
    
        if result['error']: 

            return render(request, template_name='main/Error pages/globalAlignmentError.html', context={'result': result})
        
        costs = [int(textinputs['match']), int(textinputs['transition']), int(textinputs['transversion']), int(textinputs['gap'])]
        result = alg.Algorithms.globalAlignment(ih.InputHandler.inputtypehandler(textinputs, fileinputs), costs)

        return render(request, template_name='main/Result pages/globalAlignmentResult.html', context={'result': result})


def localAlignmentPage(request):
    
    if 'back' in request.GET or 'back' in request.POST:
        return redirect ('/home/')

    if request.method == 'GET':
        if'sister' in request.GET:
            return redirect('/globalAlignment/')
        else:
            return render(request, template_name='main/algorithms pages/localAlignment.html')
 
    elif request.method == 'POST':
        result = {}
          

        textinputs = QueryDict.copy(request.POST)
        fileinputs = ih.InputHandler.fileinputhandler(request) 

        textinputs['seq1'] = str(textinputs['seq1'].upper())
        textinputs['seq2'] = str(textinputs['seq2'].upper())


        result['error'] = checkEverythingBeforeAlgorithm(textinputs, fileinputs, 4) 

        if result['error']: 

            return render(request, template_name='main/Error pages/localAlignmentError.html', context={'result': result})

        costs = [int(textinputs['match']), int(textinputs['mismatch']), int(textinputs['gap'])]
        result = alg.Algorithms.localAlignment(ih.InputHandler.inputtypehandler(textinputs, fileinputs), costs)

        return render(request, template_name='main/Result pages/localAlignmentResult.html', context={'result': result})


def editDistancePage(request):

    if 'back' in request.GET or 'back' in request.POST:
        return redirect ('/home/')

    if request.method == 'GET':
        if'sister' in request.GET:
            return redirect('/hammingDistance/')
        else:
            return render(request, template_name='main/algorithms pages/editDistance.html')
 
    elif request.method == 'POST':
        result = {}
          

        textinputs = QueryDict.copy(request.POST)
        fileinputs = ih.InputHandler.fileinputhandler(request) 

        textinputs['seq1'] = str(textinputs['seq1'].upper())
        textinputs['seq2'] = str(textinputs['seq2'].upper())

        result['error'] = checkEverythingBeforeAlgorithm(textinputs, fileinputs, 'edit distance') 
    
        if result['error']: 

            return render(request, template_name='main/Error pages/editDistanceError.html', context={'result': result})
        
        result = alg.Algorithms.editDistance(ih.InputHandler.inputtypehandler(textinputs, fileinputs))

        return render(request, template_name='main/Result pages/editDistanceResult.html', context={'result': result})


def hammingDistancePage(request):

    if 'back' in request.GET or 'back' in request.POST:
        return redirect ('/home/')

    if request.method == 'GET':
        if'sister' in request.GET:
            return redirect('/editDistance/')
        else:
            return render(request, template_name='main/algorithms pages/hammingDIstance.html')
 
    elif request.method == 'POST':
        result = {}
          

        textinputs = QueryDict.copy(request.POST)
        fileinputs = ih.InputHandler.fileinputhandler(request) 

        textinputs['seq1'] = str(textinputs['seq1'].upper())
        textinputs['seq2'] = str(textinputs['seq2'].upper())

        result['error'] = checkEverythingBeforeAlgorithm(textinputs, fileinputs, 6) 
    
        if result['error']: 

            return render(request, template_name='main/Error pages/hammingDistanceError.html', context={'result': result})
        
        result = alg.Algorithms.hammingDistance(ih.InputHandler.inputtypehandler(textinputs, fileinputs))

        return render(request, template_name='main/Result pages/hammingDistanceResult.html', context={'result': result})


def hashmapKmerIndexingPage(request):

    if 'back' in request.GET or 'back' in request.POST:
        return redirect ('/home/')

    if request.method == 'GET':
            
        if'sister' in request.GET:
            return redirect('/binarySearchKmerIndexing/')
        else:
            return render(request, template_name='main/algorithms pages/hashMapKmerIndexing.html')
 
    elif request.method == 'POST':
        result = {}    

        textinputs = QueryDict.copy(request.POST)
        fileinputs = ih.InputHandler.fileinputhandler(request) 

        textinputs['seq1'] = str(textinputs['seq1'].upper())
        textinputs['seq2'] = str(textinputs['seq2'].upper())
    


        result['error'] = checkEverythingBeforeAlgorithm(textinputs, fileinputs, 2) 
    
        if result['error']: 

            return render(request, template_name='main/Error pages/hashmapKmerIndexError.html', context={'result': result})

        if textinputs['seq2'] == '':
            template = fileinputs[1]
        else:
            template = textinputs['seq2']

        index = k_mer_indexing.Hashmap_Index(template, int(textinputs['index']))
        result = alg.Algorithms.hashmapIndex(ih.InputHandler.inputtypehandler(textinputs, fileinputs), index)

        return render(request, template_name='main/Result pages/hashmapKmerIndexResult.html', context={'result': result})


def binarySearchKmerIndexingPage(request, temp = ''):

    if 'back' in request.GET or 'back' in request.POST:
        return redirect ('/home/')

    if request.method == 'GET':   
        if'sister' in request.GET:
            return redirect('/hashmapKmerIndexing/')
        else:
            return render(request, template_name='main/algorithms pages/binarySearchKmerIndexing.html')
    
    if request.method == 'POST':
        

        result = {} 

        textinputs = QueryDict.copy(request.POST)
        fileinputs = ih.InputHandler.fileinputhandler(request) 

        textinputs['seq1'] = str(textinputs['seq1'].upper())
        textinputs['seq2'] = str(textinputs['seq2'].upper())

        
    
        result['error'] = checkEverythingBeforeAlgorithm(textinputs, fileinputs, 2) 

        if result['error']: 

            return render(request, template_name='main/Error pages/binarySearchKmerIndexError.html', context={'result': result})

        if textinputs['seq2'] == '':
            template = fileinputs[1]
        else:
            template = textinputs['seq2']
            

        index = k_mer_indexing.Binarysearch_Index(template, int(textinputs['index']))

        result = alg.Algorithms.binarySearchKmerIndex(ih.InputHandler.inputtypehandler(textinputs, fileinputs), index)
        
        return render(request, template_name='main/Result pages/binarySearchKmerIndexResult.html', context={'result': result})


def boyerMooreApproximatematchingPage(request):

    if 'back' in request.GET or 'back' in request.POST:
        return redirect ('/home/')

    if request.method == 'GET':
        if'sister' in request.GET:
            return redirect('/naiveApproximateMatching/')
        else:
            return render(request, template_name='main/algorithms pages/boyerMooreApproximateMatching.html')
 
    elif request.method == 'POST':
        result = {}
          
        textinputs = QueryDict.copy(request.POST)
        fileinputs = ih.InputHandler.fileinputhandler(request) 

        textinputs['seq1'] = str(textinputs['seq1'].upper())
        textinputs['seq2'] = str(textinputs['seq2'].upper())

        
        

        result['error'] = checkEverythingBeforeAlgorithm(textinputs, fileinputs, 'boyer moore approximate matching') 
    
        if result['error']: 

            return render(request, template_name='main/Error pages/boyerMooreApproximateMatchingError.html', context={'result': result})

        self = alg.Algorithms
        result = alg.Algorithms.boyerMooreApproximatematching(self, ih.InputHandler.inputtypehandler( textinputs, fileinputs), int(textinputs['max distance']))

        return render(request, template_name='main/Result pages/boyerMooreApproximateMatchingResult.html', context={'result': result})


def naiveApproximateMatchingPage(request):

    if 'back' in request.GET or 'back' in request.POST:
        return redirect ('/home/')

    if request.method == 'GET':

        if'sister' in request.GET:
            return redirect('/boyerMooreApproximatematching/')
        else:
            return render(request, template_name='main/algorithms pages/naiveApproximateMatching.html')
 
    elif request.method == 'POST':
        result = {}
          

        textinputs = QueryDict.copy(request.POST)
        fileinputs = ih.InputHandler.fileinputhandler(request) 

        textinputs['seq1'] = str(textinputs['seq1'].upper())
        textinputs['seq2'] = str(textinputs['seq2'].upper())
    
        

        result['error'] = checkEverythingBeforeAlgorithm(textinputs, fileinputs, 3) 
    
        if result['error']: 

            return render(request, template_name='main/Error pages/naiveApproximateMatchingError.html', context={'result': result})

        result = alg.Algorithms.naiveApproximateMatching(ih.InputHandler.inputtypehandler(textinputs, fileinputs), int(textinputs['max distance']))

        return render(request, template_name='main/Result pages/naiveApproximateMatchingResult.html', context={'result': result})


def boyer_MooreExactMatchingPage(request):

    if 'back' in request.GET or 'back' in request.POST:
        return redirect ('/home/')

    if request.method == 'GET':
            
        if'sister' in request.GET:
            return redirect('/naiveExactMatching/')
        else:
            return render(request, template_name='main/algorithms pages/boyerMooreExactMatching.html')
 
    elif request.method == 'POST':
        result = {}    

        textinputs = QueryDict.copy(request.POST)
        fileinputs = ih.InputHandler.fileinputhandler(request) 

        textinputs['seq1'] = str(textinputs['seq1'].upper())
        textinputs['seq2'] = str(textinputs['seq2'].upper())

    
        result['error'] = checkEverythingBeforeAlgorithm(textinputs, fileinputs, 1) 
    
        if result['error']: 

            return render(request, template_name='main/Error pages/boyerMooreExactMatchingError.html', context={'result': result})  

        result = result = alg.Algorithms.boyer_MooreExactMatching(ih.InputHandler.inputtypehandler(textinputs, fileinputs))


        return render(request, template_name='main/Result pages/boyerMooreExactMatchingResult.html', context={'result': result})


def naiveExactMatchingPage(request):

    if 'back' in request.GET or 'back' in request.POST:
        return redirect ('/home/')

    if request.method == 'GET':
            
        if "sister" in request.GET:
            return redirect('/boyerMooreExactMatching/')
        else:
            return render(request, template_name='main/algorithms pages/naiveExactMatching.html')

    elif request.method == 'POST':


        result = {}    

        textinputs = QueryDict.copy(request.POST)
        fileinputs = ih.InputHandler.fileinputhandler(request) 

        textinputs['seq1'] = str(textinputs['seq1'].upper())
        textinputs['seq2'] = str(textinputs['seq2'].upper())
    
        result['error'] = checkEverythingBeforeAlgorithm(textinputs, fileinputs, 1) 

    
        if result['error']: 

            return render(request, template_name='main/Error pages/naiveExactMatchingError.html', context={'result': result})


        result = alg.Algorithms.naiveExactMatching(ih.InputHandler.inputtypehandler(textinputs, fileinputs) )


        return render(request, template_name='main/Result pages/naiveExactMatchingResult.html', context={'result': result})


def checkEverythingBeforeAlgorithm (textinputs, fileinputs, algotype):

    result = ''
    
    if ih.InputHandler.seqEmptyInputHandler(textinputs, fileinputs):
        result = 'Empty input, Please check the input methods and make sure that there was a given input'

    elif ih.InputHandler.multipleinputhandler(textinputs, fileinputs):
        result = 'Multiple inputs where entered, Please check the input methods and make sure that there was only one input method used per sequence input'
    
    elif algotype == 'edit distance':

        if ic.InputChecker.dnacheck(ih.InputHandler.inputtypehandler(textinputs, fileinputs)):
            result = 'One or both the sequences that you entered were not a DNA sequence, please try again after reviewing the sequences'


    elif algotype == 1:

        if ic.InputChecker.patternlengthchecker(ih.InputHandler.inputtypehandler(textinputs, fileinputs)):
            result = 'The pattern the was given is longer than the template given, please edit the pattern so it is shorter than the pattern'

        elif ic.InputChecker.dnacheck(ih.InputHandler.inputtypehandler(textinputs, fileinputs)):
            result = 'One or both the sequences that you entered were not a DNA sequence, please try again after reviewing the sequences'

    elif algotype == 2:

        if ih.InputHandler.indexAndMaxDistanceEmptyInputHandler(textinputs):
            result = 'No index was given, Plaase review the input method for index'

        elif ic.InputChecker.indexAndDistanceInputTypeChecker(textinputs):
            result = 'The index that was given was not an Integer or it was a negative Integer, please enter a positive integer' 

        elif ic.InputChecker.indexAndDistanceSignChecker(textinputs, 1):
            result = 'The index that was given is not suitable, Please enter a positive integer that is no more than the length of the pattern'

        elif ic.InputChecker.patternlengthchecker(ih.InputHandler.inputtypehandler(textinputs, fileinputs)):
            result = 'The pattern the was given is longer than the template given, please edit the pattern so it is shorter than the pattern'

        elif ic.InputChecker.dnacheck(ih.InputHandler.inputtypehandler(textinputs, fileinputs)):
            result = 'One or both the sequences that you entered were not a DNA sequence, please try again after reviewing the sequences'

    elif algotype == 3:

        if ih.InputHandler.indexAndMaxDistanceEmptyInputHandler(textinputs):
            result = 'No max distance was given, Please review the input method '

        elif ic.InputChecker.indexAndDistanceInputTypeChecker(textinputs):
            result = 'The max distance that was given was not an Integer or it was a negative Integer, please enter a positive integer' 
    
        elif ic.InputChecker.indexAndDistanceSignChecker(textinputs, 2):
            result = 'The max distnace that was given is not suitable, Please enter a positive integer that is no more than the length of the pattern'

        elif ic.InputChecker.patternlengthchecker(ih.InputHandler.inputtypehandler(textinputs, fileinputs)):
            result = 'The pattern the was given is longer than the template given, please edit the pattern so it is shorter than the pattern'

        elif ic.InputChecker.dnacheck(ih.InputHandler.inputtypehandler(textinputs, fileinputs)):
            result = 'One or both the sequences that you entered were not a DNA sequence, please try again after reviewing the sequences'

    elif algotype == 4:

        if ih.InputHandler.localAlignmentScoreEmptyInputHandler(textinputs):
            result = "There is 1 or more missing score, please review the entered scores"

        elif ic.InputChecker.localAlignmentScoreTypeChecker(textinputs):
            result = 'The scores that were given are not integers or they were negative integers, please review the entered input'

        elif ic.InputChecker.localAlignmentScoresSignChecker(textinputs):
            result = 'The scores that were given are not suitable, the match score should be a positave, the mismatch score negative and the gap score negative and lower then the mismatch score'

        elif ic.InputChecker.dnacheck(ih.InputHandler.inputtypehandler(textinputs, fileinputs)):
            result = 'One or both the sequences that you entered were not a DNA sequence, please try again after reviewing the sequences'

    elif algotype == 5:

        if ih.InputHandler.globalAlignmentScoreEmptyInputHandler(textinputs):
            result = "There is 1 or more missing score, please review the entered scores"

        elif ic.InputChecker.globalAlignmentScoreTypeChecker(textinputs):
            result = 'The scores that were given are not integers or they were negative integers, please review the entered input'    

        elif ic.InputChecker.globalAlignmentScoresSignChecker(textinputs):
            result = 'The scores that were given are not suitable, the match score should be the lowest then transition score then trasnversion score then the gap' 

        elif ic.InputChecker.dnacheck(ih.InputHandler.inputtypehandler(textinputs, fileinputs)):
            result = 'One or both the sequences that you entered were not a DNA sequence, please try again after reviewing the sequences'

    elif algotype == 6:

        seqs = ih.InputHandler.inputtypehandler(textinputs, fileinputs)

        if ic.InputChecker.hammingdistanceseqlenghchecker(seqs):
            result = 'The given sequences are not suitable, the Two sequences need to be of the same length' 

        elif ic.InputChecker.dnacheck(seqs):
            result = 'One or both the sequences that you entered were not a DNA sequence, please try again after reviewing the sequences'  

    elif algotype == 'boyer moore approximate matching':

            if ih.InputHandler.indexAndMaxDistanceEmptyInputHandler(textinputs):
                result = 'No max distance was given, Please review the input method '

            elif ic.InputChecker.indexAndDistanceInputTypeChecker(textinputs):
                result = 'The max distance that was given was not an Integer or it was a negative Integer, please enter a positive integer' 
        
            elif ic.InputChecker.indexAndDistanceSignChecker(textinputs, 2):
                result = 'The max distnace that was given is not suitable, Please enter a positive integer that is no more than the length of the pattern'

            elif ic.InputChecker.boyerMooreApproximatematchingpatternlengthchecker(ih.InputHandler.inputtypehandler(textinputs, fileinputs)):
                result = 'The pattern the was given is longer than the template given or less than 6 nucleotides, please edit the pattern and try again'

            elif ic.InputChecker.dnacheck(ih.InputHandler.inputtypehandler(textinputs, fileinputs)):
                result = 'One or both the sequences that you entered were not a DNA sequence, please try again after reviewing the sequences'  

    return result
    




    


     


    
    
        

