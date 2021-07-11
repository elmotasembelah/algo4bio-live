
class ReadingFiles:

    @staticmethod
    def fastafilereader (request, key):
        
        file = (request.FILES[key].file)
        seq = str(file.read())
        seq = seq [2: len(seq)-1]
        if seq[0] == '>':
            seq = seq[seq.find('\\'): len(seq)]
        
        checkAnotherSeq = seq.find('>')
        if checkAnotherSeq != -1:
            seq = seq[:checkAnotherSeq]

        seq = seq.replace('\\r\\n', '')
        seq = seq.replace(' ','')
        
        return seq


            
