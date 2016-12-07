import xml.etree.ElementTree as ET

class FilterSequences:
  
  def __init__(self, inputXml):
     self.sequences = []
     self.AlignmentEntry = []
     self.xmlFile = inputXml 
     self.parseXML()

  def parseXML(self):
     tree = ET.parse(self.xmlFile)
     root = tree.getroot()
     alignStartIndex = []
     alignLen = []
     alignment = root.find('alignment').find('column-score').find('colsco-data').text.strip().split(' ')

     intron = 0
     exon = 1
     state = intron
     for i in range(len(alignment)):
         val = alignment[i]
         if val == '1':
             if state == intron:
                 alignStartIndex.append(i)
             state = exon
         elif val == '0':
             if state == exon:
                 alignLen.append(i - alignStartIndex[-1])
             state = intron
     if state == exon:
         alignLen.append(i - alignStartIndex[-1])

     #print(alignStartIndex)
     #print(alignLen)
     
     add_sequence = 1
     startIndex = stopIndex = 0
     for i in range(len(alignStartIndex)):
        multi_alignment = []
        if i != 0:
           add_sequence = 0      
        for sequence in root.find('alignment').findall('sequence'):
           sequenceData = str(sequence.find('seq-data').text)
           startIndex = alignStartIndex[i]
           stopIndex = startIndex + alignLen[i]
           
           if sequenceData[startIndex : stopIndex].find('-') != -1:
               multi_alignment.append(None)
           else:
               adjustIndex = sequenceData.count('-', 0, startIndex)
               multi_alignment.append(startIndex - adjustIndex)
           if add_sequence == 1:
               self.sequences.append(sequenceData.replace("-","").upper().strip())
        #print stopIndex
        #print startIndex
        self.AlignmentEntry.append((alignLen[i], multi_alignment))
               
    
  def printSequences(self):
     for i in range(len(self.sequences)):
       print('Sequence {0}:{1} '.format(i, self.sequences[i]))
    
     for length,indices in self.AlignmentEntry :
       print('len {0}'.format(length))
       print('alignment{0}'.format(indices))
       for index, sequence in zip(indices, self.sequences):
           if index is not None:
               print(sequence[index:index+length])

       
def test():
     filename = "data/BBS12034.xml"
     Filter = FilterSequences(filename)
     Filter.printSequences()         


#test()
