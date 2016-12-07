import xml.etree.ElementTree as ET

class FilterSequences:
  sequences = []
  AlignmentEntry = []
  
  def __init__(self, inputXml):
    self.xmlFile = inputXml 

  def add_entry(self, alignmentSize, malign):
     entry = []
     entry.append(alignmentSize)
     entry.append(malign)
     self.AlignmentEntry.append(entry)

  def parseXML(self):
     tree = ET.parse(self.xmlFile)
     root = tree.getroot()
     alignStartIndex = []
     alignLen = []
     alignment = root.find('alignment').find('column-score').find('colsco-data').text
     alignment = alignment.split(' ')
     algn = count = 0
     prev = '-1'
     for i in range(len(alignment)):
          if alignment[i] == '1':
             if algn == 1:
               alignStartIndex.append(i)
               count = 0
               algn = 0
             count = count + 1 
          elif alignment[i] == '0':
             algn = 1   
             if prev == '1':
               alignLen.append(count)
          prev = alignment[i]
     
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
               multi_alignment.append(startIndex - adjustIndex + 1)
           if add_sequence == 1:
               self.sequences.append(sequenceData.replace("-",""))   
        print stopIndex 
        print startIndex
        self.add_entry(alignLen[i], multi_alignment)
               
    
  def printSequences(self):
     for i in range(len(self.sequences)):
       print('Sequqnece {0}:{1} '.format(i, self.sequences[i]))
    
     for j in self.AlignmentEntry :
       print('len {0}'.format(j[0]))
       print('alignment{0}'.format(j[1]))
       
def main():
     Filter = FilterSequences("data/BB11002.xml")
     Filter.parseXML()
     Filter.printSequences()         
