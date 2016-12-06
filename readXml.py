import xml.etree.ElementTree as ET

class FilterSequences:
  sequences = []
  alignOrgIndex = []
  alignIndex = []
  
  def __init__(self, inputXml):
    self.xmlFile = inputXml 

  def add(self, seqData, orgAlignIdx, sIndex):
    self.sequences.append(seqData.replace("-",""))
    self.alignIndex.append(sIndex)
    self.alignOrgIndex.append(orgAlignIdx)

  def parseXML1(self):
     tree = ET.parse(self.xmlFile)
     root = tree.getroot()
     for sequence in root.find('alignment').findall('sequence'):
          sequenceData = str(sequence.find('seq-data').text)
          is_first = 1
          without_gap = 0
          for alignment in sequence.iter('fitem'):
              startIndex = int(alignment.find('fstart').text)
              stopIndex = int(alignment.find('fstop').text)
              if sequenceData[startIndex : stopIndex].find('-') != -1: 
                 without_gap = 1
              elif is_first == 1:
                 is_first = 0
                 rcrdStart = startIndex
                 adjustIndex = sequenceData.count('-',startIndex,stopIndex)
                 
          if without_gap == 0:
              self.add(sequenceData, rcrdStart, rcrdStart - adjustIndex)  

  def parseXML2(self):
     tree = ET.parse(self.xmlFile)
     root = tree.getroot()
     alignStartIndex = []
     alignEndIndex = []
     alignment = root.find('alignment').find('column-score').find('colsco-data').text
     alignment = alignment.split(' ')
     algn = count = 0
     prev = -1
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
               alignEndIndex.append(count)
          prev = alignment[i]

     for sequence in root.find('alignment').findall('sequence'):
          sequenceData = str(sequence.find('seq-data').text)
          is_first = 1
          without_gap = 0
          rcrdStart = adjustIndex = 0
          for i in range(len(alignStartIndex)):
              startIndex = alignStartIndex[i]
              stopIndex = alignEndIndex[i]
              if sequenceData[startIndex : stopIndex].find('-') != -1:
                 without_gap = 1
              elif is_first == 1:
                 is_first = 0
                 adjustIndex = sequenceData.count('-', rcrdStart, startIndex)
                 rcrdStart = startIndex

          if without_gap == 0:
              self.add(sequenceData, rcrdStart, rcrdStart - adjustIndex)

  def printSequences(self):
     for i in range(len(self.sequences)):
       print('{0} is now aligned at {1} than {other}'.format(self.sequences[i], \
              self.alignIndex[i], other = self.alignOrgIndex[i]))  


def main():
    Filter = FilterSequences("data/BB11002.xml")
    Filter.parseXML2()
    Filter.printSequences()

main()
