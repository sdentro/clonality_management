import sys
import numpy as np

class SampleSheet(object):
    
    def __init__(self, sample2normals, sample2tumours, sample2sex):
        assert len(sample2sex.keys()) == len(sample2normals.keys()) and len(sample2sex.keys()) == len(sample2tumours.keys()), "SampleSheet: Received mappings do not contain all samples"
        for samplename in sample2sex.keys():
            assert samplename in sample2normals.keys() and samplename in sample2tumours.keys(), "SampleSheet: Received mappings do not contain all samples"
        
        self._sample2normals = sample2normals
        self._sample2tumours = sample2tumours
        self._sample2sex = sample2sex

    def getSamplenames(self):
        return self._sample2tumours.keys()
    
    def getAllTumoursList(self):
        '''
        Returns a list of tumour names mentioned in the sample sheet
        '''
        tumours_list = list()
        for v in  self._sample2tumours.values():
            tumours_list.append(v)
        return tumours_list
    
    def getAllNormalsList(self):
        '''
        Returns a list of normal names mentioned in the sample sheet
        '''
        normals_list = list()
        for v in  self._sample2normals.values():
            normals_list.append(v)
        return normals_list
    
    def getTumours(self, samplename):
        return self._sample2tumours[samplename] 
    
    def getNormals(self, samplename):
        return self._sample2normals[samplename]
    
    def getSex(self, samplename):
        return self._sample2sex[samplename]  
    
    def isMale(self, samplename):
        return (self.getSex(samplename) == 'male') or (self.getSex(samplename) == 'Male')

def read_sample_infile(infile):
    '''
    Reads in a table: samplename\tmale\tnormalsubsample1,normalsubsample2,etc\ttumoursubsample1,tumoursubsample2,etc
    
    Creates and returns a SampleSheet object.
    '''
    f = open(infile, 'r')
    
    normals = dict()
    tumours = dict()
    sex = dict()
    
    for line in f:
        l = line.strip()
        if l.startswith('#'): continue
        
        c1, c2, c3, c4 = l.split("\t")
        
        if c1 in normals.keys():
            print("Item "+c1+" found more than once in input file")
            sys.exit(1)
        sex[c1] = c2
        normals[c1] = c3.split(",")
        tumours[c1] = c4.split(",")
    
    f.close()
    
    ss = SampleSheet(normals, tumours, sex)
        
    return ss

# def read_item_list(infile):
#     '''
#     Reads in a list of items specified one per line
#     '''
#     f = open(infile, 'r')
# 
#     items = list()
# 
#     for line in f:
#         l = line.strip()
#         items.append(l)
# 
#     return(items)

# def match_tumour_normal(samplenames, tumour_id, normal_id):
#     '''
#     Match the tumour and normal samples up. This currently only works when the supplied samplenames
#     end with the tumour_id or the normal_id.
#     
#     Returns a list of tuples where the first item is the tumour, second item is the normal.
#     '''
#     samplenames = np.array(samplenames)
#     tumours = samplenames[np.array([tumour_id in item for item in samplenames])]
#     normals = samplenames[np.array([normal_id in item for item in samplenames])]
#     samples = np.array([item.strip(normal_id) for item in normals])
#     
#     matches = []
#     for sample in samples:
#         t = tumours[np.array([(sample in item) for item in tumours])][0]
#         n = normals[np.array([(sample in item) for item in normals])][0]
#         matches.append((t,n))
#         
#     return matches
    
def match_sample_to_file(samplenames, list_of_files):
    '''
    Matches the samplenames to files listed in the list_of_files. Samples with no files
    are not mentioned in the mapping. The second returned item is a list of samples that
    were not mapped
    '''
    list_of_files = np.array(list_of_files)
    mapping = dict()
    unmapped = []
    for sample in samplenames:
        files_for_sample = np.array([sample in infile for infile in list_of_files])
        if (sum(files_for_sample) > 0):
            mapping[sample] = list_of_files[files_for_sample]
        else:
            unmapped.append(sample)
    return mapping,unmapped