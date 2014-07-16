import numpy as np

def read_item_list(infile):
    '''
    Reads in a list of items specified one per line
    '''
    f = open(infile, 'r')

    items = list()

    for line in f:
        l = line.strip()
        items.append(l)

    return(items)

def match_tumour_normal(samplenames, tumour_id, normal_id):
    '''
    Match the tumour and normal samples up. This currently only works when the supplied samplenames
    end with the tumour_id or the normal_id.
    
    Returns a list of tuples where the first item is the tumour, second item is the normal.
    '''
    samplenames = np.array(samplenames)
    tumours = samplenames[np.array([tumour_id in item for item in samplenames])]
    normals = samplenames[np.array([normal_id in item for item in samplenames])]
    samples = np.array([item.strip(normal_id) for item in normals])
    
    matches = []
    for sample in samples:
        t = tumours[np.array([(sample in item) for item in tumours])][0]
        n = normals[np.array([(sample in item) for item in normals])][0]
        matches.append((t,n))
        
    return matches
    
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