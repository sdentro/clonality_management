import os, stat, re
import numpy as np
from path import path
from util import merge_items

class SampleSheet(object):
    
    def __init__(self, sample2normals, sample2tumours, sample2sex, sample2normal_bam, sample2tumour_bam, sample2bb_dir, sample2variants, tumour_normal_pairs_id, tumour_normal_pairs_bam, tumour_bam2tumour_id, normal_bam2normal_id):
       # assert len(sample2sex.keys()) == len(sample2normals.keys()) and len(sample2sex.keys()) == len(sample2tumours.keys()) and \
       #     len(sample2sex.keys()) == len(sample2normal_bam.keys()) and len(sample2sex.keys()) == len(sample2tumour_bam.keys()) and \
       #     len(sample2sex.keys()) == len(sample2bb_dir.keys()) and len(sample2sex.keys()) == len(sample2variants.keys()) and \
       #     len(sample2sex.keys()) == len(tumour_normal_pairs_bam.keys()) and len(sample2sex.keys()) == len(tumour_normal_pairs_id.keys()) and \
       #     len(sample2sex.keys()) == len(normal_bam2normal_id.keys()), "SampleSheet: Received mappings do not contain all samples"
        for samplename in sample2sex.keys():
            assert samplename in sample2normals.keys() and samplename in sample2tumours.keys() and samplename in sample2normal_bam.keys() and \
            samplename in sample2tumour_bam.keys() and samplename in sample2bb_dir.keys() and samplename in sample2variants.keys() and \
            samplename in tumour_normal_pairs_id.keys() and samplename in tumour_normal_pairs_bam.keys(), "SampleSheet: Received mappings do not contain all samples"

        self._sample2normals = sample2normals
        self._sample2tumours = sample2tumours
        self._sample2sex = sample2sex
        self._sample2normal_bam = sample2normal_bam
        self._sample2tumour_bam = sample2tumour_bam
        self._sample2bb_dir = sample2bb_dir
        self._sample2variants = sample2variants
        self._tumour_normal_pairs_id = tumour_normal_pairs_id
        self._tumour_normal_pairs_bam = tumour_normal_pairs_bam
        self._tumour_bam2tumour_id = tumour_bam2tumour_id
        self._normal_bam2normal_id = normal_bam2normal_id

    def getSamplenames(self):
        return self._sample2tumours.keys()
    
    def getAllTumoursList(self):
        '''
        Returns a list of tumour names mentioned in the sample sheet
        '''
        tumours_list = list()
        for v in  self._sample2tumours.values():
            tumours_list.extend(v)
        return tumours_list
    
    def getAllNormalsList(self):
        '''
        Returns a list of normal names mentioned in the sample sheet
        '''
        normals_list = list()
        for v in  self._sample2normals.values():
            normals_list.extend(v)
        return normals_list
    
#     def getAllbbDirsList(self):
#         '''
#         Returns a list of bb_dirs mentioned in the sample sheet
#         '''
#         bb_dir_list = list()
#         for v in self._sample2bb_dir.values():
#             bb_dir_list.extend(v)
#         return bb_dir_list
    
    def getTumours(self, samplename):
        return self._sample2tumours[samplename] 
    
    def getNormals(self, samplename):
        return self._sample2normals[samplename]
    
    def getSex(self, samplename):
        return self._sample2sex[samplename]  
    
    def isMale(self, samplename):
        return (self.getSex(samplename) == 'male') or (self.getSex(samplename) == 'Male')
    
    def getTumour2NormalPairingBam(self, samplename):
        return(self._tumour_normal_pairs_bam[samplename])
    
    def getIdByTumourBam(self, bam_file):
        return(self._tumour_bam2tumour_id[bam_file])
    
    def getIdByNormalBam(self, bam_file):
        return(self._normal_bam2normal_id[bam_file])
    
    def getBbDirByTumourId(self, sample, tumourid):
        for s in self.getSamplenames():
            for ti,bb_dir in self._sample2bb_dir[s]:
                if ti==tumourid:
                    return bb_dir
        print("Warning: Did not find bb_dir in samplesheet for "+sample+", "+tumourid+" combo")
        return ""
        

def read_sample_infile(infile):
    '''
    Reads in a table: samplename\ttumour_id\ttumour_bam\tnormal_id\tnormal_bam\tbb_dir\tgender\tvariants
    Headers should be commented out with a #
    
    It checks whether a sample name is already known, if that is the case additional tumour ids and bams are 
    saved, normals not, bb_dirs not
    
    Creates and returns a SampleSheet object.
    '''
    f = open(infile, 'r')
    tumour_ids = dict()
    tumour_bam = dict()
    normal_ids = dict()
    normal_bam = dict()
    bb_dir = dict()
    sex = dict()
    variants = dict()
    tumour_normal_pairs_id = dict()
    tumour_normal_pairs_bam = dict()
    tumour_bam2tumour_id = dict()
    normal_bam2normal_id = dict()
    
    for line in f:
        l = line.strip()
        if l.startswith('#'): continue
        
        c1, c2, c3, c4, c5, c6, c7, c8 = l.split("\t")
        
        if c1 in normal_ids.keys():
            # Case add new tumour/normal to existing sample
            #print("Item "+c1+" found more than once in input file")
            #sys.exit(1)
            tumour_ids[c1] = tumour_ids[c1] + c2
            tumour_bam[c1] = tumour_bam[c1] + c3
            normal_ids[c1] = normal_ids[c1] + c4
            normal_bam[c1] = normal_bam[c1] + c5
            bb_dir[c1] = bb_dir[c1] + (c2, c6)
            sex[c1] = c7
            variants[c1] = variants[c1] + c8
            tumour_normal_pairs_id[c1] = tumour_normal_pairs_id[c1] + (c2, c4)
            tumour_normal_pairs_bam[c1] = tumour_normal_pairs_bam[c1] + (c3, c5)
        else:
            # Case new sample
            tumour_ids[c1] = [c2]
            tumour_bam[c1] = [c3]
            normal_ids[c1] = [c4]
            normal_bam[c1] = [c5]
            bb_dir[c1] = [(c2, c6)]
            sex[c1] = c7
            variants[c1] = [c8]
            tumour_normal_pairs_id[c1] = [(c2, c4)]
            tumour_normal_pairs_bam[c1] = [(c3, c5)]
            
        tumour_bam2tumour_id[c3] = c2
        normal_bam2normal_id[c5] = c4
    
    f.close()
    
    ss = SampleSheet(normal_ids, tumour_ids, sex, normal_bam, tumour_bam, bb_dir, variants, tumour_normal_pairs_id, tumour_normal_pairs_bam, tumour_bam2tumour_id, normal_bam2normal_id)
        
    return ss

def read_basic_sample_infile(infile, bb_run_dir):
    '''
    Reads in a table: samplename\ttumour_bam\tnormal_bam\tgender
    Headers should be commented out with a #
    
    It checks whether a sample name is already known, if that is the case additional tumour ids and bams are 
    saved, normals not, bb_dirs not
    
    Creates and returns a SampleSheet object.
    '''
    f = open(infile, 'r')
    tumour_ids = dict()
    tumour_bam = dict()
    normal_ids = dict()
    normal_bam = dict()
    bb_dir = dict()
    sex = dict()
    variants = dict()
    tumour_normal_pairs_id = dict()
    tumour_normal_pairs_bam = dict()
    tumour_bam2tumour_id = dict()
    normal_bam2normal_id = dict()

    for line in f:
        l = line.strip()
        if l.startswith('#'): continue
        
        # Only four of the regular columns available
        c1, c3, c5, c7 = l.split("\t")
        if (c3.endswith(".bam")):
            c2 = re.sub("\.bam", "", path(c3).basename())
            c4 = re.sub("\.bam", "", path(c5).basename())
        elif (c3.endswith(".CEL")):
            c2 = re.sub("\.CEL", "", path(c3).basename())
            c4 = re.sub("\.CEL", "", path(c5).basename())
        else:
            c2 = path(c3).basename()
            c4 = path(c5).basename()
        c6 = path.joinpath(bb_run_dir, c1)
        c8 = "placeholder"
        
        if c1 in normal_ids.keys():
            # Case add new tumour/normal to existing sample
            #print("Item "+c1+" found more than once in input file")
            #sys.exit(1)
            tumour_ids[c1] = tumour_ids[c1] + c2
            tumour_bam[c1] = tumour_bam[c1] + c3
            normal_ids[c1] = normal_ids[c1] + c4
            normal_bam[c1] = normal_bam[c1] + c5
            bb_dir[c1] = bb_dir[c1] + (c2, c6)
            sex[c1] = c7
            variants[c1] = variants[c1] + c8
            tumour_normal_pairs_id[c1] = tumour_normal_pairs_id[c1] + (c2, c4)
            tumour_normal_pairs_bam[c1] = tumour_normal_pairs_bam[c1] + (c3, c5)
        else:
            # Case new sample
            tumour_ids[c1] = [c2]
            tumour_bam[c1] = [c3]
            normal_ids[c1] = [c4]
            normal_bam[c1] = [c5]
            bb_dir[c1] = [(c2, c6)]
            sex[c1] = c7
            variants[c1] = [c8]
            tumour_normal_pairs_id[c1] = [(c2, c4)]
            tumour_normal_pairs_bam[c1] = [(c3, c5)]
            
        tumour_bam2tumour_id[c3] = c2
        normal_bam2normal_id[c5] = c4
        
    f.close()
    
    ss = SampleSheet(normal_ids, tumour_ids, sex, normal_bam, tumour_bam, bb_dir, variants, tumour_normal_pairs_id, tumour_normal_pairs_bam, tumour_bam2tumour_id, normal_bam2normal_id)
        
    return ss

def generateBsubCmd(jobname, logdir, cmd, queue="normal", mem=1, depends=None, isArray=False, threads=None):
    '''
    Transforms the cmd into a bsub command with the supplied parameters.
    '''
    bcmd = merge_items(["bsub","-q", queue, "-J \""+jobname+"\""])
    
    if isArray:
        bcmd = merge_items([bcmd, "-o", path.joinpath(logdir, jobname)+".%J.%I.out", "-e", path.joinpath(logdir, jobname+".%J.%I.err")])
    else:
        bcmd = merge_items([bcmd, "-o", path.joinpath(logdir, jobname)+".%J.out", "-e", path.joinpath(logdir, jobname+".%J.err")])

    mem = str(mem)+"000"
    bcmd = merge_items([bcmd, "-M", mem, "-R", "'span[hosts=1] select[mem>" + mem + "] rusage[mem=" + mem + "]'"])

    if depends is not None:
        depends_str = map(lambda x: "done("+x+")", depends)
        depends_str = "&&".join(depends_str)    
        bcmd = merge_items([bcmd, "-w\""+depends_str+"\""])
        
    if threads is not None:
        bcmd = merge_items([bcmd, "-n", str(threads)])

    bcmd = merge_items([bcmd, "'"+cmd+"'"])

    return(bcmd)

def writeSimpleShellScript(rundir, scriptname, cmds):
    '''
    Creates a simple script with the commands specified in cmds contained within.
    This script works with jobarrays. 
    Note: It returns the status of the last run command.
    '''
    #scriptfile = path.joinpath(rundir, 'GetAlleleFrequenciesFromBAMByChromosome'+samplename+'.sh')
    scriptfile = path.joinpath(rundir, scriptname)
    samplecommands = open(scriptfile,'w')
    samplecommands.write('#$LSB_JOBINDEX\n')
    for item in cmds:
        samplecommands.write(item+"\n")

    samplecommands.write('exit $?\n')
    samplecommands.close()
    st = os.stat(scriptfile)
    os.chmod(scriptfile, st.st_mode | stat.S_IEXEC)
    
    return(scriptfile)

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
