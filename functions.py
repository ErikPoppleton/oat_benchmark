from collections import namedtuple
import os
import numpy as np
import mmap

# chunk file into blocks of given size
Chunk = namedtuple('Chunk', ['block','offset', 'is_last','file_size'])
def blocks(file, fsize, size=10000000):
    current_chunk = 0  
    while True:
        b = file.read(size)
        if not b: break
        yield Chunk(b,current_chunk*size, current_chunk * size + size > fsize, fsize)
        current_chunk+=1

###############################################################################
#      Different implementations of the search for configuration starts       #
###############################################################################

# Current and mmap
# mmap is only faster when the chunk size is small, which is less efficient overall.
def find_all(a_str, sub):
    #https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring
    start = 0
    idxs = []
    while True:
        start = a_str.find(sub, start)
        if start == -1: return idxs
        idxs.append(start)
        start += len(sub) # use start += 1 to find overlapping matches

# Previous one started searching on the very next byte from the last start
# That's the responsible thing to do, but we know that oxDNA files have some amount of stuff between configuration starts
# So we can begin the search from something that we predict might be ~ 1/2 of the way into the next configuration and can search a lot less.
# This turns out to be very negligiblly faster
def find_all_greedy(data, substr):
    start = 0
    idxs = []
    while True:
        start = data.find(substr, start)
        if start == -1: return idxs
        idxs.append(start)
        if len(idxs) > 2:
            start += int((idxs[-1] - idxs[-2])/2)
        else:
            start += len(substr)

def find_all_bytearray(data, substr):
    pass

###############################################################################
#                              Indexing methods                               #
###############################################################################

#current 
ConfInfo = namedtuple('ConfInfo', ['offset', 'size','id'])
def index(traj_file, chunk_size): 
    val = b"t" 
    conf_starts = []
    counter = 0
    fsize = os.stat(traj_file).st_size
    for chunk in blocks(open(traj_file, 'rb'), fsize, size=chunk_size):
        idxs = np.array(find_all(chunk.block,val)) + chunk.offset # find all offsets of t in file
        conf_starts.extend(idxs)

    conf_starts = np.array(conf_starts)
    #generate a list of confs for all but the last one
    idxs = [ConfInfo(conf_starts[i], conf_starts[i+1] - conf_starts[i],i) 
                                            for i in range(len(conf_starts)-1)]
    #handle last offset
    idxs.append(ConfInfo(conf_starts[-1], fsize - conf_starts[-1], len(conf_starts)-1))
    return idxs

#greedy indexer used
def index_greedy(traj_file, chunk_size): 
    val = b"t" 
    conf_starts = []
    counter = 0
    fsize = os.stat(traj_file).st_size
    for chunk in blocks(open(traj_file, 'rb'), fsize, size=chunk_size):
        idxs = np.array(find_all_greedy(chunk.block,val)) + chunk.offset # find all offsets of t in file
        conf_starts.extend(idxs)

    conf_starts = np.array(conf_starts)
    #generate a list of confs for all but the last one
    idxs = [ConfInfo(conf_starts[i], conf_starts[i+1] - conf_starts[i],i) 
                                            for i in range(len(conf_starts)-1)]
    #handle last offset
    idxs.append(ConfInfo(conf_starts[-1], fsize - conf_starts[-1], len(conf_starts)-1))
    return idxs

#mmap 
def index_mmap(traj_file, chunk_size): 
    val = b"t" 
    conf_starts = []
    counter = 0
    fsize = os.stat(traj_file).st_size
    with open(traj_file, 'rb') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        for chunk in blocks(mm, fsize, size=chunk_size):
            idxs = np.array(find_all(chunk.block,val)) + chunk.offset # find all offsets of t in file
            conf_starts.extend(idxs)

    conf_starts = np.array(conf_starts)
    #generate a list of confs for all but the last one
    idxs = [ConfInfo(conf_starts[i], conf_starts[i+1] - conf_starts[i],i) 
                                            for i in range(len(conf_starts)-1)]
    #handle last offset
    idxs.append(ConfInfo(conf_starts[-1], fsize - conf_starts[-1], len(conf_starts)-1))
    return idxs

###############################################################################
#                                Parsing                                      #
###############################################################################


def parse_conf(lines,nbases):
    lines = lines.split('\n')
    if lines[-1] == '' and len(lines) -4 != nbases:
        raise Exception("Incorrect number of bases in topology file")
    elif lines[-1] != '' and len(lines) -3 != nbases:
        raise Exception("Incorrect number of bases in topology file")   
    #setup dummy conf 
    conf = base_array(
        0, np.zeros(3), 0,
        np.zeros([nbases, 3], dtype=float),
        np.zeros([nbases, 3], dtype=float),
        np.zeros([nbases, 3], dtype=float),
    )
    # populate our dummy conf by data from the conf
    conf.time = float(lines[0][lines[0].index("=")+1:])
    conf.box = np.array(lines[1].split("=")[1].split(), dtype=float)
    conf.energy = np.array(lines[2].split("=")[1].split(), dtype=float)
    # parse out the pos and a's 
    for i in range(nbases):
        line = lines[3+i].split()
        conf.positions[i] = np.array(line[0:3], dtype=float)
        conf.a1s[i] = np.array(line[3:6], dtype=float)
        conf.a3s[i] = np.array(line[6:9], dtype=float)
    return conf