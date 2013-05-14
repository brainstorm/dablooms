import pydablooms
import sys
import os
from itertools import islice
import readfq

capacity = 100000
error_rate = 0.005
kmer_size = 20

print("pydablooms version: %s" % pydablooms.__version__)

if len(sys.argv) != 4:
    sys.stderr.write("Usage: %s <bloom_file> <reference> <fastq_to_query>\n" % sys.argv[0])
    sys.exit(1)

bloom_fname = sys.argv[1]
reference = sys.argv[2]
fastq = sys.argv[3]

bloom = pydablooms.Dablooms(capacity=capacity,
                           error_rate=error_rate,
                           filepath=bloom_fname)

# http://scipher.wordpress.com/2010/12/02/simple-sliding-window-iterator-in-python/
def sliding_window(sequence, winSize, step=1):
    """Returns a generator that will iterate through
    the defined chunks of input sequence.  Input sequence
    must be iterable."""
 
    # Verify the inputs
    try: it = iter(sequence)
    except TypeError:
        raise Exception("**ERROR** sequence must be iterable.")
    if not ((type(winSize) == type(0)) and (type(step) == type(0))):
        raise Exception("**ERROR** type(winSize) and type(step) must be int.")
    if step > winSize:
        raise Exception("**ERROR** step must not be larger than winSize.")
    if winSize > len(sequence):
        raise Exception("**ERROR** winSize must not be larger than sequence length.")
 
    # Pre-compute number of chunks to emit
    numOfChunks = ((len(sequence)-winSize)/step)+1
 
    # Do the work
    for i in range(0,numOfChunks*step,step):
        yield sequence[i:i+winSize]

with open(reference, 'rb') as fh:
    for _, seq, _ in readfq.readfq(fh):
        for i, kmer in enumerate(sliding_window(seq, kmer_size)):
            bloom.add(kmer, i)

i = 0
with open(reference, 'rb') as fh:
    for _, seq, _ in readfq.readfq(fh):
        for kmer in sliding_window(seq, kmer_size):
            if i % 5 == 0:
                bloom.delete(kmer, i)
            i = i + 1


bloom.flush()
del bloom

bloom = pydablooms.load_dabloom(capacity=capacity,
                                error_rate=error_rate,
                                filepath=bloom_fname)

true_positives = 0
true_negatives = 0
false_positives = 0
false_negatives = 0

i = 0
with open(reference) as fh:
    for _, seq, _ in readfq.readfq(fh):
        for kmer in sliding_window(seq, kmer_size):
            exists = bloom.check(kmer)
            contains = kmer in bloom
            assert exists == contains, \
                "ERROR: %r from 'bloom.check(x)', %i from 'x in bloom'" \
                % (exists, contains)
    
            if i % 5 == 0:
                if exists:
                    false_positives += 1
                else:
                    true_negatives += 1
            else:
                if exists:
                    true_positives += 1
                else:
                    false_negatives += 1
                    sys.stderr.write("ERROR: False negative: '%s'\n" % seq.rstrip())
            i += 1

del bloom

false_discovery_rate = float(false_positives) / (false_positives + true_positives)

print('''
Elements Added:   %6d
Elements Removed: %6d

True Positives:   %6d
True Negatives:   %6d
False Positives:  %6d
False Negatives:  %6d

False discovery rate: %.4f
Total Size: %d KiB''' % (
                         i, i/5,
                         true_positives,
                         true_negatives,
                         false_positives,
                         false_negatives,
                         false_discovery_rate,
                         os.stat(bloom_fname).st_size / 1024
                        )
)

if false_negatives > 0:
    print("TEST FAIL (false negatives exist)")
elif false_discovery_rate > error_rate:
    print("TEST WARN (false discovery rate too high)")
else:
    print("TEST PASS")
print("")

if false_negatives > 0:
    sys.exit(1)
else:
    sys.exit(0)
