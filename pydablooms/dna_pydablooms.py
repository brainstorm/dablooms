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

def window(seq, n):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

i = 0
with open(reference, 'rb') as fh:
    for name, seq, qual in readfq.readfq(fh):
        i = 0;
        while(i <= (len(seq) - kmer_size)):
            # adding k_mers to a bloom filter with sliding window.
            for w in window([seq], kmer_size):
                print w
            #bloom.add(w, i)
            #i += 1

i = 0
with open(reference, 'rb') as fh:
    for name, seq, qual in readfq.readfq(fh):
        if i % 5 == 0:
            bloom.delete(seq.rstrip(), i)
        i += 1


bloom.flush()
del bloom

bloom = pydablooms.load_dabloom(capacity=capacity,
                                error_rate=error_rate,
                                filepath=bloom_fname)

true_positives = 0
true_negatives = 0
false_positives = 0
false_negatives = 0

with open(fastq) as fh:
    for name, seq, qual in readfq.readfq(fh):
        exists = bloom.check(seq.rstrip())
        contains = seq.rstrip() in bloom
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

false_positive_rate = float(false_positives) / (false_positives + true_negatives)

print('''
Elements Added:   %6d
Elements Removed: %6d

True Positives:   %6d
True Negatives:   %6d
False Positives:  %6d
False Negatives:  %6d

False positive rate: %.4f
Total Size: %d KiB''' % (
                         i, i/5,
                         true_positives,
                         true_negatives,
                         false_positives,
                         false_negatives,
                         false_positive_rate,
                         os.stat(bloom_fname).st_size / 1024
                        )
)

if false_negatives > 0:
    print("TEST FAIL (false negatives exist)")
elif false_positive_rate > error_rate:
    print("TEST WARN (false positive rate too high)")
else:
    print("TEST PASS")
print("")

if false_negatives > 0:
    sys.exit(1)
else:
    sys.exit(0)
