#
# basic_math.py
#

def sign(perm):
    copy = perm[:]
    copy.sort()
    assert copy == range(len(perm))
    out = 1
    for i in range(len(perm)):
        for j in range(i):
            if perm[i] < perm[j]:
                out *= -1
    return out

if __name__ == '__main__':
    print sign([0,1,2,3])
    print sign([0,2,1,3])
    print sign([0,1,2,2])
