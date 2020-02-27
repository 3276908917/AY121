## Change these values as you run
offset = 1600000
sample_size = 16000

def bit(chunk):
    c = []
    for i in range(0, offset, sample_size):
        real = chunk[i: i + sample_size]
        imag = chunk[offset + i: offset + i + sample_size]
        c.append(real + 1j*imag)
    return c

# e.g. off = load_saves('../on.npz')['raw_off']
def pipeline(case, reduction):
    compressor = []
    for a in case:
        c_chunk = bit(a)
        P = power_barrage(c_chunk)
        compressor.append(reduction(P))
    return reduction(compressor)

def attempt():
    return samples_median, power_barrage
