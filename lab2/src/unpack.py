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

def thermal():
    block = []
    for file in glob.glob("./*"):
        samples = []
        reader = open(file, 'r')
        for line in reader:
            trim = line.strip()
            samples.append(float(trim))
        block.append(samples)
    return block

def unpickle_folder():
    re_block, im_block = [], []
    for file in glob.glob("./*"):
        data = pickle.load(open(file, "rb"))
        re_block.append(data['real'])
        im_block.append(data['imaginary'])
    return re_block, im_block

