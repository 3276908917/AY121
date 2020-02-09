exec(open('synthesis.py').read())                                            
redos = load_saves('redos.npz')
redos = load_saves('data/redos.npz')
r1 = redos['re1']                                                            
r1c = (r1[100:]).astype(float)
p1 = power_plot(r1c, 658) 
i1 = invf(p1)
