import multiprocessing

Y = 5

def f(val):
    return Y*val, 2

def get_value(val):
    ret = val*val
    return val, f(ret)

N = 10000
pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
outputs = pool.map(get_value, range(N))
pool.close()