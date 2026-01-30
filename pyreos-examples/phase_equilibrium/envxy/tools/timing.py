import time

def timing(f):
    
    def wrap(*args, **kwargs):

        time1 = time.time()
        ret = f(*args, **kwargs)
        time2 = time.time()
        dt = (time2-time1)

        return ret, dt
    
    return wrap
