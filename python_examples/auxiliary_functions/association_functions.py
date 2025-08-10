import numpy as np
def get_non_bondend_sites_from_states(STATES):

    N=len(STATES)

    return np.array([STATES[i].non_bonded_sites() for i in range(N)])






        
