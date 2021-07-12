import numpy as np
import ctypes

LIBRARY_FILE = "./libReweight.so"
libf         = ctypes.CDLL(LIBRARY_FILE)

emat         = np.array([[0, 1, 0, 0, 0], [1, 0, 0, 0.5, 0], [0, 0, 0, 1, 0], [0, 0, 1, 0, 0], [0.5, 0, 0, 0, 1]], dtype = np.double)
edges        = np.array([[1, 2], [4, 2], [4, 3]], dtype = np.intc)
out          = np.zeros((3, ), dtype = np.double)
emat_f       = emat.flatten()
edges_f      = edges.flatten()

libf.predict_edge(ctypes.c_void_p(emat_f.ctypes.data),
                  ctypes.c_void_p(edges_f.ctypes.data),
                  ctypes.c_int(5),
                  ctypes.c_int(3),
                  ctypes.c_void_p(out.ctypes.data))

print(out)
