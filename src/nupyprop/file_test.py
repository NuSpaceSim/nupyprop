import numpy as np


def interploation(x,fname):
    
    e_format = "{:5.2f}"
    result = (x + 1)
    result2 = result/2.0
    result_tup = (x,result,result2)
    fname.write(str(result_tup[2]) + '\n')
    
    return result_tup



x = np.linspace(0,100,101)

file_name = 'testing.dat'
# result is going to be a tuple (1,2,1.5)
with open(file_name, 'w') as e_file:
    
    for i in range(len(x)):
        result = interploation(x[i], e_file)

print(result[0],result[1])