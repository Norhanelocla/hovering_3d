
from cmath import sin
import ctypes
from operator import mod
import numpy
import glob
import math

# https://nesi.github.io/perf-training/python-scatter/ctypes
# https://stackoverflow.com/questions/14887378/how-to-return-array-from-c-function-to-python-using-ctypes
# find the shared library, the path depends on the platform and Python version
libfile = glob.glob('out/build/linux-default/libinterface_simulink.so')[0]

input_length=3
output_length=12

# 1. open the shared library
mylib = ctypes.CDLL(libfile)

# 2. tell Python the argument and result types of function mysum
mylib.LoadModel()
mylib.InitModel()

mylib.run_a_step.restype = numpy.ctypeslib.ndpointer(dtype=numpy.float64,shape=(output_length,))
mylib.run_a_step.argtypes = [numpy.ctypeslib.ndpointer(dtype=numpy.float64),ctypes.c_int]
inputs = numpy.zeros(input_length,dtype=numpy.float64)
inputs[0]=0.5
inputs[1]=0.5
#inputs[2]=0.2


mylib.set_model_parameters.argtypes = [ctypes.c_double,ctypes.c_int]
null=mylib.set_model_parameters(5.0,0)
null=mylib.set_model_parameters(5.0,1)
null=mylib.set_model_parameters(-5.0,2)
null=mylib.set_model_parameters(-5.0,3)
null=mylib.set_model_parameters(1,4)
null=mylib.set_model_parameters(1,5)


time_step=0.0005 # in s

sim_time=50 # in s 

sim_num_steps=int(sim_time/time_step)

# 3. call function mysum
for i in range(sim_num_steps):
    #inputs[2]=float(math.sin(time_step*i))
    outputs = mylib.run_a_step(inputs,len(inputs))
    if i % 1000 ==1:
        print('x pos: ' + str(outputs[0]))
        print('x vel: ' + str(outputs[3]))

#print('sum of array: {}'.format(numpy.ctypeslib.as_array(outputs,(1,))))