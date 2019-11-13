import numpy as np
import math 

'''
ROwf/zero,zero,0.06562d0/
RH1wf/0.7557d0,zero,-0.5223d0/
RH2wf/-0.7557d0,zero,-0.5223d0/
RMwf/0.d0,0.d0,-0.08438d0/
'''

angleHOH=104.52*math.pi/180.0
dOH=0.9572 #experiment
dOM=0.1546
zO=0.0 #fixed
zM=-0.1546
zH=zO-math.sqrt(0.5*dOH*dOH*(1.0+math.cos(angleHOH)))
xH=math.sqrt(dOH*dOH-(zO-zH)*(zO-zH))
a = np.array([ xH, 0.0, zH])
b = np.array([0.0, 0.0, zO])
c = np.array([-xH, 0.0, zH])
d = np.array([0.0, 0.0, zM])

ba = b-a
bc = b-c
bd = b-d

cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
angle = np.arccos(cosine_angle)*180.0/math.pi

print("d_OH = "+str(np.linalg.norm(ba)))
print("d_OM = "+str(np.linalg.norm(bd)))
print(angle)
print(a)
print(b)
print(c)
print(d)
