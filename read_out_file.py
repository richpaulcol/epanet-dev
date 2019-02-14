import numpy as np 
import pylab as pp
from ast import literal_eval
from struct import *

No_nodes = 3
No_links = 2
No_tanks = 0
No_pumps = 0
No_periods = 3

f = open("Example1.out",'rb')
Prolog = f.read(84)
EnergyUse = f.read(4)
EPS = f.read()

PrologData = []

for i in range(21):
	PrologData.append(unpack('i',Prolog[4*i:4*i+4])[0])
NoNodes = PrologData[6]
NoLinks = PrologData[7]
NoPumps = PrologData[8]

NodeDataSize = PrologData[-3]
LinkDataSize = PrologData[-2]
PumpDataSize = PrologData[-1]

dT = PrologData[17]

print 'Dave'


EPS_data_points = len(EPS)/4
Data = []
for i in range(len(EPS)/4):
	Data.append( unpack('f',EPS[4*i:4*i+4]))
Data = np.array(Data)

DataTimeSize = NoNodes*NodeDataSize + NoLinks*LinkDataSize + NoPumps*PumpDataSize
Times = EPS_data_points / DataTimeSize
Data = Data.reshape(Times,DataTimeSize)

TimeList = np.linspace(0,Times*60,Times)/(60*60)
NodalData = {}
for i in range(NoNodes):
	NodalData[str(i)] = Data[:,i*NodeDataSize:i*NodeDataSize+NodeDataSize]

fig,axs = pp.subplots(3,1,sharex=True)
axs[0].plot(TimeList,NodalData['1'][:,2])
axs[0].plot(TimeList,NodalData['0'][:,2])

axs[1].plot(TimeList,NodalData['1'][:,0])
axs[1].plot(TimeList,NodalData['0'][:,0])



axs[2].plot(TimeList,NodalData['1'][:,5])
axs[2].plot(TimeList,NodalData['0'][:,5])

axs[0].set_ylabel('Flow')
axs[1].set_ylabel('Head')
axs[2].set_ylabel('Turbidity')
pp.show()
### Node:: Head, Pressure, demand, balance, outflow, quality
#print unpack('hh',f.read(4))

