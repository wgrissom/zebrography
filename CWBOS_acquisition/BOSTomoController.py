'''
Control background pattern on iPad by sending commands to BOSTomoDisplay App remotely.

Will Grissom, Vanderbilt University, 2018
'''

import socket
import struct
import sys

pattern = str(sys.argv[1])
pattern = pattern.encode('utf-8')
x = int(sys.argv[2])
y = int(sys.argv[3])
z = int(sys.argv[4])
k = int(sys.argv[5])
ip = str(sys.argv[6])
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((ip, 5002))
s.sendall(struct.pack('ciiii',pattern,x,y,z,k))
data = s.recv(32)
foo, bar, hi, bop, eoi = struct.unpack('ciiii',data)
s.close()
print ('Received'+repr(foo)+repr(bar)+repr(hi)+repr(bop)+repr(eoi))
