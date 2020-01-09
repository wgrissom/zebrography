'''

On the experiment computer, 
control background pattern on iPad by sending commands to BOSTomoDisplay APP remotely.

Will Grissom, Vanderbilt University, 2018
'''

import socket
import struct
import sys


# pattern: "p": pins; "s": Stripes; "r": RGB; "R": Red; "G": Green; "B": Blue; "A": Random
# x: pattern width
# y: pattern spacing
# z: offset
# k: if the black background is needed. 0: White background, 1: Black background.
# ip: IP address displayed on the iPad 

pattern = str(sys.argv[1])
pattern = pattern.encode('utf-8')
x = int(sys.argv[2]) 
y = int(sys.argv[3])
z = int(sys.argv[4])
k = int(sys.argv[5])
ip = str(sys.argv[6])
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((ip, 5002)) # create the connection to the iPad. 
s.sendall(struct.pack('ciiii',pattern,x,y,z,k)) # send commands to iPad to switch the background pattern.
data = s.recv(32)
foo, bar, hi, bop, eoi = struct.unpack('ciiii',data) # interpret data received by iPad and check if it works well. 
s.close()
print ('Received'+repr(foo)+repr(bar)+repr(hi)+repr(bop)+repr(eoi))
