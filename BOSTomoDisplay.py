'''
An app for displaying background patterns for BOS tomography
Will Grissom, Vanderbilt University, 2018
'''

import ui
import console
import numpy as np
from PIL import Image
import io
import socket
import select
import struct

def get_local_ip():
	s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
	try:
		s.connect(('google.com', 80))
		ip = s.getsockname()[0]
		s.close()
	except:
		ip = 'N/A'
	return ip

def pil2ui(pil_img):
    with io.BytesIO() as buffer:
        pil_img.save(buffer, format='PNG')
        return ui.Image.from_data(buffer.getvalue())

class SketchView (ui.View):
	def __init__(self, width=1024, height=1024):
		self.bg_color = 'white'
		self.imgWidth = width
		self.imgHeight = height
		self.imgEmpty = Image.fromarray(np.zeros([width,height]), 'L')
		self.imgFull = Image.fromarray(np.ones([width,height],dtype=np.uint8)*255, 'L')
		iv = ui.ImageView(frame=(0,0,width,height), image=pil2ui(self.imgEmpty))
		self.add_subview(iv)
		red_button = ui.ButtonItem() # all red
		red_button.title = 'Red'
		red_button.tint_color = 'red'
		red_button.action = self.togglePattern
		green_button = ui.ButtonItem() # all green
		green_button.title = 'Green'
		green_button.tint_color = 'green'
		green_button.action = self.togglePattern
		blue_button = ui.ButtonItem() # all blue
		blue_button.title = 'Blue'
		blue_button.tint_color = 'blue'
		blue_button.action = self.togglePattern
		#white_button = ui.ButtonItem()
		#white_button.title = 'White'
		#white_button.tint_color = 'black'
		#white_button.action = self.togglePattern
		random_button = ui.ButtonItem() # random black + white
		random_button.title = 'Random'
		random_button.action = self.togglePattern
		RGB_button = ui.ButtonItem() # RGB pins, black or white background
		RGB_button.title = 'RGB'
		RGB_button.action = self.togglePattern
		pins_button = ui.ButtonItem() # white pins on black background or vice versa
		pins_button.title = 'Pins'
		pins_button.action = self.togglePattern
		stripes_button = ui.ButtonItem() # white stripes on black background or vice versa
		stripes_button.title = 'Stripes'
		stripes_button.action = self.togglePattern
		self.right_button_items = [blue_button, green_button, red_button, RGB_button, pins_button, stripes_button, random_button]
		self.image_view = iv
		self.blackBackground = True
		self.pattern = 'Pins' # white pins, black or white background
		self.pinWidth = 1
		self.pinSpacing = 4
		self.stripeWidth = 2
		self.stripeSpacing = 8
		self.offset = 2
		self.keepServerOn = True;
		self.updatePattern()
		#self.updateImg()
		
	def togglePattern(self, sender):
		if self.pattern == sender.title:
			self.blackBackground = not self.blackBackground # toggle black background
		self.pattern = sender.title;
		self.updatePattern()
		#self.updateImg()
		#self.setPattern(self, sender.title)
		
	def updateImg(self):	
		#if (self.blackBackground == False) & (self.patternColor == 'White'):
		#	imgPattern = Image.fromarray(255-self.mtxPattern, 'L')
		#	self.imgMerged = Image.merge('RGB',(imgPattern,imgPattern,imgPattern))
		#else:
		#	if self.patternColor == 'White':
		#		self.imgMerged = Image.merge('RGB', (self.imgPattern,self.imgPattern,self.imgPattern))	
		self.image_view.image = pil2ui(self.imgMerged)
			
	def updatePattern(self):
		stripeWidth = np.int(self.stripeWidth)
		stripeSpacing = np.int(self.stripeSpacing)
		if self.pattern == 'Random': # random black and white
			mtxPattern = np.array(np.random.random((self.imgWidth,self.imgHeight))*2, dtype=np.uint8)*255
			imgPattern = Image.fromarray(mtxPattern, 'L')
			self.imgMerged = Image.merge('RGB',(imgPattern,imgPattern,imgPattern))
		if self.pattern == 'Pins': # white pins on black background or vice versa
			spacing = np.int(self.pinSpacing)
			width = np.int(self.pinWidth)
			offset = np.int(self.offset) % (spacing*spacing)
			xOffset = offset % spacing
			yOffset = np.int(np.float(offset)/np.float(spacing))
			mtxPattern = np.zeros([self.imgWidth,self.imgHeight],dtype=np.uint8)
			for i in range(width):
				for j in range(width):
					mtxPattern[i+xOffset::spacing,j+yOffset::spacing] = 255
					#mtxPattern[i+spacing/2::spacing,j+spacing/2::spacing] = 255
			if self.blackBackground == False:
				imgPattern = Image.fromarray(255-mtxPattern, 'L')
				self.imgMerged = Image.merge('RGB',(imgPattern,imgPattern,imgPattern))
			else:
				imgPattern = Image.fromarray(mtxPattern, 'L')
				self.imgMerged = Image.merge('RGB',(imgPattern,imgPattern,imgPattern))
		if self.pattern == 'Stripes':
			spacing = np.int(self.stripeSpacing)
			width = np.int(self.stripeWidth)
			offset = np.int(self.offset) % spacing
			mtxPattern = np.zeros([self.imgWidth,self.imgHeight],dtype=np.uint8)
			for i in range(stripeWidth):
				mtxPattern[:,i+offset::stripeSpacing] = 255
			if self.blackBackground == False:
				imgPattern = Image.fromarray(255-mtxPattern, 'L')
				self.imgMerged = Image.merge('RGB',(imgPattern,imgPattern,imgPattern))
			else:
				imgPattern = Image.fromarray(mtxPattern, 'L')
				self.imgMerged = Image.merge('RGB',(imgPattern,imgPattern,imgPattern))
		if self.pattern == 'RGB': # three colors interleaved with black or white background
			spacing = np.int(self.pinSpacing)
			width = np.int(self.pinWidth)
			if self.blackBackground == True:
				mtxPattern = np.zeros([self.imgWidth,self.imgHeight],dtype=np.uint8)
				mtxPattern2 = np.zeros([self.imgWidth,self.imgHeight],dtype=np.uint8)
				mtxPattern3 = np.zeros([self.imgWidth,self.imgHeight],dtype=np.uint8)
				for i in range(width):
					for j in range(width):
						mtxPattern[i::spacing,j::spacing] = 255
						#mtxPattern[i+spacing/2::spacing,j+spacing/2::spacing] = 255
						mtxPattern2[i+spacing/2::spacing,j+spacing/2::spacing] = 255
						#mtxPattern2[i+spacing/2+spacing/2::spacing,j+spacing/2::spacing] = 255
						mtxPattern3[i+spacing/2::spacing,j::spacing] = 255
						mtxPattern3[i::spacing,j+spacing/2::spacing] = 255
				imgPattern = Image.fromarray(mtxPattern, 'L')
				imgPattern2 = Image.fromarray(mtxPattern2, 'L')
				imgPattern3 = Image.fromarray(mtxPattern3, 'L')
				self.imgMerged = Image.merge('RGB',(imgPattern,imgPattern2,imgPattern3))
			else:
				mtxPattern = np.ones([self.imgWidth,self.imgHeight],dtype=np.uint8)*255
				mtxPattern2 = np.ones([self.imgWidth,self.imgHeight],dtype=np.uint8)*255
				mtxPattern3 = np.ones([self.imgWidth,self.imgHeight],dtype=np.uint8)*255
				for i in range(width):
					for j in range(width):
						mtxPattern[i::spacing,j::spacing] = 0
						mtxPattern2[i::spacing,j::spacing] = 0
						#mtxPattern[i+spacing/2::spacing,j+spacing/2::spacing] = 255
						mtxPattern2[i+spacing/2::spacing,j+spacing/2::spacing] = 0
						mtxPattern3[i+spacing/2::spacing,j+spacing/2::spacing] = 0
						#mtxPattern2[i+spacing/2+spacing/2::spacing,j+spacing/2::spacing] = 255
						mtxPattern[i+spacing/2::spacing,j::spacing] = 0
						mtxPattern[i::spacing,j+spacing/2::spacing] = 0
						mtxPattern3[i+spacing/2::spacing,j::spacing] = 0
						mtxPattern3[i::spacing,j+spacing/2::spacing] = 0
				imgPattern = Image.fromarray(mtxPattern, 'L')
				imgPattern2 = Image.fromarray(mtxPattern2, 'L')
				imgPattern3 = Image.fromarray(mtxPattern3, 'L')
				self.imgMerged = Image.merge('RGB',(imgPattern,imgPattern2,imgPattern3))
		if self.pattern == 'Red':
			self.imgMerged = Image.merge('RGB', 	(self.imgFull,self.imgEmpty,self.imgEmpty))
		if self.pattern == 'Green':
			self.imgMerged = Image.merge('RGB', (self.imgEmpty,self.imgFull,self.imgEmpty))
		if self.pattern == 'Blue':
			self.imgMerged = Image.merge('RGB', (self.imgEmpty,self.imgEmpty,self.imgFull))
		self.updateImg()

	#def changeColor(self, sender):
	#	self.patternColor = sender.title
	#	self.updateImg()

		#			console.hud_alert('Saved')
	#	else:
	#		console.hud_alert('No Image', 'error')
	
	@ui.in_background
	def runserver(self):
		s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
		s.bind(('',50000))
		s.listen(1)
		#s.settimeout(0.2)
		while self.keepServerOn:
			conn, addr = s.accept() # this blocks execution. Need a timeout strategy that doesnt throw an error, so we can keep checking the keepServerOn method
			while 1:
				data = conn.recv(32)
				if not data:
					break
				conn.sendall(data)
				pattern, width, spacing, offset, blackBackground = struct.unpack('ciiii',data)
				if ord(pattern) == ord('p'): # not sure why I have to use ord, which converts to ascii number
					self.pattern = 'Pins'
					self.pinWidth = width
					self.pinSpacing = spacing
				if ord(pattern) == ord('s'):
					self.pattern = 'Stripes'
					self.stripeWidth = width
					self.stripeSpacing = spacing
				if ord(pattern) == ord('r'):
					self.pattern = 'RGB'
					self.pinWidth = width
					self.pinSpacing = spacing
				if ord(pattern) == ord('R'):
					self.pattern = 'Red'
				if ord(pattern) == ord('G'):
					self.pattern = 'Green'
				if ord(pattern) == ord('B'):
					self.pattern = 'Blue'
				if ord(pattern) == ord('A'):
					self.pattern = 'Random'
				if blackBackground == 0:
					self.blackBackground = False
				else:
					self.blackBackground = True
				self.offset = offset
				self.updatePattern()
			conn.close()
		s.close()
			
	def will_close(self):
		self.keepServerOn = False;
	
	
# We use a square canvas, so that the same image
# can be used in portrait and landscape orientation.
w, h = ui.get_screen_size()
canvas_size = max(w, h)

sv = SketchView(canvas_size, canvas_size)
sv.name = 'Zebrography'
#s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#s.bind(('',8080))
#s.listen(5)
#x = 
sv.name = get_local_ip() + ':5000' #socket.gethostbyname(socket.getfqdn() + '.local')
sv.present(style = 'fullscreen', orientations = ['landscape'])
sv.runserver()

