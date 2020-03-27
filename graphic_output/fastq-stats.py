import argparse
from PIL import Image, ImageDraw, ImageFont, ImageOps


parser = argparse.ArgumentParser()
parser.add_argument('-fasq', action="store", dest = 'File', required = "True")
parser.add_argument('-out', action="store", dest = 'out', required = "True")

args = parser.parse_args()
reverse = "no"
files = args.File
out = args.out


#phred_not="J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,[,\,],^,_,`,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,{,|,},~"

def character_to_ASCII(string):
	st = []
	phred_result = 0

	for items in string:
		if len(items)>1:
			for i in items:
				ascii = ord(i)
				if int(ascii) >= 75:
					phred_result = 1
				ascii = ascii-33
				
		else:
			ascii = ord(items)
			ascii = ascii-33
			if ascii >= 75:
				phred_result = 1#if ascii not in -33 result = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		st.append(int(ascii))

	return st,phred_result

def fasq_preprocess(fil):
	fastaq_list = []
	i = 0
	with open(fil) as fastaq:
		p = 1
		for lines in fastaq:
			lines = lines.rstrip()
			if p%4 == 0:
				fastaq_list.append(lines)
			i += 1
			p +=1
			if i > 20000:
				break
		ascii_check = character_to_ASCII(fastaq_list)[1]
		return fastaq_list, ascii_check
def fasq_process(fil,x,group):
	dic = {}
	dic[str(x)+"_"+str(x+group)] = []
	i = 0
	p = 1
	with open(fil) as fastaq:
		for lines in fastaq:
			lines = lines.rstrip()
			if p%4== 0:
				for items in range(x,x+group):
					#print x, x+group,items
					try:
						dic[str(x)+"_"+str(x+group)].append(lines[items])
					except:
						h = "h"
				i += 1
			
			if i >= 5000:
				break
			p += 1
	qual_dic= {}
	for position in dic:
		qual_dic[position] = character_to_ASCII(dic[position])[0]

	return qual_dic
def average(lista):
	n = 0 
	l = []
	for items in lista:
		l.append(float(items))
		n += 1
	try:
		average = sum(l)/n
	except:
		average = 0
	return average 
def boxplot_stats(lista):
	lista = sorted(lista)
	lenght = len(lista)
	position1 = float(lenght/2)
	if position1.is_integer():
		per_50 = lista[int(position1)]
	else:
		position1 = position1.split(".")[0]
		position2 = position1 + 1
		position1 =  lista[int(position1)]
		position2 = lista[int(position2)]
		per_50 = (postion1 + position2)/2 

	position2 = float((lenght+1)/4)
	if position2.is_integer():
		per_25 = lista[int(position2)]
	else:
		position2 = position2.split(".")[0]
		position3 = position2 + 1
		position2 =  lista[int(position2)]
		position3 = lista[int(position3)]
		position4 = position2 + position3*position2.split(".")[1]*((lenght+1)/4)
		per_25 = int(position4) 

	position3 = float((lenght+1)*3/4)
	if position3.is_integer():
		per_75 = lista[int(position3)]
	else:
		position3 = position3.split(".")[0]
		position4 = position4 + 1
		position3 =  lista[int(position3)]
		position4 = lista[int(position4)]
		position5 = position3 + position4*position3.split(".")[1]*((lenght+1)*3/4)
		per_75 = int(position5) 
	
	IQR = int(per_75) - int(per_25)
	Extreme_positive = 1.5*IQR + float(per_75)
	Extreme_negative =  float(per_25)-1.5*IQR 
	
	if Extreme_negative < 0:
		Extreme_negative = 0
	
	if float(min(lista))> float(Extreme_negative): Extreme_negative = min(lista)
	
	if float(max(lista))< float(Extreme_positive): Extreme_positive = max(lista)
	
	
	return float(Extreme_negative), float(per_25), float(per_50), float(per_75), float(Extreme_positive)

def calculations(dic_pos): 	
	for item in dic_pos:
		average_val = average(dic_pos[item])
		stats = boxplot_stats(dic_pos[item])

	return average_val,stats[0],stats[1],stats[2],stats[3],stats[4]

def lenght_reads_cal(pre_dic):
	value_list = []
	biggest = 0
	for reads in pre_dic:	
		if len(reads) > biggest:
			biggest = len(reads) 
	return biggest

def Draw_box_plot(table,out):
	fnt1 = ImageFont.truetype('./fonts/VeraMono.ttf', 14)
	#Size of the window

	a = 60
	b = 70
	c = 20
	x_window = len(table)*10+ 10+ a + b 
	y_window = 460
	#Generation of the file
	im = Image.new("RGB", (x_window, y_window), (255,255,255))	
	draw = ImageDraw.Draw(im)
	#Creation of the axiss: axiss will start with an indentation of 20 above, below and beside. Each Phred quality score will be in a X% proportion of the y_window  px
	size_y_axis = y_window-80 #Total size minus above and below indentations
	position_y_axis= size_y_axis+20
	draw.line(((a, c) + (a, position_y_axis)), fill=(0, 0, 0, 0), width=1)
	size_x_axis = len(table)*10 +10 #number of positions*10 pxls which is one will take + 10 for the position 1. 
	draw.line(((a, position_y_axis) + (a+size_x_axis, position_y_axis)), fill=(0, 0, 0, 0), width=1) 
	#Close chart
	draw.line(((a, c) + (size_x_axis+a, c)), fill=(0, 0, 0, 0), width=1)
	draw.line(((size_x_axis+a, c) + (size_x_axis+a, position_y_axis)), fill=(0, 0, 0, 0), width=1)
	
	#Vertical values
	step = float(size_y_axis)/42
	j = 0
	for values in range(42,-1,-1):
		
		draw.line(((a,20+abs(values-42)*step) + (a-4,20+abs(values-42)*step)), fill=(0, 0, 0, 0), width=1)
		
		if values%5 == 0:
	
			draw.line(((a,20+abs(values-42)*step) + (a-6,20+abs(values-42)*step)), fill=(0, 0, 0, 0), width=1)
			text = abs(values-42)
			w, h = draw.textsize(str(text))
			draw.text((a-25,20+text*step-h/2-3), str(values), font=fnt1, fill=(0,0,0,0))
		j +=1	
			
	i = 10 + a #indentation + space for the first box (same space as in size_x_axis)
	for position in table:
		name = position[0]
		position = position[1:]

		#write the position in the x axis
		draw.line(((i, position_y_axis) + (i, position_y_axis+4)), fill=(0, 0, 0, 0), width=1)

		if (i-a)%50 == 0 :
			draw.line(((i, position_y_axis) + (i, position_y_axis+6)), fill=(0, 0, 0, 0), width=1)
			w, h = draw.textsize(str(name))
			if len(name)==3:draw.text((i-w/2-4+1 , position_y_axis+10), name, font=fnt1, fill=(0,0,0,0))
			if len(name) == 4: draw.text((i-w/2-3 +4, position_y_axis+10), name, font=fnt1, fill=(0,0,0,0))				
			elif len(name) == 5: draw.text((i-w/2-4, position_y_axis+10), name, font=fnt1, fill=(0,0,0,0))
			elif len(name) == 6: draw.text((i-w/2-9, position_y_axis+10), name, font=fnt1, fill=(0,0,0,0))
			elif len(name) == 7: draw.text((i-w/2-14, position_y_axis+10), name, font=fnt1, fill=(0,0,0,0))
		#Create a line from the begining to the end of the parameters
		beg = float(position[1]) * step
		end = float(position[-1]) * step
		draw.line(((i, position_y_axis-beg) + (i, position_y_axis-end)), fill=(0, 0, 0, 0), width=1)
		#Close the whiskers
		draw.line(((i-1, position_y_axis-beg)+(i+1, position_y_axis-beg)),fill=(0, 0, 0, 0), width=1)
		draw.line(((i-1, position_y_axis-end)+(i+1, position_y_axis-end)),fill=(0, 0, 0, 0), width=1)


		#Create the boxplot 
		beg = float(position[2]) * step
		end = float(position[-2]) * step
		draw.rectangle([(i-3, position_y_axis-beg), (i+3, position_y_axis-end)], fill=(24, 56, 214), outline= None)

		#Draw the average and the MEDIANA?
		av = float(position[0]) * step
		med = float(position[3]) * step
		draw.line(((i-3, position_y_axis-med) + (i+3, position_y_axis-med)), fill=(191, 17, 54), width=1)
		draw.line(((i, position_y_axis-av) + (i, position_y_axis-av)), fill=(50, 214, 25), width=1)
		i +=10

	#axiss draw:
	x_name = "Position (bp)"
	draw.text((size_x_axis/2, position_y_axis+35), x_name , font=fnt1, fill=(0,0,0,0)) #Horizontal
	y_name ="Quality (phred)"






	#label = Image.new("RGB", (140, 20), (255,255,255))
	#draw2 = ImageDraw.Draw(label)
	#draw2.text((1, 1), y_name, font=fnt1, fill=(0,0,0))
	#label = label.rotate(90)
	#im.paste(label, (2,size_y_axis/2-50))

	txt=Image.new('L', (500, 50))
	d = ImageDraw.Draw(txt)
	d.text( (0, 0), y_name,  font=fnt1, fill=255)
	w=txt.rotate(90,  expand=1)
	im.paste( ImageOps.colorize(w, (0,0,0), (0,0,0)), (2,-(size_y_axis/2+20)),  w)


	#Vertical
	#save image, specifying the format with the extension
	im.save(out)


pre_dic = fasq_preprocess(files)[0]
phred_result = fasq_preprocess(files)[1]
if phred_result == 1:
	print phred_result
	exit()
lenght_reads = int(lenght_reads_cal(pre_dic))

if lenght_reads<80: group = 1
if 100>=lenght_reads>80: group = 2
if lenght_reads>100: group = 5


x =0
final_dic= {}
#print range((lenght_reads/group)-1)
#exit()
final_list=[]

for repetitions in range((lenght_reads/group)-1):
	lis = []
	position_table=fasq_process(files,x,group)
	result = calculations(position_table)
	lis.append(str(x+1)+"-"+str(x+group))
	lis.extend(result)
	final_list.append(lis)
	
	#final_dic[str(x+1)+"-"+str(x+group)]= result
	x += group

Draw_box_plot(final_list,out)

print phred_result


