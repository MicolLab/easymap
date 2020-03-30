import argparse
from PIL import Image, ImageDraw, ImageFont, ImageOps

#Si no esta en orden: ../samtools1/samtools sort alignment1.sam > alignment1.bam
#para obtener el archivo de profundidad: ../samtools1/samtools depth -a  alignment1.bam > coverage.txt

parser = argparse.ArgumentParser()
parser.add_argument('-coverages', action="store", dest = 'cov', required = "True")
parser.add_argument('-out', action="store", dest = 'out', required = "True")

args = parser.parse_args()


def read_file(f):
	dic = {}
	with open(f,"r") as coverage_file:
		for lines in coverage_file:
			lines = lines.rstrip()
			column = lines.split("\t")
			if int(column[2]) in dic:
				if not int(column[2]) == 0:
					dic[int(column[2])]+=1
			else:
				if not int(column[2]) == 0:
					dic[int(column[2])] = 1


	list_of_values= dic.values()
	total_values= float(sum(list_of_values))
	
	dic_percentage= {}
	for values in dic:
		dic_percentage[values] = (dic[values]/total_values)
	
	sort_position = sorted(dic.keys())
	return dic_percentage, sort_position


def draw(dic,sort_positions,out):
	fnt1 = ImageFont.truetype('./fonts/VeraMono.ttf', 14)
	group = 2
	#Size of the window
	a = 80 #right width space
	c = 80 # left width space
	b= 80 #up sapace
	d = 80 #down space

	size_x_window = 480 
	size_y_window = 480

	x_window = size_y_window+ a +c
	y_window = size_y_window+ b+ d
	
	#Generation of the file
	im = Image.new("RGB", (x_window, y_window), (255,255,255))	
	draw = ImageDraw.Draw(im)

	draw.line(((a, b) + (a, size_y_window+b)), fill=(0, 0, 0, 0), width=1)

	draw.line(((a, size_y_window+b) + (a+size_x_window, a+size_y_window)), fill=(0, 0, 0, 0), width=1) 
	
	#Close chart
	draw.line(((a, c) + (size_x_window+a, c)), fill=(0, 0, 0, 0), width=1)
	draw.line(((size_x_window+a, c) + (size_x_window+a, size_y_window+b)), fill=(0, 0, 0, 0), width=1)



	#Vertical values
	
	step = float(size_y_window)/100
	v = max(dic.values())
	v = v+ 0.1*v
	upper = size_y_window
	step2 = upper/v


	for values in range(100,-1,-1):
		draw.line(((a,b+abs(values-100)*step) + (a-4,b+abs(values-100)*step)), fill=(0, 0, 0, 0), width=1)
		if values%10 == 0:
			draw.line(((a,b+abs(values-100)*step) + (a-6,b+abs(values-100)*step)), fill=(0, 0, 0, 0), width=1) 
			write = str(float(values * v))
			w = write.split(".")
			write = str(w[0])+"."+str(w[1][0])

			draw.text((a-40,b+abs(values-100)*step-5), str(write), font=fnt1, fill=(0,0,0,0))


	#Horizontal values
	step = size_x_window/120
	for values in range(121):
		
		draw.line(((c+values*step, size_y_window+b) + (c+values*step, size_y_window+4+b)), fill=(0, 0, 0, 0), width=1)
		if values%10 == 0:
			draw.line(((c+values*step, size_y_window+b) + (c+values*step, size_y_window+6+b)), fill=(0, 0, 0, 0), width=1)
			draw.text((c+values*step-3, size_y_window+b+10), str(values), font=fnt1, fill=(0,0,0,0))
	order = range(len(sort_positions))

	#Vertical axis name:
	#axis draw:
	x_name = "Read depth (X)"
	draw.text((size_x_window/2+10, size_y_window+b+35), x_name , font=fnt1, fill=(0,0,0,0)) #Horizontal
	y_name ="Frequency (%)"
	#label = Image.new("RGB", (140, 20), (255,255,255))
	#draw2 = ImageDraw.Draw(label)
	#draw2.text((1, 1), y_name, font=fnt1, fill=(0,0,0))
	#label = label.rotate(90)
	#im.paste(label, (2,size_y_window/2+c-80))
	txt=Image.new('L', (500, 50))
	d = ImageDraw.Draw(txt)
	d.text( (0, 0), y_name,  font=fnt1, fill=255)
	w=txt.rotate(90,  expand=1)
	im.paste( ImageOps.colorize(w, (0,0,0), (0,0,0)), (2,-size_y_window/3+20),  w)



	#Graph draw
	#draw.line(((a, size_y_window+b),(c+sort_positions[0]*step,  size_y_window+b)), fill=(0, 0, 0, 0), width=1)
	draw.line(((c+sort_positions[0]*step,  size_y_window+b),(c+sort_positions[0]*step,  size_y_window+b-dic[sort_positions[0]]*step2)), fill=(0, 0, 0, 0), width=1)
	draw.line(((c+sort_positions[-1]*step,  size_y_window+b),(c+sort_positions[-1]*step,  size_y_window+b-dic[sort_positions[-1]]*step2)), fill=(0, 0, 0, 0), width=1)
	for i in order:
		x_coverage= sort_positions[i]
		y_proportion = dic[x_coverage]
##OJOOO
		try:
			if c+sort_positions[i+1]*step > c+size_x_window: #When the line is bigger than the right site of the window
				draw.line(((c+x_coverage*step, size_y_window+b-y_proportion*step2),(c+size_x_window-1, size_y_window+b-dic[sort_positions[i+1]]*step2)))
				break
			#else:
			draw.line(((c+x_coverage*step, size_y_window+b-y_proportion*step2),(c+sort_positions[i+1]*step,  size_y_window+b-dic[sort_positions[i+1]]*step2)), fill=(0, 0, 0, 0), width=1)
		except: #if there is no next i
			draw.line(((c+x_coverage*step, size_y_window+b-y_proportion*step2),(c+x_coverage*step, size_y_window+b-y_proportion*step2)))
# X en try y except tenia -3 despues de step, asi como draw.line de debajo de draw.line
	#im.save(out)
	w, h = im.size
	im.crop((0, 60, w-0, h-0)).save(out)

	# Print out average RD value 
	print max(dic, key=dic.get)




cov = args.cov
file_data = read_file(cov)
draw(file_data[0],file_data[1],args.out)

 
