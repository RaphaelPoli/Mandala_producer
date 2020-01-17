#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import cairo
#import rsvg
import math
import Image

#implement virtual brush as in trace canvas
#implement outline brush
#implement round end brush
#implement disk add v
#implement wave line from a to b
#implement spirals
# implement gone-accolade (the outside of accolades come along straight segments rather than a circle v
#implement a cycle method that would allow to trace a portion of the segments.
#implement a ram version of each figure so that a full tracing procedure can be writen, and groups of dots are correct in svg
#implement a loop removal procedure
#implement a loop outline preparation (adds a point at the crossing and recreates the order)
#implement a hatching (hachure) fill procedure
#implement a filled "pointe" procedure
#implement a two disk and filled pointes to get radial rectangles






## my geometric library ---------------------------------------------------------------------------------------------------------------------------------------

def dot_line_distance(dot,a,b):
	same_point=False
	vertical=False
	horizontal=False
	error=False
	line_equation=""
	line_a=0
	line_b=0
	xa=a[0]
	ya=a[1]
	xb=b[0]
	yb=b[1]
	x=dot[0]
	y=dot[1]
	if (xa==xb and ya==yb) or (x==xa and y==ya) or (x==xb and y==yb):
		same_point=True
		error=True
		distance_to_line=9999
		#print "meme point"
	if ya==yb:
		horizontal=True
	if xa==xb:
		vertical=True
		
	if not (vertical or horizontal or same_point):
		
		line_a=(yb-ya)/float(xb-xa)
		#print ""
		#print "a de dot ",xa,ya," et ",xb,yb," ",line_a
		line_b=(yb-line_a*xb)
		#print "b de dot ",xa,ya," et ",xb,yb," ",line_b
		#print "equation de la droite ",xa,ya,xb,yb,":",line_a,"x+",line_b
		line_equation="y="+str(line_a)+"x+"+str(line_b)
		if line_a==0:
			print "line_a =0"
			horizontal==True
	if not (horizontal or vertical or same_point):
		a_perp=-1/float(line_a)
		b_perp=(y-a_perp*x)
		#print "equation de la perpendiculaire passant par",x,y," y=",a_perp,"x +",b_perp
		intersect_x=((line_a*(xa-ya))-(a_perp*(x+y))) / float(line_a-a_perp)
		#print "x intersect",intersect_x
		intersect_y=a_perp*intersect_x+b_perp
		#print "y intersect",intersect_y
		distance_to_line=math.sqrt(math.fabs(y-intersect_y)**2+math.fabs(x-intersect_x)**2)
		#print "distance to line",distance_to_line
		
	else:
		if horizontal:
			line_equation="hy="+str(xa)
			if y==ya:
				distance_to_line=0
			else:
				distance_to_line=math.fabs(y-ya)
				#print "(",x,y,")-(",xa,ya,")=",distance_to_line,"/","ya",ya,"yb",yb
		if vertical:
			line_equation="vx="+str(xa)
			if x==xa:
				distance_to_line=0
			else:
				distance_to_line=math.fabs(x-xa)
	return [error,distance_to_line,[line_equation,line_a,line_b]]
		#implement here the distance to a vertical or horizontal line
		
		#print"meme point, doite horizontale ou verticale",dot1,dot2
		#print"vertical",vertical
		#print"horizontal",horizontal
		#print"same_point",same_point
		
		
def middle(dot0,dot1):
	return ((dot0[0]+((dot1[0]-dot0[0])/2.0)), (dot0[1]+((dot1[1]-dot0[1])/2.0)))

def is_close(a,b,margin):
	if math.fabs(b-a)<margin:
		return True
	else:
		return False

def is_close_dot(a,b, margin):
	return is_close(a[0],b[0],margin) and is_close(a[1],b[1], margin)
	
def line_intersection(line1, line2):#found on stack exchange
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    
    intersection=(x,y)
    return intersection
    
    
def radius_from_prop(radius, n_sides, prop):#ai je inversé sinus et cosinus?
	angle=math.radians((360.0/float(n_sides)))
	apotheme=radius*math.cos(angle/2.0)#
	cote=(radius*math.sin(angle/2.0))*2.0
	new_radius=apotheme+cote*float(prop)
	return new_radius
	
	
def next_radius(radius, n_sides):#ai je inversé sinus et cosinus?
	angle=math.radians((360.0/float(n_sides)))
	apotheme=radius*math.cos(angle/2.0)#
	cote=(radius*math.sin(angle/2.0))*2.0
	new_radius=apotheme
	return new_radius
	
	
def apotheme(radius, n_sides):#ai je inversé sinus et cosinus?
	angle=math.radians((360.0/float(n_sides)))
	apotheme=radius*math.cos(angle/2.0)#c'est juste
	cote=(radius*math.sin(angle/2.0))*2.0
	return apotheme
	
def radius_from_apotheme(apotheme, n_sides):
	angle=math.radians((360.0/float(n_sides)))
	radius = apotheme/float(math.cos(angle/2.0))
	return radius

def cote(radius, n_sides):#ai je inversé sinus et cosinus?
	angle=math.radians((360.0/float(n_sides)))
	apotheme=radius*math.cos(angle/2.0)#c'est juste
	cote=(radius*math.sin(angle/2.0))*2.0#c'est juste
	return cote
	
def prop_from_radius(radius, n_sides, radius_new):#ai je inversé sinus et cosinus?
	angle=math.radians((360/float(n_sides)))
	apotheme=radius*math.cos(angle/2.0)#
	new_apotheme=new_radius*math.cos(angle/2.0)#
	distance=math.fabs(apotheme-new_apotheme)
	cote=(radius*math.sin(angle/2.0))
	prop=distance/float(cote)
	new_radius=apotheme+cote*prop
	return prop
	
def prop_from_hight_and_radius(demanded_hight,radius):
	prop=demanded_hight/float(radius)
	return prop

def hauteur_from_angle(cote,prop_cote,angle):#marche pour un triangle isocèle
	cote_frag=cote*float(prop_cote)
	angle_r=math.radians(angle/2.0)
	#print "angle",math.degrees(angle_r)
	cote_oppose=cote_frag/2.0
	#print "cote oppose",cote_oppose
	hypoth=float(cote_oppose)/float(math.sin(angle_r))
	#print "hypothénuse",hypoth
	hauteur=math.cos(angle_r)*hypoth
	#print "hauteur", hauteur
	return hauteur
	

def dot_distance(a,b):
	xa=a[0]
	ya=a[1]
	xb=b[0]
	yb=b[1]
	
	distance=math.sqrt(math.fabs((xb-xa)**2+math.fabs((yb-ya)**2)))
	if is_close(xa,xb,0.00001):
		distance=math.fabs(yb-ya)
	#	print "vertical distance"
	if is_close(ya,yb,0.00001):
	#	print "horizontal distance"
		distance=math.fabs(xb-xa)

	return distance

def dot_signed_distance_one_axis(a,b,axis_for_sign):
	sign=1
	xa=a[0]
	ya=a[1]
	xb=b[0]
	yb=b[1]
	distance_in_x=xb-xa
	distance_in_y=yb-ya
	distance=math.sqrt(math.fabs((xb-xa)**2+math.fabs((yb-ya)**2)))
	if xa>xb and axis_for_sign=="x":
		sign=-1
		#print "distance x is negative"
	if ya>yb and axis_for_sign=="y":
		sign=-1
		#print "distance y is negative"
	distance=distance*sign
	return distance

def dot_signed_distance(a,b):
	sign=1
	xa=a[0]
	ya=a[1]
	xb=b[0]
	yb=b[1]
	distance_in_x=xb-xa
	distance_in_y=yb-ya
	distance=math.sqrt(math.fabs((xb-xa)**2+math.fabs((yb-ya)**2)))
	if xa>xb:
		sign=-1
		#print "distance x is negative"
	if ya>yb:
		sign=-1
		#print "distance y is negative"
	distance=distance*sign
	return distance


def from_distance_and_line(a,b,dot, signed_distance):
	line_eq=line_equation(a,b)
	#circle_add(ctx,dot,signed_distance,"ff0000")
	angle=angle_from_line(line_eq,dot)#outputs radians
	
	#print "from_distance",line_eq
	#print "angle:",math.degrees(angle)
	x=(math.cos(math.radians(angle))*float(signed_distance))#difference en x trouvée avec cosinus
	#print "x",x
	x=dot[0]+x
	y=line_eq[1]*x+line_eq[2]
	return (x,y)
	
def angle_from_line(line_eq,dot1):#line eq is a list with string(equation) in 0, a in 1 and b in 2
	x_further=dot1[0]+1
	y_further=line_eq[1]*x_further+line_eq[2]
	angle_r=(math.atan((y_further-dot1[1])/float(x_further-dot1[0])))
	return math.degrees(angle_r)

def angle (dot_a, dot_b, dot_c):#will return the smallest angle of the crossing if in range 90 is demanded
	def back_in_range_90(angle_r):
		while angle_r>math.pi/2.0:
			angle_r=angle_r-math.pi/2.0
		while angle_r<-math.pi/2.0:
			angle_r=angle_r+math.pi/2.0
		if angle_r<0:
			angle_r=math.pi/2.0+angle_r
		return angle_r
	def back_in_range_180(angle_r):
		while angle_r>math.pi:
			angle_r=angle_r-math.pi
		while angle_r<-math.pi:
			angle_r=angle_r+math.pi
		if angle_r<0:
			angle_r=math.pi+angle_r
		return angle_r
	#b is the peak of the angle.
	#print "calculating composing angles"
	line_ab_eq=line_equation(dot_a, dot_b)
	line_bc_eq=line_equation(dot_b, dot_c)
	angle_line_ab=angle_from_line(line_ab_eq, dot_b)
	angle_line_bc=angle_from_line(line_bc_eq, dot_b)
	if is_close(dot_a[0],dot_b[0],0.00001):
		angle_line_ab=math.radians(90)
	if is_close(dot_b[0],dot_c[0],0.00001):
		angle_line_bc=math.radians(90)
	if is_close(dot_a[1],dot_b[1],0.00001):
		angle_line_ab=math.radians(0)
	if is_close(dot_b[1],dot_c[1],0.00001):
		angle_line_bc=math.radians(0)
	#making all angles positive
	#print"making positive"
	if angle_line_ab<0:
		angle_line_ab=angle_line_ab+(math.pi*2)
	if angle_line_bc<0:
		angle_line_bc=angle_line_bc+(math.pi*2)
	#calculating angle (narrowest of crossing)
	result_angle=angle_line_ab-angle_line_bc
	result_angle=back_in_range_180(result_angle)

	return result_angle#in radians
	
def angle_from_line(line_eq,dot1):#line eq is a list with string(equation) in 0, a in 1 and b in 2 
	# dot1 is a dot in the line, the other line for the angle is the horizon
	x_further=dot1[0]-0.2
	
	#print "line_eq angle",line_eq
	y_further=line_eq[1]*x_further+line_eq[2]
	#disk_add(ctx,(x_further,y_further),line_thickness*5,"00ff00")
	
	#lin(ctx,dot1,(x_further,y_further),line_thickness,"ff0000")
	angle_r=(math.atan((y_further-dot1[1])/float(x_further-dot1[0])))
	return angle_r#en radian


def line_equation(a,b):
	same_point=False
	vertical=False
	horizontal=False
	error=False
	line_equation=""
	#print "a:",a," b:",b
	line_a=0
	line_b=0
	xa=a[0]
	ya=a[1]
	xb=b[0]
	yb=b[1]
	if (xa==xb and ya==yb):
		same_point=True
		error=True
		distance_to_line=9999
		return ("same point",0,0)
		#print "meme point"
	if ya==yb:
		horizontal=True
	if xa==xb:
		#print "vertical"
		vertical=True
		
	if not (vertical or horizontal or same_point):
		
		line_a=(yb-ya)/float(xb-xa)
		#print ""
		#print "a de dot ",xa,ya," et ",xb,yb," ",line_a
		line_b=(yb-line_a*xb)
		#print "b de dot ",xa,ya," et ",xb,yb," ",line_b
		#print "equation de la droite ",xa,ya,xb,yb,":",line_a,"x+",line_b
		line_equation="y="+str(line_a)+"x+"+str(line_b)
		if line_a==0:
			#print "line_a =0"
			horizontal==True	
	else:
		if horizontal:
			line_equation="y="+str(ya)
			line_a=0
			line_b=ya
			#print "(",x,y,")-(",xa,ya,")=",distance_to_line,"/","ya",ya,"yb",yb
		if vertical:
			line_equation="x="+str(xa)
			line_a=0
			line_b=0
			
	return [line_equation,line_a,line_b]
	

def gon_dot(dot_center,nsides, radius, dot1_angle, direction, i):
	dot1_angle=math.radians(dot1_angle)
	dot=(dot_center[0] + radius * math.cos((dot1_angle+(2 * math.pi * i / float(nsides)*direction))),
		dot_center[1] + radius * math.sin((dot1_angle+(2 * math.pi * i / float(nsides)*direction))))
	return dot
	
def dot_belongs_to(a,b,dot):
	a_dot=dot_distance(dot,a)
	b_dot=dot_distance(dot,b)
	ab=dot_distance(a,b)
	if a_dot<ab and b_dot<ab:
		return True
	else:
		return False

def dot_is_in(figure,dot):#used to check between two points which one is closer to the barycenter of a figure
	def mean(distance_list):
		result=sum(distance_list)/float(len(distance_list))
		return result
	distance_list=[]
	for line in figure:
		for dot_figure in line:
			distance_list.append(dot_distance(dot_figure,dot))
	return (max(distance_list),mean(distance_list))#mean is currently unused but I was wondering if in some case it could be useful


# geometric motions -------------------------------------
	
def dot_rotate(dot,center,angle_deg):
	angle=(angle_deg/180.0*math.pi)
	x = (((dot[0]-center[0]) * math.cos(angle)) + ((dot[1]-center[1]) * math.sin(angle)))+ center[0]
	y = ((-(dot[0]-center[0]) * math.sin(angle)) + ((dot[1]-center[1]) * math.cos(angle)))+ center[1] 
	new_dot=(x,y)
	return new_dot
	
def translate_absolute(a,b,distance):#here distance is the canvas length, not a percent of a b like in aniation chase
	line_eq=line_equation(a,b)[0]
	line_a=line_equation(a,b)[1]
	line_b=line_equation(a,b)[2]
	#circle_add(ctx,a,distance,"ff0000")
	linear=True
	sign=1
	#if a[0]>b[0] and not unsigned:
	#	sign=-1
	#else:
	#	sign=1
	vertical=False
	horizontal=False
	x=0
	y=0
	#why in adding to x no need to invert and in y yes?
	#print "line equation:",line_eq
	if "x="in line_eq and not "+" in line_eq:#for vertical or horizontal lines do both direction work?
		x=a[0]
		vertical=True
		#print "vertical translate"
		if linear:
			if a[1]>b[1]:				
				y=a[1]-distance*sign
			if b[1]>a[1]:
				y=a[1]+distance*sign
	if "y=" in line_eq and not ("+" in line_eq):
		horizontal=True
		#print "horizontal translate"
		y=a[1]
		if linear:
			if a[0]>b[0]:				
				x=a[0]+distance*sign
			if b[0]>a[0]:
				x=a[0]+distance*sign
	if not(horizontal or vertical):
		#print "slope is true"
		angle_r=angle_from_line(line_equation(a,b),b)
		#print "angle in translate",angle_r
		x_distance=distance*math.cos(angle_r)
		y_distance=distance*math.sin(angle_r)
		#print "distance in translate",x_distance
		if linear:
			x=a[0]+x_distance
			y=line_a*x+line_b
	dot=(x,y)
	#print "dot",dot
	return dot
	
	
def translate_absolute_2(a,b,dot,distance):#here distance is the canvas length, not a percent of a b like in aniation chase
	line_eq=line_equation(a,b)[0]
	line_a=line_equation(a,b)[1]
	line_b=line_equation(a,b)[2]
	#circle_add(ctx,dot,distance,"ff0000")
	linear=True
	
	#sign=1
	#if a[0]>b[0] and not unsigned:
	#	sign=-1
	#else:
	#	sign=1
	vertical=False
	horizontal=False
	x=0
	y=0
	#why in adding to x no need to invert and in y yes?
	#print "line equation:",line_eq
	if is_close(a[0],b[0],0.00001):#for vertical or horizontal lines do both direction work?
		x=dot[0]
		vertical=True
		#print "vertical translate"
		if linear:
			y=dot[1]+distance
	if is_close(a[1],b[1],0.00001):
		horizontal=True
		#print "horizontal translate"
		y=dot[1]
		if linear:
			x=dot[0]+distance
	if not(horizontal or vertical):
		#print "slope is true"
		angle_r=angle_from_line(line_equation(a,b),b)
		#line(line_equation(a,b))
		#print "angle in translate",math.degrees(angle_r)
		x_distance=distance*math.cos(angle_r)*-1
		y_distance=distance*math.sin(angle_r)
		#print "distance x in translate",x_distance
		#print "distance y in translate",y_distance
		if linear:
			x=dot[0]+x_distance
			#y=dot[1]+y_distance
			y=line_a*x+line_b
			#x=dot[0]-x_distance
			#x=(y-line_b)/float(line_a)
	dot=(x,y)
	#print "dot",dot
	return dot
	
	
	
def get_further(a,b,distance):
	line_eq=line_equation(a,b)[0]
	line_a=line_equation(a,b)[1]
	line_b=line_equation(a,b)[2]
	linear=True
	vertical=False
	horizontal=False
	x=0
	y=0
	#print "a",a
	#print "b",b
	#print "line equation:",line_eq
	if "x="in line_eq and not "+" in line_eq:#for vertical or horizontal lines do both direction work?
		x=a[0]
		vertical=True
		#print "vertical"
		if linear:
			if a[1]>b[1]:				
				y=a[1]+distance
			if b[1]>a[1]:
				y=a[1]-distance
		#print "x",x
	#	print "y",y
	if "y=" in line_eq and not ("+" in line_eq):
		horizontal=True
		#print "horizontal"
		y=a[1]
		if linear:
			if a[0]>b[0]:				
				x=a[0]+distance
			if b[0]>a[0]:
				x=a[0]-distance
	if not(horizontal or vertical):
		#print "slope is true"
		if linear:
			#print "line eq",line_eq
			angle_droite=angle_from_line(line_equation(a,b),b)
			if a[0]>b[0]:
				x=dot_rotate((a[0]+distance,a[1]),a,angle_droite)[0]
			if a[0]<b[0]:
				x=dot_rotate((a[0]-distance,a[1]),a,angle_droite)[0]
			y=line_a*x+line_b
	dot=(x,y)
	#print "dot",dot
	return dot
	
def perp_dot(line_a,dot):#does not take vertical or horizontal
	#returns a and b for the perpendicular line that crosses dot (x,y)
	x=dot[0]
	y=dot[1]
	a_perp=-1/float(line_a)
	b_perp=(y-a_perp*x)
	return(a_perp, b_perp)
	
	
def perp_dot2(a,b,dot):# does not allow dot to be a or b (this should be corrected)
	x=dot[0]
	y=dot[1]
	xa=a[0]
	ya=a[1]
	xb=b[0]
	yb=b[1]
	same_point=False
	vertical=False
	horizontal=False
	error=False
	line_equation=""
	line_a=0
	line_b=0
		
	if (xa==xb and ya==yb) or (x==xa and y==ya) or (x==xb and y==yb):
		same_point=True
		error=True
		distance_to_line=9999
		projection_distance=0
		a_perp=0
		b_perp=0
		#print "meme point"
	if is_close(xa,xb,0.00001):
		horizontal=True
	if is_close(ya,yb,0.00001):
		vertical=True
		
	if not (vertical or horizontal or same_point):
		
		line_a=(yb-ya)/float(xb-xa)
		#print ""
		#print "a de dot ",xa,ya," et ",xb,yb," ",line_a
		line_b=(yb-line_a*xb)
		#print "b de dot ",xa,ya," et ",xb,yb," ",line_b
		#print "equation de la droite ",xa,ya,xb,yb,":",line_a,"x+",line_b
		line_equation="y="+str(line_a)+"x+"+str(line_b)
		if line_a==0:
		#	print "line_a =0"
			horizontal==True
	if not (horizontal or vertical or same_point):
		a_perp=-1/float(line_a)
		b_perp=(y-a_perp*x)
		#print "equation de la perpendiculaire passant par",x,y," y=",a_perp,"x +",b_perp
		intersect_x=((line_a*(xa-ya))-(a_perp*(x+y))) / float(line_a-a_perp)
		#print "x intersect",intersect_x
		intersect_y=a_perp*intersect_x+b_perp
		dot_intersect=(intersect_x,intersect_y)
		projection_distance=dot_distance(a,dot_intersect)
		#print "y intersect",intersect_y
		distance_to_line=math.sqrt(math.fabs(y-intersect_y)**2+math.fabs(x-intersect_x)**2)
		#print "distance to line",distance_to_line
		
	else:
		if vertical:# les y sont différents, les x sont les mêmes. (les mots horizontal et vertical sont inversés)
			line_equation="hy="+str(xa)
			projection_distance=dot_distance(a,dot)
			if y==ya:
				distance_to_line=0
				return (999999999999999999,-99999999999999999,vertical,horizontal)
				
			else:#ça n'arrive pas ça si?
				distance_to_line=math.fabs(y-ya)
				return (999999999999999999,-99999999999999999,vertical,horizontal)
				#print "(",x,y,")-(",xa,ya,")=",distance_to_line,"/","ya",ya,"yb",yb
		if horizontal:# les x sont différents, les y sont les mêmes. (les mots horizontal et vertical sont inversés)
			line_equation="vx="+str(xa)
			projection_distance=dot_distance(a,dot)
			if x==xa:
				distance_to_line=0
				a_perp=0
				b_perp=y
			else:#ça n'arrive pas ça si?
				distance_to_line=math.fabs(x-xa)
				a_perp=0
				b_perp=y
	return (a_perp,b_perp,vertical,horizontal)
	
# linear interpolation procedure -------------------------------------------------------------------------------------------------------
# depends on line_equation procedure in my geometric library

	
def retreive_coordinates(curve_array, x):#x est un float ici, alors que array ne contient que des clés integer
	lower_limit=0
	upper_limit=0
	def find_lower_limit_x_to_y(curve_array,xx):
		result=0
		done=False
		for i in range(len(curve_array)+1):
			#print "i",i
			#print "xx",xx
			if i>xx and not done:
				result=i-1
				done=True
				#print "result",result
			else:
				pass
		return result
	lower_limit=find_lower_limit_x_to_y(curve_array,x)
	#print "lower_limit",lower_limit
	upper_limit=find_lower_limit_x_to_y(curve_array,x)+1
	if upper_limit>len(curve_array)-1:
		upper_limit=len(curve_array)-1
		y=curve_array[len(curve_array)-1]
		return y
	#print "upper_limit",upper_limit
	a=(lower_limit,curve_array[lower_limit])
	b=(upper_limit,curve_array[upper_limit])
	line_eq=line_equation(a,b)
	y=line_eq[1]*x+line_eq[2]
	#print "y inside",y
	return y

	
# smart tracing library---------------------------------------------------------------------------------------------------------------------
	
	
def lin (ctx,dot1,dot2,line_thickness,color):
	#print (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2])
	ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
	ctx.move_to(dot1[0], dot1[1])
	ctx.line_to (dot2[0], dot2[1]) # Line to (x,y)
	ctx.set_line_width (line_thickness)
	ctx.stroke ()
		
	ctx.arc(dot1[0], dot1[1],line_thickness/4.0, 0, 2*math.pi)
	ctx.set_line_width (line_thickness/2.0)
	ctx.stroke()
	
	ctx.arc(dot2[0], dot2[1],line_thickness/4.0, 0, 2*math.pi)
	ctx.set_line_width (line_thickness/2.0)
	ctx.stroke()
	
	
def dotted_line(ctx,dot1,dot2, dash_length, rapport, line_thickness, color):
	#dash_length contains empty and full 
	relative=True
	if rapport==0:
		absolute=0.005
		relative=False
	line_eq=line_equation(dot1, dot2)
	line_length=dot_distance(dot1,dot2)
	n_of_dot=line_length/(dash_length)
	#print "n_of_dot",n_of_dot
	
	#lin (dot1,dot2)
	#Vertical lines
	if "x="in line_eq[0]:
		#print"vertical line detected"
		value_of_x=float(line_eq[0].split("=")[1])
		if dot1[1]<dot2[1]:
			for j in range(int(n_of_dot)):
				i=j
				dot1_t=[dot1[0],dot1[1]+i*dash_length]
				if relative:
					dot2_t=[dot1[0],dot1[1]+(i*dash_length)+dash_length*rapport]
				else:
					dot2_t=[dot1[0],dot1[1]+(i*dash_length)+absolute]
				lin (ctx,dot1_t,dot2_t,line_thickness,color)	
		if dot1[1]>dot2[1]:
			for j in range(int(n_of_dot)):
				i=j
				dot1_t=[dot2[0],dot2[1]+i*dash_length]
				if relative:
					dot2_t=[dot2[0],dot2[1]+(i*dash_length)+dash_length*rapport]
				else:
					dot2_t=[dot2[0],dot2[1]+(i*dash_length)+absolute]
				lin (ctx,dot1_t,dot2_t,line_thickness,color)	
	#Horizontal lines
	if line_eq[1]==0:
		#print"horizontal line detected"
		value_of_y=line_eq[2]
		if dot1[0]<dot2[0]:
			for j in range(int(n_of_dot)):
				i=j
				dot1_t=[dot1[0]+i*dash_length,dot1[1]]
				if relative:
					dot2_t=[dot1[0]+(i*dash_length)+dash_length*rapport,dot1[1]]
				else:
					dot2_t=[dot1[0]+(i*dash_length)+absolute,dot1[1]]
				lin (ctx,dot1_t,dot2_t,line_thickness,color)	
		if dot1[0]>dot2[0]:
			for j in range(int(n_of_dot)):
				i=j
				dot1_t=[dot2[0]+i*dash_length,dot2[1]]
				if relative:
					dot2_t=[dot2[0]+(i*dash_length)+dash_length*rapport,dot2[1]]
				else:
					dot2_t=[dot2[0]+(i*dash_length)+absolute,dot2[1]]
				lin (ctx,dot1_t,dot2_t,line_thickness,color)	
	
	if line_eq[1]!=0:
		#print"not vertical or horizontal line detected"
		#put that into a method (get an angle from a line equation)
		x_further=dot1[0]+1
		y_further=line_eq[1]*x_further+line_eq[2]
		angle_r=(math.atan((y_further-dot1[1])/float(x_further-dot1[0])))
		value_of_y=line_eq[2]
		if dot1[0]<dot2[0]:
			for j in range(int(n_of_dot)):
				i=j
				dot1_t=[dot1[0]+((i*dash_length))*math.cos(angle_r),dot1[1]+((i*dash_length))*math.sin(angle_r)]
				if relative:
					dot2_t=[dot1[0]+((i*dash_length)+dash_length*rapport)*math.cos(angle_r),dot1[1]+((i*dash_length)+dash_length*rapport)*math.sin(angle_r)]
				else:
					dot2_t=[dot1[0]+((i*dash_length)+absolute)*math.cos(angle_r),dot1[1]+((i*dash_length)+absolute)*math.sin(angle_r)]
				lin (ctx,dot1_t,dot2_t,line_thickness,color)	
			#to finish the uncomplete dash
			dot1_t=[dot1[0]+(((int(n_of_dot))*dash_length))*math.cos(angle_r),dot1[1]+(((int(n_of_dot))*dash_length))*math.sin(angle_r)]
			dot2_t=[dot1[0]+(((n_of_dot)*dash_length))*math.cos(angle_r),dot1[1]+(((n_of_dot)*dash_length))*math.sin(angle_r)]
			lin (ctx,dot1_t,dot2_t,line_thickness,color)	
		if dot1[0]>dot2[0]:
			for j in range(int(n_of_dot)):
				i=j
				dot1_t=[dot2[0]+((i*dash_length))*math.cos(angle_r),dot2[1]+((i*dash_length))*math.sin(angle_r)]
				if relative:
					dot2_t=[dot2[0]+((i*dash_length)+dash_length*rapport)*math.cos(angle_r),dot2[1]+((i*dash_length)+dash_length*rapport)*math.sin(angle_r)]
				else:
					dot2_t=[dot2[0]+((i*dash_length)+absolute)*math.cos(angle_r),dot2[1]+((i*dash_length)+absolute)*math.sin(angle_r)]
				lin (ctx,dot1_t,dot2_t,line_thickness,color)	
			#to finish the uncomplete dash
			dot1_t=[dot2[0]+(((int(n_of_dot))*dash_length))*math.cos(angle_r),dot2[1]+(((int(n_of_dot))*dash_length))*math.sin(angle_r)]
			dot2_t=[dot2[0]+(((n_of_dot)*dash_length))*math.cos(angle_r),dot2[1]+(((n_of_dot)*dash_length))*math.sin(angle_r)]
			lin (ctx,dot1_t,dot2_t,line_thickness,color)	
	


def line(line_eq):#line eq is a list with in 0 string of line equation, in 1 equation a, in 2 equation b
	#line_eq is output by the procedure line_equation
	
	color=("ffb709")
	dot1=(0,line_eq[2])
	dot2=(1,line_eq[1]+line_eq[2])
	lin (ctx,dot1,dot2,line_thickness,color)
		
	
def trace_figure(line_list,color,close):
	i=0
	for line in line_list:
		i+=1
		#print i
		lin(ctx,line[0],line[1],line_thickness,color)
	if close:
		lin(ctx,line_list[len(line_list)-1][1],line_list[0][0],line_thickness,color)
	
	
def trace_figure_dot(dot_list,color,close):
	for i in range(len(dot_list)-1):
		lin(ctx,dot_list[i],dot_list[i+1],line_thickness,color)
	if close:
		lin(ctx,dot_list[len(dot_list)-1],dot_list[0],line_thickness,color)
	
	
def disk_add(ctx, center, radius, color):#color is hex
	ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
	ctx.arc(center[0], center[1],radius, 0, 2*math.pi)
	ctx.set_line_width (line_thickness*2)
	#ctx.stroke()
	ctx.fill()

def circle_add(ctx, center, radius, color):#color is hex
	ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
	ctx.arc(center[0], center[1],radius, 0, 2*math.pi)
	ctx.set_line_width (line_thickness*2)
	ctx.stroke()
	#ctx.fill()
	
def hexcolor(string):
	red_s=string[0:2]
	green_s=string[2:4]
	blue_s=string[4:6]
	red=int(red_s,16)/255.0
	green=int(green_s,16)/255.0
	blue=int(blue_s,16)/255.0
	return [red,green,blue]
	
	
	
def dotted_line_rnd(ctx,dot1,dot2, dash_length, rapport, line_thickness, color):
	def approx(quantity, variation):
		demi_var=variation/2.0
		variation=-demi_var+random.random()*variation
		output=quantity+variation
		
	first_length=approx(dash_length,variation)
	last_length=approx(dash_length,variation)
	#create prop list
	#create dot list
	#join dots
	
# almost useless ones:

def green_big_circle(dot):
	color="005500"
	ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
	ctx.arc(dot[0], dot[1],line_thickness*2, 0, 2*math.pi)
	ctx.set_line_width (line_thickness*2)
	ctx.stroke()

def purple_circle(dot):
	color="850085"
	ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
	ctx.arc(dot[0], dot[1],line_thickness*2, 0, 2*math.pi)
	ctx.set_line_width (line_thickness*2)
	ctx.stroke()
	
def red_big_circle(dot):
	color="ff0000"
	ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
	ctx.arc(dot[0], dot[1],line_thickness*2, 0, 2*math.pi)
	ctx.set_line_width (line_thickness*2)
	ctx.stroke()
def blue_big_circle(dot):
	color="0000ff"
	ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
	ctx.arc(dot[0], dot[1],line_thickness*2, 0, 2*math.pi)
	ctx.set_line_width (line_thickness*2)
	ctx.stroke()

def red_circle(dot):
	color="ff0000"
	ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
	ctx.arc(dot[0], dot[1],line_thickness, 0, 2*math.pi)
	ctx.set_line_width (line_thickness*2)
	ctx.stroke()
	
def green_circle(dot):
	color="00ff00"
	ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
	ctx.arc(dot[0], dot[1],line_thickness, 0, 2*math.pi)
	ctx.set_line_width (line_thickness*2)
	ctx.stroke()

def red_big_circle(dot):
	color="ff0000"
	ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
	ctx.arc(dot[0], dot[1],line_thickness*2, 0, 2*math.pi)
	ctx.set_line_width (line_thickness*2)
	ctx.stroke()	
	
# library for mandala application ------------------------------------------------------------------------------------------------------------------------------------
	
	
def cycle_dotted(ctx,center, n_sides, radius,dash_length,prop_dash, jump, rotate):
	Dot_list=[]

	for i in range(n_sides):
		Dot_list.append(gon_dot(center, n_sides, radius, rotate, 1,i))


	color="000000"
	for i in range(n_sides):
		j=i+jump
		if j>n_sides-1: 
			j=j-n_sides
		dotted_line(ctx, Dot_list[i], Dot_list[j], dash_length,prop_dash, line_thickness, color)

def cycle(ctx,center, n_sides, radius, jump, rotate,color):
	Dot_list=[]

	for i in range(n_sides):
		Dot_list.append(gon_dot(center, n_sides, radius, rotate, 1,i))


	for i in range(n_sides):
		j=i+jump
		if j>n_sides-1: 
			j=j-n_sides
		lin(ctx, Dot_list[i], Dot_list[j], line_thickness, color)
		

def ram_cycle(center, n_sides, radius, jump, rotate):
	Dot_list=[]

	for i in range(n_sides):
		Dot_list.append(gon_dot(center, n_sides, radius, rotate, 1,i))


	color="000000"
	k=0
	last_k=0
	line_list=[]
	for i in range(n_sides):
		if k<=n_sides-1:
			
			if k>0:
				line_list.append((Dot_list[last_k],Dot_list[k]))
			last_k=k
				
		else:
			g=k-n_sides
			line_list.append((Dot_list[last_k],Dot_list[g]))
			last_k=g
		
		k=k+jump
		j=i+jump
		if j>n_sides-1: 
			j=j-n_sides
		if k>n_sides-1: 
			k=k-n_sides
		#lin(ctx, Dot_list[i], Dot_list[j], line_thickness, color)
		#those lines are not traced in a radial order so cannot be outlined as such
	
	return line_list
	
def pointes_expon(ctx,center, n_sides, radius, width, depth, rotate, close_triangle):
	#width est une fraction du côté
	#depth est une fraction du rayon
	
	Dot_list=[]
	Dot_list0=[]
	Dot_list1=[]
	Dot_list2_=[]
	Dot_list2_0=[]
	Dot_list2_1=[]

	for i in range(n_sides):
		sign=1
		if depth<0:
			sign=-1
		Dot_list.append(gon_dot(center, n_sides, (radius-((radius-apotheme(radius,n_sides))*width))+(radius*depth), 0, 1,i))

	for i in range(n_sides):
		
		
		#the pin is on the gon dot unless rotated
		Dot_list0.append(dot_rotate(gon_dot(center, n_sides, radius, 0, 1,i),center,((360.0/float(n_sides))*width/2.0)))

	for i in range(n_sides):
		Dot_list1.append(dot_rotate(gon_dot(center, n_sides, radius, 0, 1,i),center,-((360.0/float(n_sides))*width/2.0)))
	for i in range(n_sides):
		Dot_list2_1.append(dot_rotate(Dot_list1[i],center,rotate))
		Dot_list2_0.append(dot_rotate(Dot_list0[i],center,rotate))
		Dot_list2_.append(dot_rotate(Dot_list[i],center,rotate))


	color="000000"
	for i in range(n_sides):
		#j=i+jump
		#if j>n_sides-1: 
	#		j=j-n_sides
		#radius=radius_from_side(distance(Dot_list2_0[i],Dot_list2_[i]),n_sides)
		trace_accolade("exponential", Dot_list2_0[i],Dot_list2_1[i],depth,12,1,color)
		if close_triangle:
			lin(ctx, Dot_list2_0[i], Dot_list2_1[i], line_thickness, color)
	
def pointes(ctx,center, n_sides, radius, width, depth, rotate, close_triangle):
	#width est une fraction du côté
	#depth est une fraction du rayon
	
	Dot_list=[]
	Dot_list0=[]
	Dot_list1=[]
	Dot_list2_=[]
	Dot_list2_0=[]
	Dot_list2_1=[]

	for i in range(n_sides):
		sign=1
		if depth<0:
			sign=-1
		Dot_list.append(gon_dot(center, n_sides, (radius-((radius-apotheme(radius,n_sides))*width))+(radius*depth), 0, 1,i))

	for i in range(n_sides):
		
		
		#the pin is on the gon dot unless rotated
		Dot_list0.append(dot_rotate(gon_dot(center, n_sides, radius, 0, 1,i),center,((360.0/float(n_sides))*width/2.0)))

	for i in range(n_sides):
		Dot_list1.append(dot_rotate(gon_dot(center, n_sides, radius, 0, 1,i),center,-((360.0/float(n_sides))*width/2.0)))
	for i in range(n_sides):
		Dot_list2_1.append(dot_rotate(Dot_list1[i],center,rotate))
		Dot_list2_0.append(dot_rotate(Dot_list0[i],center,rotate))
		Dot_list2_.append(dot_rotate(Dot_list[i],center,rotate))


	color="000000"
	for i in range(n_sides):
		#j=i+jump
		#if j>n_sides-1: 
	#		j=j-n_sides
		lin(ctx, Dot_list2_0[i], Dot_list2_[i], line_thickness, color)
		lin(ctx, Dot_list2_[i], Dot_list2_1[i], line_thickness, color)
		if close_triangle:
			lin(ctx, Dot_list2_0[i], Dot_list2_1[i], line_thickness, color)
		
def pointes_ram(center, n_sides, radius, width, depth, rotate, close_triangle):
	#width est une fraction du côté
	#depth est une fraction du rayon
	
	Dot_list=[]
	Dot_list0=[]
	Dot_list1=[]
	Dot_list2_=[]
	Dot_list2_0=[]
	Dot_list2_1=[]

	for i in range(n_sides):
		sign=1
		if depth<0:
			sign=-1
		Dot_list.append(gon_dot(center, n_sides, (radius-((radius-apotheme(radius,n_sides))*width))+(radius*depth), 0, 1,i))

	for i in range(n_sides):
		
		
		#the pin is on the gon dot unless rotated
		Dot_list0.append(dot_rotate(gon_dot(center, n_sides, radius, 0, 1,i),center,((360.0/float(n_sides))*width/2.0)))

	for i in range(n_sides):
		Dot_list1.append(dot_rotate(gon_dot(center, n_sides, radius, 0, 1,i),center,-((360.0/float(n_sides))*width/2.0)))
	for i in range(n_sides):
		Dot_list2_1.append(dot_rotate(Dot_list1[i],center,rotate))
		Dot_list2_0.append(dot_rotate(Dot_list0[i],center,rotate))
		Dot_list2_.append(dot_rotate(Dot_list[i],center,rotate))


	color="000000"
	line_list=[]
	#list_figures=[]
	for i in range(n_sides):
		
		#j=i+jump
		#if j>n_sides-1: 
	#		j=j-n_sides
		line_list.append((Dot_list2_0[i], Dot_list2_[i]))
		line_list.append((Dot_list2_[i], Dot_list2_1[i]))
		#lin(ctx, Dot_list2_0[i], Dot_list2_[i], line_thickness, color)
		#lin(ctx, Dot_list2_[i], Dot_list2_1[i], line_thickness, color)
		if close_triangle:
			#lin(ctx, Dot_list2_0[i], Dot_list2_1[i], line_thickness, color)
			line_list.append((Dot_list2_0[i], Dot_list2_1[i]))
		#list_figures.append(line_list)
	return line_list


def pointes_ram_separate(center, n_sides, radius, width, depth, rotate, close_triangle):
	#width est une fraction du côté
	#depth est une fraction du rayon
	
	Dot_list=[]
	Dot_list0=[]
	Dot_list1=[]
	Dot_list2_=[]
	Dot_list2_0=[]
	Dot_list2_1=[]

	for i in range(n_sides):
		sign=1
		if depth<0:
			sign=-1
		Dot_list.append(gon_dot(center, n_sides, (radius-((radius-apotheme(radius,n_sides))*width))+(radius*depth), 0, 1,i))

	for i in range(n_sides):
		
		
		#the pin is on the gon dot unless rotated
		Dot_list0.append(dot_rotate(gon_dot(center, n_sides, radius, 0, 1,i),center,((360.0/float(n_sides))*width/2.0)))

	for i in range(n_sides):
		Dot_list1.append(dot_rotate(gon_dot(center, n_sides, radius, 0, 1,i),center,-((360.0/float(n_sides))*width/2.0)))
	for i in range(n_sides):
		Dot_list2_1.append(dot_rotate(Dot_list1[i],center,rotate))
		Dot_list2_0.append(dot_rotate(Dot_list0[i],center,rotate))
		Dot_list2_.append(dot_rotate(Dot_list[i],center,rotate))


	color="000000"
	line_list=[]
	list_figures=[]
	for i in range(n_sides):
		line_list=[]
		#j=i+jump
		#if j>n_sides-1: 
		#	j=j-n_sides
		line_list.append((Dot_list2_0[i], Dot_list2_[i]))
		line_list.append((Dot_list2_[i], Dot_list2_1[i]))
		#lin(ctx, Dot_list2_0[i], Dot_list2_[i], line_thickness, color)
		#lin(ctx, Dot_list2_[i], Dot_list2_1[i], line_thickness, color)
		if close_triangle:
			#lin(ctx, Dot_list2_0[i], Dot_list2_1[i], line_thickness, color)
			line_list.append((Dot_list2_0[i], Dot_list2_1[i]))
		list_figures.append(line_list)
	return list_figures
	
def outline_ram(figure,outline_dist,close):#outlining one figure (closes the dot list if close is True)
	#there is still a problem with wide angles: distance is wrong (calculated with cosinus)
	#I tried four angle modification without success (the four folowing methods)
	#should also trunc looping segments
	global line_thickness
	def complementary(angle_r):
		#print "angle_r",math.degrees(angle_r)
		#print "angle complementary", math.degrees(math.pi-angle_r)
		return math.pi-angle_r
	def back_in_range_2(angle_r):
		while angle_r>math.pi/64.0:
			angle_r=angle_r-math.pi/64.0
		while angle_r<-math.pi/64.0:
			angle_r=angle_r+math.pi/64.0
		if angle_r<0:
			angle_r=math.pi/64.0+angle_r
		return angle_r
	def back_in_range_11(angle_r):
		while angle_r>math.pi/16.0:
			angle_r=angle_r-math.pi/16.0
		while angle_r<-math.pi/16.0:
			angle_r=angle_r+math.pi/16.0
		if angle_r<0:
			angle_r=math.pi/16.0+angle_r
		return angle_r
	def back_in_range_45(angle_r):
		while angle_r>math.pi/4.0:
			angle_r=angle_r-math.pi/4.0
		while angle_r<-math.pi/4.0:
			angle_r=angle_r+math.pi/4.0
		if angle_r<0:
			angle_r=math.pi/4.0+angle_r
		return angle_r
	def dot_bissect(a,b,c,outline_dist):
		flat=False
		#print
		#print "a:",a
		#print "b:",b
		#print "c:",c
		reject=False
		if a==b or a==c or b==a or b==c:
			reject=True
			return None
		angle_fl=0.0
		angle_fl=angle(a,b,c)
		#print "angle",math.degrees(angle_fl)
		#disk_add(ctx,b,line_thickness*2,"ff0009")
		alignement=dot_line_distance(a,b,c)
		#print "alignement",alignement
		if alignement<0.01 or angle_fl==0:
			#print "detected flat angle"
			p_only=perp_dot2(a,c,b)
			intersect_point=(b[0]+0.1,p_only[0]*(b[0]+0.1)+p_only[1])
			#lin(ctx,b,intersect_point, line_thickness*1,"999900")
			flat=True
		if outline_dist<0:
			demanded_sign=-1
		else:
			demanded_sign=1
		#creating isocel triangle
		if dot_distance(a,b)<dot_distance(b,c):#choosing hald of smaller side as distance
			distance=dot_distance(a,b)/2.0
		else:
			distance=dot_distance(b,c)/2.0
			
		perp_place1=translate_absolute_2(a,b,b,distance)
	
		#inverting perp place if it is outside figure segment
		if dot_belongs_to (a,b,perp_place1):
			pass
		else:
			perp_place1=translate_absolute_2(a,b,b,-distance)
		perp_place2=translate_absolute_2(b,c,b,distance)
		if dot_belongs_to (b,c,perp_place2):
			pass
		else:
			pass
			perp_place2=translate_absolute_2(b,c,b,-distance)
		
	#	disk_add(ctx,perp_place1,line_thickness*5,"ff0009")
	#	disk_add(ctx,perp_place2,line_thickness*5,"ff0009")
		
		
		#creating bissecrice
		p1=perp_dot2(a,b,perp_place1)
		p2=perp_dot2(b,c,perp_place2)
		p1_dot2=(perp_place1[0]+0.1,p1[0]*(perp_place1[0]+0.1)+p1[1])
		p2_dot2=(perp_place2[0]+0.1,p2[0]*(perp_place2[0]+0.1)+p2[1])
		#lin(ctx,perp_place1,p1_dot2, line_thickness*1,"999900")
		if not (flat or reject):
			intersect_point=line_intersection((perp_place1,p1_dot2),(perp_place2,p2_dot2))
		
		#print "intersect",intersect_point
		#line(line_equation(b,intersect_point))
		
		#this eroneous calculation was used for beautiful perspective mandalas.
		#if angle_fl>math.pi/2.0:
		#	angle_new=complementary(angle_fl)
		#	angle_new=back_in_range_11(angle_new)
		#angle_new=back_in_range_45(angle_fl)
		#hyp=math.cos(angle_new)*outline_dist
		
		angle_new=angle_fl/2.0
		hyp=math.sin(angle_new)*outline_dist
		dot1=translate_absolute_2(b,intersect_point,b,hyp)
		dot2=translate_absolute_2(b,intersect_point,b,-hyp)
		#disk_add(ctx,dot1,line_thickness*5,"00b709")
		#disk_add(ctx,dot2,line_thickness*5,"00b709")
			
		#if dot_belongs_to(intersect_point,b,dot):# this was to know is dot is inside the figure but does not work
		#	if demanded_sign==-1:
		#		pass
		#		print "place of dot is correct"
		#	else:
		#		print "inverting dot"
		#		dot=translate_absolute_2(b,intersect_point,b,-hyp)
		#if not dot_belongs_to(intersect_point,b,dot):# this was to know is dot is outside the figure but does not work
		#	if demanded_sign==1:
		#		pass
		#		print "place of dot is correct"
		#	else:
		#		print "inverting dot"
		#		dot=translate_absolute_2(b,intersect_point,b,-hyp)
		
		if dot_is_in(figure,dot1)[0]<=dot_is_in(figure,dot2)[0]:# if this is true intersect should be inside of the figure
			if demanded_sign==-1:
				dot=dot1
				#print "default place of dot is correct"
			else:
				#print "inverting dot"
				dot=dot2	
		if dot_is_in(figure,dot1)[0]>=dot_is_in(figure,dot2)[0]:# if this is true intersect should be outside of the figure
			if demanded_sign==1:
				dot=dot1
				#print "default place of dot is correct"
			else:
				#print "inverting dot"
				dot=dot2	
	#if all dots are equidistant from the bissect dot then the outline will be made in positive sign distance
	
		#for sign of hyp check if dot is between peak and intersect
		#disk_add(ctx,dot,line_thickness*3,"ffb709")
		return dot
	#print
	#print "input_figure:",figure
	dot_list=[]
	#lin(ctx,figure[1][0],figure[1][1],line_thickness,"0000ff")
	#lin(ctx,figure[2][0],figure[2][1],line_thickness,"ffff00")
	
	#first angle (beware if first and last points are the same!)
	if close:
		dot=dot_bissect(figure[len(figure)-1][1],figure[0][0],figure[0][1], outline_dist)
		#disk_add(ctx,figure[0][0],line_thickness*10,"ffb709")
		if dot!=None:
			dot_list.append(dot)
	#line(("",p1[0],p1[1]))
	#line(("",p2[0],p2[1]))
	#print "dist demand", distance, "calc:",dot_distance(b,perp_place2)
	if close:
		count=len(figure)-1
	else:
		count=len(figure)-1
	for i in range(count):#cycling thru lines (minus closing angles (start and end))
		#print "angle number",i
		#if i ==59:
		#	disk_add(ctx,figure[i][1],line_thickness*4,"ffb7ff")
		dot=dot_bissect(figure[i][0],figure[i][1],figure[i+1][1], outline_dist)
		if dot!=None:
			dot_list.append(dot)
		#print "angle",i,":",math.degrees(angle_i)
	#final dot in closed figures,
	if not close:
		dot=dot_bissect(figure[len(figure)-1][0],figure[0][0],figure[0][1],outline_dist)
		if dot!=None:
			dot_list.append(dot)
	if close:
		dot=dot_bissect(figure[len(figure)-1][0],figure[len(figure)-1][1],figure[0][0], outline_dist)
		#disk_add(ctx,figure[0][0],line_thickness*6,"10ff09")
		if dot!=None:
			dot_list.append(dot)
	line_list=[]
	for i in (range(int(len(dot_list)-1))):
		line_list.append((dot_list[int(i)],dot_list[int(i)+1]))
	#print "last angle ",i,":",math.degrees(last_angle)
	return line_list
	


def cycle_cycle(ctx, center,radius_main, radius_secondary, n_side_main, n_side_secondary, jump_secondary, rotate_main, rotate_secondary,color):
	local_angle=0
	Dot_list=[]
	for i in range(n_side_main):
		Dot_list.append(gon_dot(center, n_side_main, radius_main, rotate_main, 1,i))
		
	for dot in Dot_list:
		local_angle=local_angle+360.0/float(n_side_main)
		cycle(ctx,dot, n_side_secondary, radius_secondary, jump_secondary, rotate_secondary+local_angle,color)


def ram_cycle_accolade(center,radius, curve, n_of_sides, resolution, ratio_width,depth, direction, rotate):
	#this procedure traces a circle of accolade and outputs each line in a line list
	#please note that lines are not in a radial order, so a special outline procedure is needed to outline such a list
	Dots2=[]
	accolade_list=[]
	Dots=[]
	Dot_list=[]
	Dot_list1=[]
	Dot_list2_=[]
	Dot_list2_0=[]
	Dot_list2_1=[]
	
	for i in range(n_of_sides):
		Dot_list.append(gon_dot(center, n_of_sides, radius, 0, 1,i))

	for i in range(n_of_sides):
	
		if i<n_of_sides-1:
			dot=middle(Dot_list[i],Dot_list[i+1])
			if dot[1]<origin[1]:
				inv=-1#used to invert the upper half
			else:
				inv=1
			angle=math.degrees(angle_from_line(line_equation(origin,dot),origin))
		

			Dots=accolade(curve,Dot_list[i],Dot_list[i+1], depth,radius,resolution,direction)
		else:
			dot=middle(Dot_list[i],Dot_list[0])
			if dot[1]<origin[1]:
				inv=-1#used to invert the upper half
			else:
				inv=1
			angle=math.degrees(angle_from_line(line_equation(dot,origin),origin))		
			

			Dots=accolade(curve,Dot_list[i],Dot_list[0], depth,radius,resolution,direction)
		#rotating dots
		#print "dots",Dots
		Dots2=[]
		#print "len dots",len(Dots)
		for i in range(len(Dots)):#cycling thru strokes
			Dots2.append([])
			#print "len stroke",i,":",len(Dots[i])
			for j in range(len(Dots[i])):#cycling thru dots
			#	print "Dot_extracted",Dots[i][j]
				Dots2[i].append(dot_rotate(Dots[i][j],center,rotate))#rotating
			#	print "new_dot",dot_rotate(Dots[i][j],center,rotate)
			#print "2",len(Dots2[i]),"1", len(Dots[i])
		#the pin is on the middle of the segment
		#adding accolade to side list
		accolade_list.append(Dots2)
		#print "Dots",Dots
		#print"len dots2",len(Dots2)
		#print "Dots2",Dots2
	#print "accolade_list",accolade_list

	colors=["000000","ff0000","00ff00","0000ff","ffff00","00ffff","090909"]
	line_list=[]
	for side in range (n_of_sides):# cycing thru sides
		#side=
		#print "side",side
		for dot_index in range(len(accolade_list[side][0])-1):
			#for j in range (resolution-1):
		#	print accolade_list[side][0][dot_index]#cycling thru steps
			#color=colors[0]
			#lin(ctx, accolade_list[side][0][dot_index], accolade_list[side][0][dot_index+1], line_thickness, color)
			line_list.append((accolade_list[side][0][dot_index], accolade_list[side][0][dot_index+1]))
	
		dot_index=0
		for dot_index in range(len(accolade_list[side][1])-1):
			#print"otherside",dot_index
			#print accolade_list[side][1][dot_index]#cycling thru steps
		#for j in range (resolution-1):
			#color=colors[0]
		#	red_big_circle(accolade_list[side][1][dot_index])
			#lin(ctx, accolade_list[side][1][dot_index], accolade_list[side][1][dot_index+1], line_thickness, color)
			line_list.append((accolade_list[side][1][dot_index], accolade_list[side][1][dot_index+1]))
	
		for dot_index in range (2):
			pass
			#color="0000ff"
			#lin(ctx, accolade_list[side][2][dot_index], accolade_list[side][2][dot_index+1], line_thickness, color)
	return line_list
	
	
def cycle_accolade(ctx, center,radius, curve, n_of_sides, resolution, ratio_width,depth, direction, rotate):
	Dots2=[]
	accolade_list=[]
	Dots=[]
	Dot_list=[]
	Dot_list1=[]
	Dot_list2_=[]
	Dot_list2_0=[]
	Dot_list2_1=[]
	
	for i in range(n_of_sides):
		Dot_list.append(gon_dot(center, n_of_sides, radius, 0, 1,i))

	for i in range(n_of_sides):
	
		if i<n_of_sides-1:
			dot=middle(Dot_list[i],Dot_list[i+1])
			if dot[1]<origin[1]:
				inv=-1#used to invert the upper half
			else:
				inv=1
			angle=math.degrees(angle_from_line(line_equation(origin,dot),origin))
		

			Dots=accolade(curve,Dot_list[i],Dot_list[i+1], depth,radius,resolution,direction)
		else:
			dot=middle(Dot_list[i],Dot_list[0])
			if dot[1]<origin[1]:
				inv=-1#used to invert the upper half
			else:
				inv=1
			angle=math.degrees(angle_from_line(line_equation(dot,origin),origin))		
			

			Dots=accolade(curve,Dot_list[i],Dot_list[0], depth,radius,resolution,direction)
		#rotating dots
		#print "dots",Dots
		Dots2=[]
		#print "len dots",len(Dots)
		for i in range(len(Dots)):#cycling thru strokes
			Dots2.append([])
			#print "len stroke",i,":",len(Dots[i])
			for j in range(len(Dots[i])):#cycling thru dots
			#	print "Dot_extracted",Dots[i][j]
				Dots2[i].append(dot_rotate(Dots[i][j],center,rotate))#rotating
			#	print "new_dot",dot_rotate(Dots[i][j],center,rotate)
			#print "2",len(Dots2[i]),"1", len(Dots[i])
		#the pin is on the middle of the segment
		#adding accolade to side list
		accolade_list.append(Dots2)
		#print "Dots",Dots
		#print"len dots2",len(Dots2)
		#print "Dots2",Dots2
	#print "accolade_list",accolade_list

	colors=["000000","ff0000","00ff00","0000ff","ffff00","00ffff","090909"]
#	line_list=[]
	for side in range (n_of_sides):# cycing thru sides
		#side=
		#print "side",side
		for dot_index in range(len(accolade_list[side][0])-1):
			#for j in range (resolution-1):
		#	print accolade_list[side][0][dot_index]#cycling thru steps
			color=colors[0]
			lin(ctx, accolade_list[side][0][dot_index], accolade_list[side][0][dot_index+1], line_thickness, color)
		#	line_list.append((accolade_list[side][0][dot_index], accolade_list[side][0][dot_index+1]))
	
		dot_index=0
		for dot_index in range(len(accolade_list[side][1])-1):
			#print"otherside",dot_index
			#print accolade_list[side][1][dot_index]#cycling thru steps
		#for j in range (resolution-1):
			color=colors[0]
		#	red_big_circle(accolade_list[side][1][dot_index])
			lin(ctx, accolade_list[side][1][dot_index], accolade_list[side][1][dot_index+1], line_thickness, color)
			#line_list.append((accolade_list[side][1][dot_index], accolade_list[side][1][dot_index+1]))
	
		for dot_index in range (2):
			color="0000ff"
			#lin(ctx, accolade_list[side][2][dot_index], accolade_list[side][2][dot_index+1], line_thickness, color)



def hauteur_pointe(radius,n_sides,prop_width, prop_depth):
	#base=cote(radius, n_sides)*prop_width
	hauteur=radius*prop_depth
	return hauteur
	
def losanges(ctx,center,radius, n_sides, prop_width, prop_depth, prop_middle, rotate):
	
	Hauteur_1=(radius*prop_depth)*prop_middle
	prop_depth_1=prop_depth*prop_middle
	pointes(ctx, center,n_sides,radius+Hauteur_1, prop_width,-prop_depth_1, rotate, False)
	pointes(ctx, center,n_sides,radius+Hauteur_1, prop_width,prop_depth*(1-prop_middle), rotate, False)
	

def ram_losanges(center,radius, n_sides, prop_width, prop_depth, prop_middle, rotate):
	new_figure_list=[]
	new_figure_list_2=[]
	Hauteur_1=(radius*prop_depth)*prop_middle
	prop_depth_1=prop_depth*prop_middle
	#creating losanges halfs with pointes
	part_1=pointes_ram_separate(center,n_sides,radius+Hauteur_1, prop_width,-prop_depth_1, rotate, False)
	part_2=pointes_ram_separate(center,n_sides,radius+Hauteur_1, prop_width,prop_depth*(1-prop_middle), rotate, False)
	#joining into one figure list
	for i in range(len(part_1)):
		new_figure=[]
		for line in part_1[i]:
			new_figure.append(line)
		for line in part_2[i]:
			new_figure.append(line)
		new_figure_list.append(new_figure)
	
	# now rearanging each losange
	for i in range(len(new_figure_list)):
		new_figure=[]
		new_figure=re_order(new_figure_list[i])
		new_figure_list_2.append(new_figure)
	return new_figure_list_2
		
	
			
	
	
	
def trace_accolade(curve, dot0,dot1,depth, steps,direction,color):#direction is minus one or one, depth is a distance
	#the origin of the graph is taken as a landmark for vector direction
	
	dot_list=[]
	global origin
	y_end_choose=((0,0),(0,0))
	# radial method has been temporarily abandonned
	
	dictionary_method={"petal":"parallels",
				"exponential":"parallels",
				"reverse exponential":"parallels"}
	if curve=="petal":
		curve_array=petal
	if curve=="exponential":
		curve_array=exponential
	if curve=="reverse exponential":
		curve_array=reverse_exponential
		
	def trace_curve(array,steps,x_max_array,dot_origin,dot_x_end,dot_y_end):# cette pocédure trace la courbe selon un repère donné c'est la procédure la plus élémentale
		sign_y=1
		#purple_circle(dot_origin)
		#print
		if is_close(dot_origin[0],dot_y_end[0], 0.00001):
		
			#print "vertical side detected"
			
			if dot_y_end[1]<=dot_origin[1]:# this if is maybe useless
				#print "side toward up"
				#print "introducing sign - in trace"
				sign_y=-1*sign_y
			else:
				#print "side toward down"
				#print "introducing sign - in trace"
				sign_y=1*sign_y
		
		if is_close(dot_origin[1],dot_y_end[1], 0.00001):
			#print
			#print "horizontal side detected"
			
			if dot_y_end[1]<=dot_origin[1]:# this if is maybe useless
			#	print "side toward right"
			#	print "introducing sign - in trace"
				sign_y=-1*sign_y
			else:
			#	print "side toward left"
			#	print "introducing sign - in trace"
				sign_y=1*sign_y
	
	#	print "sign",sign_y
		
		y_max_array=max(array)
		trace=[]
		full_x_length=dot_signed_distance_one_axis(dot_origin, dot_x_end,"x")
		#print "full_x_length",full_x_length
		if full_x_length<0:
			sign_x=-1
		else:
			sign_x=1
		full_y_length=dot_signed_distance_one_axis(dot_origin, dot_y_end,"y")
		#print "full_y_length",full_y_length
	
		inc=full_x_length/float(steps)
		#print "inc",inc
		inc_array=x_max_array/float(steps)
		#purple_circle(dot_origin)
		for i in range(steps+1):
			dot_a=translate_absolute(dot_origin,dot_x_end, i*inc)#a is the dot on axe x
			
			
				
			
			distance_x=inc_array*i
			
			y=retreive_coordinates(array,distance_x)
			y_ratio=(y*(full_y_length/float(y_max_array)))*sign_y
			

			dot_b=translate_absolute(dot_origin,dot_y_end, y_ratio)
			if not is_close(dot_origin[1],dot_y_end[1], 0.00001) and not is_close(dot_origin[0],dot_y_end[0],0.00001):
					a_coef_axis_y=line_equation(dot_origin,dot_y_end)[1]
					b2=dot_a[1]-a_coef_axis_y*dot_a[0]
					#print a_coef_axis_y
					dot_a_prime=(dot_a[0]+0.2,a_coef_axis_y*(dot_a[0]+0.2)+b2)
					#green_circle(dot_a_prime)
			else:
				if is_close(dot_origin[1],dot_y_end[1], 0.00001):
					if dot_origin[0]<dot_x_end[0]:
						dot_a_prime=(dot_a[0]-0.2,dot_a[1])
					else:
						dot_a_prime=(dot_a[0]+0.2,dot_a[1])
					
				if is_close(dot_origin[0],dot_y_end[0],0.00001):
					if dot_origin[1]<dot_y_end[1]:
						dot_a_prime=(dot_a[0],dot_a[1]-0.2)
					else:
						dot_a_prime=(dot_a[0],dot_a[1]+0.2)
			#print "y",y
			dot_found=translate_absolute(dot_a,dot_a_prime, y_ratio)
			#print "adding dot---------"
			trace.append(dot_found)
		return trace# outputs one list of dots
		
		
		
		
	if steps%2!=0: steps=steps-1
	exponential_segment=(middle(dot0,dot1),middle(middle(dot0,dot1),dot1))
	line_a=line_equation(dot0,dot1)[1]
	perp=perp_dot2(dot0,dot1, exponential_segment[0])#extracting line equation of perpendicular crossing the center of the segment
	distance=math.fabs(dot_signed_distance(dot0,dot1))*depth
	#print perp
	#print "ACCOLDE"
	
	# ici on trouve les points du repère de destination de la courbe
	
	
	#radial y extremity
	
	if perp[2]:
		#print "perp horizontal"
		if dot1[0]<dot0[0]:
			sign=-1
		else:
			sign=1
		
		sign=sign*direction
		dot_high=(exponential_segment[0][0],exponential_segment[0][1]+distance*sign)
		purple_circle(dot_high)
	
		y_end1=(dot1[0],dot1[1]+distance*sign)
		y_end2=(dot0[0],dot0[1]+distance*sign)
		
	if perp[3]:
		#print "perp vertical "
		if dot1[1]<dot0[1]:
			sign=-1
		else:
			sign=1
			
			
		sign=sign*direction
		dot_high=(exponential_segment[0][0]+distance*sign,exponential_segment[0][1])
		y_end1=(dot1[0]+distance*sign,dot1[1])
		y_end2=(dot0[0]+distance*sign,dot0[1])
	
	if not perp[2] and not perp[3]:#lorsque la ligne est en biais
		
		if exponential_segment[0][0]<origin[0]:
			sign=-1
		else:
			sign=1
		
		
		sign=sign*direction
		
		x=exponential_segment[0][0]+0.1#arbitrary new x for prolonging perpendicular
		perp_other_dot=(x,perp[0]*x+perp[1])
		dot_high=translate_absolute(exponential_segment[0], perp_other_dot, distance*sign)
		b2=dot1[1]-perp[0]*dot1[0]
		b3=dot0[1]-perp[0]*dot0[0]
		y_end_step1=(dot1[0]+0.2,(perp[0]*(dot1[0]+0.2))+b2)
		y_end_step2=(dot0[0]+0.2,(perp[0]*(dot0[0]+0.2))+b3)
	
		y_end1=translate_absolute(dot1,y_end_step1,distance*sign)
		y_end2=translate_absolute(dot0,y_end_step2,distance*sign)
	
	#print curve
	#print dictionary_method[curve]
	if dictionary_method[curve]=="parallels":
		#print "parallels method chosen"
		y_end_choose=(y_end1,y_end2)
	if dictionary_method[curve]=="radial":
		y_end_choose=(y_end_1_radial,y_end_2_radial)
		#print "radial method chosen"
	
	Dots_1=trace_curve(curve_array, steps, len(curve_array)-1,dot_high,middle(dot0,dot1),y_end_choose[0])
	Dots_2=trace_curve(curve_array, steps, len(curve_array)-1,dot_high,middle(dot0,dot1),y_end_choose[1])
	
	
	#print "dots accolade",Dots_1
	for dot_index in range(len(Dots_1)-1):
		#print"tracing-----------------------------"
		lin(ctx, Dots_1[dot_index], Dots_1[dot_index+1], line_thickness, color)
	dot_index=0
	for dot_index in range(len(Dots_2)-1):
		lin(ctx, Dots_2[dot_index], Dots_2[dot_index+1], line_thickness, color)



def exponential_concentric(ctx, center,start_radius,fraction_of_radius, segment_de_reference, base_radius_outline, iterations, jump,rotate):
# start_radius is the outside radius
#segment of reference is the width of the whole rythm it's expressed in a fraction of the start radius
#base_radius_outline is an increment of start_radius expressed in the same fraction (can be practical to change the radius with a convenient unit (typically this value may be 0)
#iterations is the number of lines
#jump is what dots are joined by lines
#rotate is a global rotate to align peaks

	span=len(exponential)/float(iterations+2)
	y_array=[]
	
	
	
	for i in range(iterations):
		offset=(i+1)*span
		if offset>len(exponential)-1:
			offset=len(exponential)-1
		y_array.append(retreive_coordinates(exponential,offset))
	def h(y,maxi,apply_to):
		h=y*apply_to/maxi
		return h
	
	hauteur=[]
	for i in range(iterations):
		hauteur.append(segment_de_reference-h(y_array[i],max(exponential),segment_de_reference))

	cycle(ctx, center, 124, start_radius+base_radius_outline*start_radius*fraction_of_radius,jump,rotate,"000000")#base_cycle
	
	
	for i in range(iterations):
		cycle(ctx, center, 124, start_radius+(base_radius_outline+hauteur[i])*start_radius*fraction_of_radius,jump,rotate)


	
	
def accolade(curve, dot0,dot1,depth, radius, steps,direction):#direction is minus one or one, depth is between zero and one it is a fraction of segment length
	#radius is used for radial method
	
	dot_list=[]

	y_end_choose=((0,0),(0,0))
	
	dictionary_method={"petal":"parallels",
				"exponential":"parallels",
				"reverse exponential":"parallels"}
	
	if curve=="petal":
		curve_array=petal
	if curve=="exponential":
		curve_array=exponential
	if curve=="reverse exponential":
		curve_array=reverse_exponential
		
	global origin
	
	def green_big_circle(dot):
		color="005500"
		ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
		ctx.arc(dot[0], dot[1],line_thickness*2, 0, 2*math.pi)
		ctx.set_line_width (line_thickness*2)
		ctx.stroke()
	
	def purple_circle(dot):
		color="850085"
		ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
		ctx.arc(dot[0], dot[1],line_thickness*2, 0, 2*math.pi)
		ctx.set_line_width (line_thickness*2)
		ctx.stroke()
		
	def red_big_circle(dot):
		color="ff0000"
		ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
		ctx.arc(dot[0], dot[1],line_thickness*2, 0, 2*math.pi)
		ctx.set_line_width (line_thickness*2)
		ctx.stroke()
	def blue_big_circle(dot):
		color="0000ff"
		ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
		ctx.arc(dot[0], dot[1],line_thickness*2, 0, 2*math.pi)
		ctx.set_line_width (line_thickness*2)
		ctx.stroke()
	
	def red_circle(dot):
		color="ff0000"
		ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
		ctx.arc(dot[0], dot[1],line_thickness, 0, 2*math.pi)
		ctx.set_line_width (line_thickness*2)
		ctx.stroke()
		
	def green_circle(dot):
		color="00ff00"
		ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
		ctx.arc(dot[0], dot[1],line_thickness, 0, 2*math.pi)
		ctx.set_line_width (line_thickness*2)
		ctx.stroke()
		
	def trace_curve(array,steps,x_max_array,dot_origin,dot_x_end,dot_y_end):# cette pocédure trace la courbe selon un repère donné c'est la procédure la plus élémentale
		sign_y=1
		#purple_circle(dot_origin)
		#print
		if is_close(dot_origin[0],dot_y_end[0], 0.00001):
		
			#print "vertical side detected"
			
			if dot_y_end[1]<=dot_origin[1]:# this if is maybe useless
				#print "side toward up"
				#print "introducing sign - in trace"
				sign_y=-1*sign_y
			else:
				#print "side toward down"
				#print "introducing sign - in trace"
				sign_y=1*sign_y
		
		if is_close(dot_origin[1],dot_y_end[1], 0.00001):
			#print
			#print "horizontal side detected"
			
			if dot_y_end[1]<=dot_origin[1]:# this if is maybe useless
				#print "side toward right"
				#print "introducing sign - in trace"
				sign_y=-1*sign_y
			else:
				#print "side toward left"
				#print "introducing sign - in trace"
				sign_y=1*sign_y
	
		#print "sign",sign_y
		
		y_max_array=max(array)
		trace=[]
		full_x_length=dot_signed_distance_one_axis(dot_origin, dot_x_end,"x")
		#print "full_x_length",full_x_length
		if full_x_length<0:
			sign_x=-1
		else:
			sign_x=1
		full_y_length=dot_signed_distance_one_axis(dot_origin, dot_y_end,"y")
		#print "full_y_length",full_y_length
	
		inc=full_x_length/float(steps)
		#print "inc",inc
		inc_array=x_max_array/float(steps)
		#purple_circle(dot_origin)
		for i in range(steps+1):
		
			#print
			#print "i",i
			#place coords on axis
			dot_a=translate_absolute(dot_origin,dot_x_end, i*inc)#a is the dot on axe x
			
			
				
			
			distance_x=inc_array*i
			
			y=retreive_coordinates(array,distance_x)
			y_ratio=(y*(full_y_length/float(y_max_array)))*sign_y
			
			#mener une parallèle
			
			
			
			#red_big_circle(dot_origin)
			
			
			
			#a_coef_axis_y=line_equation(dot_x_end,dot_y_end)[1]
			#b1=dot_y_end[1]-a_coef_axis_y*dot_y_end[0]
			#dot_y_end_parallel=(dot_origin[0]+0.2,a_coef_axis_y*(dot_origin[0]+0.2)+b1)
			dot_b=translate_absolute(dot_origin,dot_y_end, y_ratio)
			
			
			
			
		#	red_circle(dot_b)
			
			
			
			if not is_close(dot_origin[1],dot_y_end[1], 0.00001) and not is_close(dot_origin[0],dot_y_end[0],0.00001):
				a_coef_axis_y=line_equation(dot_origin,dot_y_end)[1]
				b2=dot_a[1]-a_coef_axis_y*dot_a[0]
				#print a_coef_axis_y
				dot_a_prime=(dot_a[0]+0.2,a_coef_axis_y*(dot_a[0]+0.2)+b2)
				#green_circle(dot_a_prime)
			else:
				if is_close(dot_origin[1],dot_y_end[1], 0.00001):
					if dot_origin[0]<dot_x_end[0]:
						dot_a_prime=(dot_a[0]-0.2,dot_a[1])
					else:
						dot_a_prime=(dot_a[0]+0.2,dot_a[1])
					
				if is_close(dot_origin[0],dot_y_end[0],0.00001):
					if dot_origin[1]<dot_y_end[1]:
						dot_a_prime=(dot_a[0],dot_a[1]-0.2)
					else:
						dot_a_prime=(dot_a[0],dot_a[1]+0.2)
			#print "y",y
			dot_found=translate_absolute(dot_a,dot_a_prime, y_ratio)
			#red_circle(dot_found)
			#dot_found=translate_absolute(dot_a_prime,dot_a, y_ratio)
		
				
			trace.append(dot_found)
		return trace# outputs one list of dots
		
		
		
		
	if steps%2!=0: steps=steps-1
	exponential_segment=(middle(dot0,dot1),middle(middle(dot0,dot1),dot1))
	line_a=line_equation(dot0,dot1)[1]
	perp=perp_dot2(dot0,dot1, exponential_segment[0])#extracting line equation of perpendicular crossing the center of the segment
	distance=math.fabs(dot_signed_distance(dot0,dot1))*depth
	#print perp
	
	
	# ici on trouve les points du repère de destination de la courbe
	
	
	#radial y extremity
	if dot0[0]>dot1[0]:
		distance_radial=radius+distance
	else:
		distance_radial=radius+distance

	#print "distance radial",distance_radial
	#print "defining radial extremity"
	
	if dot1[0]<origin[0]:
		#print "inverting radius distance"
		distance_dot1=-distance_radial
	else:
		distance_dot1=distance_radial
	if dot0[0]<origin[0]:
		#print "inverting radius distance"
		distance_dot0=-distance_radial
	else:
		distance_dot0=distance_radial
		
	
	sign=direction
	y_end_1_radial=translate_absolute(origin,dot1,distance_dot1*sign)
	y_end_2_radial=translate_absolute(origin,dot0,distance_dot0*sign)
	
	
	#print "distance made (should be equal to distance demanded",dot_distance(origin,y_end_1_radial)
	#purple_circle(y_end_1_radial)
	#purple_circle(y_end_2_radial)
	#blue_big_circle(dot1)
	#blue_big_circle(dot0)
	#green_big_circle(origin)
	#print "radial extremity defined"
	
	
	
	if perp[2]:
		#print "perp horizontal"
		if dot1[0]<dot0[0]:
			sign=-1
		else:
			sign=1
		
		sign=sign*direction
		dot_high=(exponential_segment[0][0],exponential_segment[0][1]+distance*sign)
		#dot_high=translate_absolute(exponential_segment[0], perp_other_dot, distance)
		#purple_circle(dot_high)
	
		y_end1=(dot1[0],dot1[1]+distance*sign)
		y_end2=(dot0[0],dot0[1]+distance*sign)
		
	if perp[3]:
		#print "perp vertical "
		if dot1[1]<dot0[1]:
			sign=-1
		else:
			sign=1
			
			
		sign=sign*direction
		dot_high=(exponential_segment[0][0]+distance*sign,exponential_segment[0][1])
		#dot_high=translate_absolute(exponential_segment[0], perp_other_dot, distance)
		#red_circle(dot_high)
		#purple_circle(dot_high)
	
		y_end1=(dot1[0]+distance*sign,dot1[1])
		y_end2=(dot0[0]+distance*sign,dot0[1])
		#blue_big_circle(dot0)
		
	if not perp[2] and not perp[3]:#lorsque la ligne est en biais
		
		if exponential_segment[0][0]<origin[0]:
			sign=-1
		else:
			sign=1
		
		
		sign=sign*direction
		
		x=exponential_segment[0][0]+0.1#arbitrary new x for prolonging perpendicular
		perp_other_dot=(x,perp[0]*x+perp[1])
		#purple_circle(perp_other_dot)
		dot_high=translate_absolute(exponential_segment[0], perp_other_dot, distance*sign)
		#purple_circle(dot_high)
		# arrived here we have three points (one is along the bisectrice of the segment)
		b2=dot1[1]-perp[0]*dot1[0]
		b3=dot0[1]-perp[0]*dot0[0]
		y_end_step1=(dot1[0]+0.2,(perp[0]*(dot1[0]+0.2))+b2)
		y_end_step2=(dot0[0]+0.2,(perp[0]*(dot0[0]+0.2))+b3)
	
		y_end1=translate_absolute(dot1,y_end_step1,distance*sign)
		y_end2=translate_absolute(dot0,y_end_step2,distance*sign)
	
#	print curve
#	print dictionary_method[curve]
	if dictionary_method[curve]=="parallels":
		#print "parallels method chosen"
		y_end_choose=(y_end1,y_end2)
	if dictionary_method[curve]=="radial":
		y_end_choose=(y_end_1_radial,y_end_2_radial)
		#print "radial method chosen"
	#red_big_circle(y_end1)
	#red_big_circle(y_end2)
		
	#perp2=perp_dot2(dot0,dot1, dot1)
	#if type(perp2)==int:
	#	print "vertical"
	#x=dot1[0]+0.1
	#perp_other_dot2=(x,perp2[0]*x+perp2[1])
	#y_end=translate_absolute(dot1[0], perp_other_dot2, distance)
	return (trace_curve(curve_array, steps, len(curve_array)-1,dot_high,middle(dot0,dot1),y_end_choose[0]),
			trace_curve(curve_array, steps, len(curve_array)-1,dot_high,middle(dot0,dot1),y_end_choose[1]),
	(dot_high, middle(dot0,dot1),dot1))

	
def trace_accolade(curve, dot0,dot1,depth, steps,direction,color):#direction is minus one or one, depth is a distance
	#the origin of the graph is taken as a landmark for vector direction
	
	dot_list=[]
	global origin
	y_end_choose=((0,0),(0,0))
	# radial method has been temporarily abandonned
	
	dictionary_method={"petal":"parallels",
				"exponential":"parallels",
				"reverse exponential":"parallels"}
	if curve=="petal":
		curve_array=petal
	if curve=="exponential":
		curve_array=exponential
	if curve=="reverse exponential":
		curve_array=reverse_exponential
		
	def trace_curve(array,steps,x_max_array,dot_origin,dot_x_end,dot_y_end):# cette pocédure trace la courbe selon un repère donné c'est la procédure la plus élémentale
		sign_y=1
		#purple_circle(dot_origin)
		print
		if is_close(dot_origin[0],dot_y_end[0], 0.00001):
		
		#	print "vertical side detected"
			
			if dot_y_end[1]<=dot_origin[1]:# this if is maybe useless
				#print "side toward up"
				#print "introducing sign - in trace"
				sign_y=-1*sign_y
			else:
				#print "side toward down"
				#print "introducing sign - in trace"
				sign_y=1*sign_y
		
		if is_close(dot_origin[1],dot_y_end[1], 0.00001):
		#	print
		#	print "horizontal side detected"
			
			if dot_y_end[1]<=dot_origin[1]:# this if is maybe useless
				#print "side toward right"
				#print "introducing sign - in trace"
				sign_y=-1*sign_y
			else:
				#print "side toward left"
				#print "introducing sign - in trace"
				sign_y=1*sign_y
	
		#print "sign",sign_y
		
		y_max_array=max(array)
		trace=[]
		full_x_length=dot_signed_distance_one_axis(dot_origin, dot_x_end,"x")
		#print "full_x_length",full_x_length
		if full_x_length<0:
			sign_x=-1
		else:
			sign_x=1
		full_y_length=dot_signed_distance_one_axis(dot_origin, dot_y_end,"y")
		#print "full_y_length",full_y_length
	
		inc=full_x_length/float(steps)
		#print "inc",inc
		inc_array=x_max_array/float(steps)
		#purple_circle(dot_origin)
		for i in range(steps+1):
			dot_a=translate_absolute(dot_origin,dot_x_end, i*inc)#a is the dot on axe x
			
			
				
			
			distance_x=inc_array*i
			
			y=retreive_coordinates(array,distance_x)
			y_ratio=(y*(full_y_length/float(y_max_array)))*sign_y
			

			dot_b=translate_absolute(dot_origin,dot_y_end, y_ratio)
			if not is_close(dot_origin[1],dot_y_end[1], 0.00001) and not is_close(dot_origin[0],dot_y_end[0],0.00001):
					a_coef_axis_y=line_equation(dot_origin,dot_y_end)[1]
					b2=dot_a[1]-a_coef_axis_y*dot_a[0]
					#print a_coef_axis_y
					dot_a_prime=(dot_a[0]+0.2,a_coef_axis_y*(dot_a[0]+0.2)+b2)
					#green_circle(dot_a_prime)
			else:
				if is_close(dot_origin[1],dot_y_end[1], 0.00001):
					if dot_origin[0]<dot_x_end[0]:
						dot_a_prime=(dot_a[0]-0.2,dot_a[1])
					else:
						dot_a_prime=(dot_a[0]+0.2,dot_a[1])
					
				if is_close(dot_origin[0],dot_y_end[0],0.00001):
					if dot_origin[1]<dot_y_end[1]:
						dot_a_prime=(dot_a[0],dot_a[1]-0.2)
					else:
						dot_a_prime=(dot_a[0],dot_a[1]+0.2)
			#print "y",y
			dot_found=translate_absolute(dot_a,dot_a_prime, y_ratio)
			#print "adding dot---------"
			trace.append(dot_found)
		return trace# outputs one list of dots
		
		
		
		
	if steps%2!=0: steps=steps-1
	exponential_segment=(middle(dot0,dot1),middle(middle(dot0,dot1),dot1))
	line_a=line_equation(dot0,dot1)[1]
	perp=perp_dot2(dot0,dot1, exponential_segment[0])#extracting line equation of perpendicular crossing the center of the segment
	distance=math.fabs(dot_signed_distance(dot0,dot1))*depth
	#print perp
	#print "ACCOLDAE"
	
	# ici on trouve les points du repère de destination de la courbe
	
	
	#radial y extremity
	
	if perp[2]:
		#print "perp horizontal"
		if dot1[0]<dot0[0]:
			sign=-1
		else:
			sign=1
		
		sign=sign*direction
		dot_high=(exponential_segment[0][0],exponential_segment[0][1]+distance*sign)
		#purple_circle(dot_high)
	
		y_end1=(dot1[0],dot1[1]+distance*sign)
		y_end2=(dot0[0],dot0[1]+distance*sign)
		
	if perp[3]:
		#print "perp vertical "
		if dot1[1]<dot0[1]:
			sign=-1
		else:
			sign=1
			
			
		sign=sign*direction
		dot_high=(exponential_segment[0][0]+distance*sign,exponential_segment[0][1])
		y_end1=(dot1[0]+distance*sign,dot1[1])
		y_end2=(dot0[0]+distance*sign,dot0[1])
	
	if not perp[2] and not perp[3]:#lorsque la ligne est en biais
		
		if exponential_segment[0][0]<origin[0]:
			sign=-1
		else:
			sign=1
		
		
		sign=sign*direction
		
		x=exponential_segment[0][0]+0.1#arbitrary new x for prolonging perpendicular
		perp_other_dot=(x,perp[0]*x+perp[1])
		dot_high=translate_absolute(exponential_segment[0], perp_other_dot, distance*sign)
		b2=dot1[1]-perp[0]*dot1[0]
		b3=dot0[1]-perp[0]*dot0[0]
		y_end_step1=(dot1[0]+0.2,(perp[0]*(dot1[0]+0.2))+b2)
		y_end_step2=(dot0[0]+0.2,(perp[0]*(dot0[0]+0.2))+b3)
	
		y_end1=translate_absolute(dot1,y_end_step1,distance*sign)
		y_end2=translate_absolute(dot0,y_end_step2,distance*sign)
	
	#print curve
	#print dictionary_method[curve]
	if dictionary_method[curve]=="parallels":
		#print "parallels method chosen"
		y_end_choose=(y_end1,y_end2)
	if dictionary_method[curve]=="radial":
		y_end_choose=(y_end_1_radial,y_end_2_radial)
		#print "radial method chosen"
	
	Dots_1=trace_curve(curve_array, steps, len(curve_array)-1,dot_high,middle(dot0,dot1),y_end_choose[0])
	Dots_2=trace_curve(curve_array, steps, len(curve_array)-1,dot_high,middle(dot0,dot1),y_end_choose[1])
	
	
	#print "dots accolade",Dots_1
	for dot_index in range(len(Dots_1)-1):
		#print"tracing-----------------------------"
		lin(ctx, Dots_1[dot_index], Dots_1[dot_index+1], line_thickness, color)
	dot_index=0
	for dot_index in range(len(Dots_2)-1):
		lin(ctx, Dots_2[dot_index], Dots_2[dot_index+1], line_thickness, color)

	
def circle_cascade_2(iterations, center, start_circle_radius, end_circle_radius, curve, space_prop, radius_start_segment, n_sides, rotate,color):
	if start_circle_radius<end_circle_radius: 
		direction=1
	else:
		direction =-1
				
	print"direction",direction
	iterations=iterations+1
	radius_start=radius_start_segment+start_circle_radius/2.0
	#current_radius=radius_start
	prop_curve_list=[]
	cumulated_radii=0
	max_curve=max(curve)
	curve_length=len(curve)
	espacement=curve_length/float(iterations)
	for i in range(iterations):#applying curve to case
		prop_curve_list.append(retreive_coordinates(curve,espacement*i)/float(max_curve))
	prec_radius=0
	print prop_curve_list
	radius_list=[]
	if direction==-1:
		temp=start_circle_radius
		start_circle_radius=end_circle_radius
		end_circle_radius=temp
	for i in range(iterations-1):#creating radius list
		
		current_s_radius=start_circle_radius*prop_curve_list[i+1]+end_circle_radius*prop_curve_list[iterations-1-i]#increasing from center
		#"increasing from center"
		radius_list.append(current_s_radius)

	print "len prop list", len(prop_curve_list)
	
	if direction==-1:
		temp=start_circle_radius
		start_circle_radius=end_circle_radius
		end_circle_radius=temp
	current_radius=radius_list[0]
	print "radii list",radius_list
	for i in range(iterations-1):
		print "i",i
		prec_radius=current_radius
		if direction == -1:
			print "decreasing from center"
			current_radius=radius_list[i]#decreasing from center
		else:
			print "increasing from center"
			current_radius=radius_list[iterations-2-i]#increasing from center
		print "current radius",current_radius
		cumulated_radii=cumulated_radii+(current_radius)+prec_radius
		if direction==1:
			space=space_prop*(start_circle_radius*2.0)#space is thought in relation to diameter
		else:
			space=space_prop*(end_circle_radius*2.0)#space is thought in relation to diameter
		print "radius_start",radius_start
		print "space",space
		print "cumulated_radii",cumulated_radii
		main_radius=radius_start+cumulated_radii+space*i
		cycle_cycle(ctx,center,main_radius, current_radius, n_sides, 64,1, rotate,0, color)
	
def trace_accolade(curve, dot0,dot1,depth, steps,direction,color):#direction is minus one or one, depth is a distance
	#the origin of the graph is taken as a landmark for vector direction
	
	dot_list=[]
	global origin
	y_end_choose=((0,0),(0,0))
	# radial method has been temporarily abandonned
	
	dictionary_method={"petal":"parallels",
				"exponential":"parallels",
				"reverse exponential":"parallels"}
	if curve=="petal":
		curve_array=petal
	if curve=="exponential":
		curve_array=exponential
	if curve=="reverse exponential":
		curve_array=reverse_exponential
		
	def trace_curve(array,steps,x_max_array,dot_origin,dot_x_end,dot_y_end):# cette pocédure trace la courbe selon un repère donné c'est la procédure la plus élémentale
		sign_y=1
		#purple_circle(dot_origin)
		#print
		if is_close(dot_origin[0],dot_y_end[0], 0.00001):
		
			#print "vertical side detected"
			
			if dot_y_end[1]<=dot_origin[1]:# this if is maybe useless
				#print "side toward up"
				#print "introducing sign - in trace"
				sign_y=-1*sign_y
			else:
			#	print "side toward down"
			#	print "introducing sign - in trace"
				sign_y=1*sign_y
		
		if is_close(dot_origin[1],dot_y_end[1], 0.00001):
			#print
			#print "horizontal side detected"
			
			if dot_y_end[1]<=dot_origin[1]:# this if is maybe useless
				#print "side toward right"
			#	print "introducing sign - in trace"
				sign_y=-1*sign_y
			else:
				#print "side toward left"
				#print "introducing sign - in trace"
				sign_y=1*sign_y
	
		#print "sign",sign_y
		
		y_max_array=max(array)
		trace=[]
		full_x_length=dot_signed_distance_one_axis(dot_origin, dot_x_end,"x")
		#print "full_x_length",full_x_length
		if full_x_length<0:
			sign_x=-1
		else:
			sign_x=1
		full_y_length=dot_signed_distance_one_axis(dot_origin, dot_y_end,"y")
		#print "full_y_length",full_y_length
	
		inc=full_x_length/float(steps)
		#print "inc",inc
		inc_array=x_max_array/float(steps)
		#purple_circle(dot_origin)
		for i in range(steps+1):
			dot_a=translate_absolute(dot_origin,dot_x_end, i*inc)#a is the dot on axe x
			
			
				
			
			distance_x=inc_array*i
			
			y=retreive_coordinates(array,distance_x)
			y_ratio=(y*(full_y_length/float(y_max_array)))*sign_y
			

			dot_b=translate_absolute(dot_origin,dot_y_end, y_ratio)
			if not is_close(dot_origin[1],dot_y_end[1], 0.00001) and not is_close(dot_origin[0],dot_y_end[0],0.00001):
					a_coef_axis_y=line_equation(dot_origin,dot_y_end)[1]
					b2=dot_a[1]-a_coef_axis_y*dot_a[0]
					#print a_coef_axis_y
					dot_a_prime=(dot_a[0]+0.2,a_coef_axis_y*(dot_a[0]+0.2)+b2)
					#green_circle(dot_a_prime)
			else:
				if is_close(dot_origin[1],dot_y_end[1], 0.00001):
					if dot_origin[0]<dot_x_end[0]:
						dot_a_prime=(dot_a[0]-0.2,dot_a[1])
					else:
						dot_a_prime=(dot_a[0]+0.2,dot_a[1])
					
				if is_close(dot_origin[0],dot_y_end[0],0.00001):
					if dot_origin[1]<dot_y_end[1]:
						dot_a_prime=(dot_a[0],dot_a[1]-0.2)
					else:
						dot_a_prime=(dot_a[0],dot_a[1]+0.2)
			#print "y",y
			dot_found=translate_absolute(dot_a,dot_a_prime, y_ratio)
			#print "adding dot---------"
			trace.append(dot_found)
		return trace# outputs one list of dots
		
		
		
		
	if steps%2!=0: steps=steps-1
	exponential_segment=(middle(dot0,dot1),middle(middle(dot0,dot1),dot1))
	line_a=line_equation(dot0,dot1)[1]
	perp=perp_dot2(dot0,dot1, exponential_segment[0])#extracting line equation of perpendicular crossing the center of the segment
	distance=math.fabs(dot_signed_distance(dot0,dot1))*depth
	#print perp
#	print "ACCOLDE"
	
	# ici on trouve les points du repère de destination de la courbe
	
	
	#radial y extremity
	
	if perp[2]:
		#print "perp horizontal"
		if dot1[0]<dot0[0]:
			sign=-1
		else:
			sign=1
		
		sign=sign*direction
		dot_high=(exponential_segment[0][0],exponential_segment[0][1]+distance*sign)
		purple_circle(dot_high)
	
		y_end1=(dot1[0],dot1[1]+distance*sign)
		y_end2=(dot0[0],dot0[1]+distance*sign)
		
	if perp[3]:
		#print "perp vertical "
		if dot1[1]<dot0[1]:
			sign=-1
		else:
			sign=1
			
			
		sign=sign*direction
		dot_high=(exponential_segment[0][0]+distance*sign,exponential_segment[0][1])
		y_end1=(dot1[0]+distance*sign,dot1[1])
		y_end2=(dot0[0]+distance*sign,dot0[1])
	
	if not perp[2] and not perp[3]:#lorsque la ligne est en biais
		
		if exponential_segment[0][0]<origin[0]:
			sign=-1
		else:
			sign=1
		
		
		sign=sign*direction
		
		x=exponential_segment[0][0]+0.1#arbitrary new x for prolonging perpendicular
		perp_other_dot=(x,perp[0]*x+perp[1])
		dot_high=translate_absolute(exponential_segment[0], perp_other_dot, distance*sign)
		b2=dot1[1]-perp[0]*dot1[0]
		b3=dot0[1]-perp[0]*dot0[0]
		y_end_step1=(dot1[0]+0.2,(perp[0]*(dot1[0]+0.2))+b2)
		y_end_step2=(dot0[0]+0.2,(perp[0]*(dot0[0]+0.2))+b3)
	
		y_end1=translate_absolute(dot1,y_end_step1,distance*sign)
		y_end2=translate_absolute(dot0,y_end_step2,distance*sign)
	
	#print curve
	#print dictionary_method[curve]
	if dictionary_method[curve]=="parallels":
		#print "parallels method chosen"
		y_end_choose=(y_end1,y_end2)
	if dictionary_method[curve]=="radial":
		y_end_choose=(y_end_1_radial,y_end_2_radial)
		#print "radial method chosen"
	
	Dots_1=trace_curve(curve_array, steps, len(curve_array)-1,dot_high,middle(dot0,dot1),y_end_choose[0])
	Dots_2=trace_curve(curve_array, steps, len(curve_array)-1,dot_high,middle(dot0,dot1),y_end_choose[1])
	
	
	#print "dots accolade",Dots_1
	for dot_index in range(len(Dots_1)-1):
		#print"tracing-----------------------------"
		lin(ctx, Dots_1[dot_index], Dots_1[dot_index+1], line_thickness, color)
	dot_index=0
	for dot_index in range(len(Dots_2)-1):
		lin(ctx, Dots_2[dot_index], Dots_2[dot_index+1], line_thickness, color)

	

def exponential_concentric(ctx, center,start_radius,fraction_of_radius, segment_de_reference, base_radius_outline, iterations, jump,rotate,color):
# start_radius is the outside radius
#segment of reference is the width of the whole rythm it's expressed in a fraction of the start radius
#base_radius_outline is an increment of start_radius expressed in the same fraction (can be practical to change the radius with a convenient unit (typically this value may be 0)
#iterations is the number of lines
#jump is what dots are joined by lines
#rotate is a global rotate to align peaks

	span=len(exponential)/float(iterations+2)
	y_array=[]
	
	
	
	for i in range(iterations):
		offset=(i+1)*span
		if offset>len(exponential)-1:
			offset=len(exponential)-1
		y_array.append(retreive_coordinates(exponential,offset))
	def h(y,maxi,apply_to):
		h=y*apply_to/maxi
		return h
	
	hauteur=[]
	for i in range(iterations):
		hauteur.append(segment_de_reference-h(y_array[i],max(exponential),segment_de_reference))

	cycle(ctx, center, 144, start_radius+base_radius_outline*start_radius*fraction_of_radius,jump,rotate,color)#base_cycle
	
	
	for i in range(iterations):
		cycle(ctx, center, 144, start_radius+(base_radius_outline+hauteur[i])*start_radius*fraction_of_radius,jump,rotate,color)

	

def re_order(figure):
	line_list=[]
	def line_not_in(line_list, tested):
		for line in line_list:
			if (line[0]==tested[0] and line[1]==tested[1]) or (line[1]==tested[0] and line[0]==tested[1]):
				return False
				#print "rejecting double line"
		return True
	def find_door(figure):
		door_list=[]
		for i in range(len(figure)):#cycling thru lines
			for j in range(len(figure[i])):#cycling thru dots
				found_bound=False
				for k in range(len(figure)):#cycling thru lines for each dot
					if is_close_dot(figure[i][j],figure[k][0],0.0001) or is_close_dot(figure[i][j],figure[k][1],0.0001):
						found_bound=True
				if not found_bound:
					door_list.append(figure[i][j])
		return door_list
		
	def find_dot(figure,dot, exclude):
		for i in range (len(figure)):
			for j in range (len(figure[i])):
				#print "exclude",exclude
				if is_close_dot(figure[i][j], dot, 0.0001) and (i,j)!=exclude:
					return (i,j)
		return ()
		
	def find_direction(dot_coords_in_line):
		result=1
		if dot_coords_in_line[1]==1:
			result=0
		else:
			result=1
		return result
		
	def invert (line):
		return (line[1],line[0])
		
	door_list=find_door(figure)
	#print "door_list",door_list
	if len(door_list)>0:
		door=door_list[0]# unbound parts of the figure are becoming bound (else, it is needed to make a figure list)
	else:
		#print "figure is closed"
		door=(0,0)
	#next step is to be able to change direction (each next dot must associate a direction)
	#print"following normal direction"
	line_list.append(figure[door[0]])# adding first line
#	print" adding first line",door[0]
	start=True
	find_bound=[]
	direction=find_direction(door)
	next_point=figure[door[0]][direction]
	count=0#panic count
	while (len(find_bound)>0 or start==True) and (next_point!=figure[door[0]][1] or start==True):#either we reach an end or the first point
		#print
		count+=1
		if start:
			next_point=figure[door[0]][direction]
			direction=find_direction(door)
			last_find_bound=door
		else:
			pass
		#	print"getting serious"
		#print "last_find_bound before",last_find_bound
		#print "excluding",(last_find_bound[0],direction)
		#print "searching for", next_point
		find_bound=find_dot(figure,next_point,(last_find_bound[0],direction))#outputs coordinates
		direction=find_direction(find_bound)
			#print "new find bound",find_bound
		#print "direction",direction
		# when a new line is found which is in the other direction there is a stall
		#print "last_find_bound after",last_find_bound
		
		start=False
		next_point=figure[find_bound[0]][direction]
		direction=find_direction(find_bound)
		#print "next point",next_point
		last_find_bound=find_bound
		line_to_append=figure[find_bound[0]]
		if len(find_bound)>0 and line_not_in(line_list,line_to_append):
			if direction==0:
				line_to_append=invert(line_to_append)
			line_list.append(line_to_append)# adding new line
			#print " adding new line",find_bound[0]
		else:
			pass
			#print "found no match"
	return line_list
				
	
	

#-------------------------------------------------------------main mandala application ---------------------------------------------------------------------------------
	
WIDTH=2500
HEIGHT=2500




# defining curve tables----------------------

exponential_2=[]

for i in range (64):
	exponential_2.append(i/10.0*i/50.0)
	
exponential=[]

for i in range (64):
	exponential.append(i/10.0*i/5.0)
	
expmax=max(exponential)
reverse_exponential=[]
for i in range (len(exponential)):
	reverse_exponential.append((expmax-exponential[len(exponential)-1-i]))

#print reverse_exponential

petal=[]
for i in range(len(exponential)+len(reverse_exponential)):
	if i<=(len(exponential)-1):
		petal.append(exponential[i])
	else:
		petal.append(expmax+reverse_exponential[i-(len(exponential))])
		
#print petal

origin=(0.5,0.5)
line_thickness=0.009/5.0#most thick line

Harmonizer_radius=1/6.0
	
#preparing canvas ------------------------------------------------------
Path_perso="/home/raphael/"
name="seven"
filename=u""+Path_perso+name+u"_Harmonizer.svg"
filename_small=u""+Path_perso+name+u"_Harmonizer_small.png"

surface = cairo.SVGSurface(filename, WIDTH, HEIGHT)
ctx = cairo.Context (surface)
ctx.scale (WIDTH, HEIGHT) # Normalizing the canvas

# creating white background
ctx.rectangle (0, 0, 1, 1) # Rectangle(x0, y0, x1, y1)
color="697189"
color="ffffff"
ctx.set_source_rgb (hexcolor(color)[0], hexcolor(color)[1], hexcolor(color)[2]) # Solid color
#ctx.set_source_rgb (1, 1, 1) # Solid color
ctx.fill () #action to fill previous rectangle with previously selected color	




#----------------------------------- here starts the drawing
#drawing to work on outlines

origin=(0.5,0.5)
line_thickness=0.009/7.0#most thick line
Harmonizer_radius=1/21.0
space=0.009/14.0

n_sides=8

#circle_cascade_2(21,origin,space*2,space*10,exponential,0.5,Harmonizer_radius,n_sides,0,"000000")
#show=ram_losanges(origin,Harmonizer_radius,n_sides,
	#											0.3,#width
	#											0.5,#depth
	#											0.5,0)

#for i in range(len(show)):
#	show2=outline_ram(show[i],space*50,False)
#	trace_figure(show[i],"000000",True)
#	trace_figure(show2,"000000",True)

stroke_color4="000000"
stroke_color3="000000"
stroke_color2="000000"
stroke_color="000000"
#radii
#print "space before",space-
#test outline
n_sides=7

show=ram_cycle_accolade(origin,Harmonizer_radius*5,"petal",n_sides,5,1,0.618,1,180+90)
#print "original figure",show
trace_figure(show,stroke_color3,False)
show3=re_order(show)#outputs opened shape
#print "re ordered",show3
for i in range(10):
	show4=outline_ram(show3,(i)*(space*6),False)#outputs an opened figure event if input is closed
	trace_figure(show4,stroke_color3,True)
	pass
for i in range(10):
	#show4=outline_ram(show3,-(i+1)*(space*6),False)#outputs an opened figure event if input is closed
	#trace_figure(show4,stroke_color3,True)
	pass

#show=ram_cycle_accolade(origin,Harmonizer_radius*3,"petal",n_sides,21,1,0.618,1,180+90)
#print "original figure",show
#show3=re_order(show)#outputs opened shape
#print "re ordered",show3
for i in range(7):
	pass
	#show4=outline_ram(show3,(i)*(space*6),False)#outputs an opened figure event if input is closed
	#trace_figure(show4,stroke_color3,True)
for i in range(7):
	pass
	#show4=outline_ram(show3,-(i+1)*(space*6),False)#outputs an opened figure event if input is closed
	#trace_figure(show4,stroke_color3,True)

#show=ram_cycle_accolade(origin,Harmonizer_radius,"exponential",n_sides,8,1,0.618,1,180+90)
#print "original figure",show
#show3=re_order(show)#outputs opened shape
#print "re ordered",show3
for i in range(14):
	pass
	#show4=outline_ram(show3,(i)*(space*6),False)#outputs an opened figure event if input is closed
	#trace_figure(show4,stroke_color3,True)
for i in range(0):
	pass
	#show4=outline_ram(show3,-(i+1)*(space*6),False)#outputs an opened figure event if input is closed
	#trace_figure(show4,stroke_color3,True)
#print "space after",space
#trace_figure(show[0],stroke_color,False)
#trace_figure(show2,stroke_color2,True)
#trace_figure(show3,stroke_color4,False)

#trace_figure(show4,stroke_color3,True)

#for i in range (72):
#	show2=outline_ram(show,-space*(i+1)*7,False)
#	trace_figure(show2,stroke_color,True)

#print "figure_list:",show
#print "figure_one:",show[0]
#print "coutour_output",show_2
# Now ends the drawing, saving file
#cycle_dotted(ctx,origin, 5, Harmonizer_radius/1.618, line_thickness*3, 1/12.0 ,2,0)


#-------------------here ends the drawing

#Path_perso="/home/raphael/"
#name="seven"
#filename=u""+Path_perso+name+u"_Harmonizer.png"
#filename_small=u""+Path_perso+name+u"_Harmonizer_small.png"

#filename_small = unidecode.unidecode(filename_small)	
#filename = unidecode.unidecode(filename)	
#print "filename type:",type(filename)
print "saving ",filename

ctx.save () # Output to svg
#im=Image.open(filename)
#im=im.resize((600, 600), Image.ANTIALIAS)
#im.save(filename_small)

