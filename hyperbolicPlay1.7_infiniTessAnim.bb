Graphics 300,300

R = 250

Type point
	Field id
	Field x#,y#
	Field u#,v#
	Field tu#, tv#
	Field r,g,b
	Field links.point[100]
	Field numlinks
	Field mouseControlled
	Field show
	Field toDelete
End Type

Type interPoint
	Field x#, y#
	Field p1x#, p1y#
	Field p2x#, p2y#
End Type

Type camera
	Field id
	Field x#, y#
	Field orient#
	Field snap.point
	Field yparity
	
	;virtual camera that ignores boundaries
	Field vx#, vy#
	Field vorient#
End Type

Type boundary
	Field x1#,y1#
	Field x2#,y2#
End Type

n = 3
k = 7

Global numFrames = 0
Global frameLimit = 38

init(n,k)
draw(R)

SetBuffer(BackBuffer())
While Not KeyHit(1)

	If getInput(R)
	
		Cls
	
		numFrames = numFrames + 1
		draw(R)
	
		Flip
	
	EndIf

Wend
End

Function init(n,k)

;	n = polygons around a point
	theta# = 360./n
	
;	k = sides per polygon
	phi# = 360./k
	
	numer# = (1+Cos(phi))*Cos(theta/2)
	denom# = Sin(phi)*Sin(theta/2)

	dis# = Sqr( (numer/denom)^2 - 1 )
	rad# = dis*(Sqr(1+dis^2)*(1-Cos(phi))*Cos(theta/2) + Sin(phi)*Sin(theta/2))

	distanceLimit = 50^2
	limit = 1000
	tolerance# = 0.1
	
	Local queue.point[10000]
	Local ao#[10000] ;angle offset
	center.point = Null
	p2.point = Null
	np.point = Null
	
	q = 2
	queue[1] = createPoint(0,0)
	queue[2] = createPoint(0,-rad)
	queue[2]\numLinks = 1
	queue[2]\links[1] = queue[1]
	ao[1] = -90
	ao[2] = 90
	
	cam.camera = New camera
	cam\id = 1
	cam\snap = First point
	cam\x = 0.01
	cam\y = 0.01
	
	i = 0
	
	While i < q
		i = i + 1
		
		center = queue[i]
		angoff# = ao[i]
		
		For j = 1 To n-1 ;skips origin point
			Local nXY#[2]
			ang# = (j*theta + angoff) Mod 360
			x2# = rad*Cos(ang)
			y2# = rad*Sin(ang)
			
			translate(x2,y2, -center\x,-center\y, nXY)
			
			coincide = 0
			For k = 1 To q
				If KeyHit(1) End
				
				p2 = queue[k]
				If Abs(p2\x - nXY[1]) < tolerance And Abs(p2\y - nXY[2]) < tolerance
					center\numLinks = center\numLinks + 1
					center\links[ center\numLinks ] = p2
					coincide = 1
					Exit
				EndIf
			Next
			
			If coincide = 0 And q < limit And (nXY[1]^2+nXY[2]^2) < distanceLimit
				np = createPoint(nXY[1],nXY[2])
				np\numLinks = 1
				np\links[1] = center
				
				Local offXY#[2]
				translate(center\x,center\y, np\x,np\y, offXY)
				
				q = q+1
				queue[q] = np
				ao[q] = ATan2(offXY[2],offXY[1])
			EndIf
			
			If KeyDown(28) Return
			
		Next
	Wend
	
	
	Local bounds.point[100]
	off# = 0.5*(1 - n Mod 2)
	
	;create boundary lines
	For i = 0 To n-1
		b.boundary = New boundary
		b\x1 = dis*Cos(90+(i+off)*theta)
		b\y1 = dis*Sin(90+(i+off)*theta)
		b\x2 = dis*Cos(90+(i+1+off)*theta)
		b\y2 = dis*Sin(90+(i+1+off)*theta)
	Next
	
End Function

Function createPoint.point(x#,y#, r=255,g=255,b=255, show=0, toDelete=0)

	p.point = New point
	p\x = x
	p\y = y
	
	p\u = transform(p\x,p\y, 1)
	p\v = transform(p\x,p\y, 2)
	
	lastPoint.point = Last point
	p\id = lastPoint\id + 1
	
	p\r = r
	p\g = g
	p\b = b
	
	p\show = show
	p\toDelete = toDelete
	
	Return p

End Function

Function getInput(R)

	active = 0

	mx = MouseX()
	my = MouseY()
	mkey = 0
		
	cam.camera = First camera
	
	tx# = 0
	ty# = 0
	scrollSpeed# = 0.10
	k# = cosh(scrollSpeed)
	
	mult = 1 - 2*cam\yparity
	
	If KeyDown(32) Or KeyDown(205)
		tx = -scrollSpeed
	ElseIf KeyDown(30) Or KeyDown(203)
		tx =  scrollSpeed
	EndIf
	If KeyDown(17) Or KeyDown(200)
		ty = -scrollSpeed * mult
	ElseIf KeyDown(31) Or KeyDown(208)
		ty =  scrollSpeed * mult
	EndIf
	
	If KeyDown(16)
		cam\orient = (cam\orient + 5*mult + 360) Mod 360
		cam\vorient = (cam\vorient + 5) Mod 360
		active = 1
	ElseIf KeyDown(18)
		cam\orient = (cam\orient - 5*mult + 360) Mod 360
		cam\vorient = (cam\vorient - 5 + 360) Mod 360
		active = 1
	EndIf
	
	If numFrames < frameLimit
		ty = -scrollSpeed * mult
	EndIf
	
	If tx <> 0 Or ty <> 0

		d# = Sqr(tx*tx + ty*ty)
		ang# = ATan2(ty,tx)
		
		cam\vx = offsetStep(cam\vx,cam\vy, d,ang+cam\vorient, 1)
		cam\vy = offsetStep(cam\vx,cam\vy, d,ang+cam\vorient, 2)
		
		nx# = offsetStep(cam\x,cam\y, d,ang+1*cam\orient, 1)
		ny# = offsetStep(cam\x,cam\y, d,ang+1*cam\orient, 2)
		
		Local out#[2]
		
		done = 0
		crossed.boundary = Null
		pick.boundary = Null
		
		While Not done
			done = 1
			max# = 10000
			pick = Null
			
			For b.boundary = Each boundary
				If b <> crossed
					intersection(b\x1,b\y1, b\x2,b\y2, cam\x,cam\y, nx,ny, out)
					
					If out[0] = 1 ;there was an intersection
						done = 0
						d# = hyperD(out[1],out[2], cam\x,cam\y)

						If pick = Null
							max = d
							pick = b
						ElseIf d < max
							max = d
							pick = b
						EndIf
					EndIf
				EndIf
			Next
			
			If Not done
				intersection(pick\x1,pick\y1, pick\x2,pick\y2, cam\x,cam\y, nx,ny, out)

				cam\x = out[1]
				cam\y = out[2]
				
				flipOverLine(nx,ny, out[1],out[2], pick\x2,pick\y2, out)
				
				nx = out[1]
				ny = out[2]
				
				ang3# = hyperAng(cam\x,cam\y, pick\x1,pick\y1)
				
				cam\orient = (2*ang3 - cam\orient + 3600) Mod 360
				cam\yparity = Not cam\yparity
				
				crossed = pick
			EndIf

		Wend
		
		If iters > 0
			FlushKeys
		EndIf
		
		cam\x = nx
		cam\y = ny
		
		active = 1
	EndIf
	
	Return active
	
End Function

Function draw(R)

	gw = GraphicsWidth()
	gh = GraphicsHeight()	

	Color 255,255,255
	circ(gw/2,gh/2, R, 0)
	
	cam.camera = First camera

	For p.point = Each point
		If KeyHit(57) Stop
		
		Local tXY#[2]
		translate(p\x,p\y, cam\x,cam\y, tXY)
		
		mult = 1 - 2*cam\yparity
		
		tu# = -transform(tXY[1],tXY[2], 1)
		tv# = mult * transform(tXY[1],tXY[2], 2)
		
		p\tu = tu*Cos(cam\orient*mult) - tv*Sin(cam\orient*mult)
		p\tv = tu*Sin(cam\orient*mult) + tv*Cos(cam\orient*mult)
		
		If p\show
			Color p\r,p\g,p\b
			circ(R*p\tu+gw/2,R*p\tv+gh/2, 3, 1)
		EndIf
		If p\toDelete
			Delete p
		EndIf
	Next
	
	Local vals#[3]
	
	For p1.point = Each point
		For idx = 1 To p1\numlinks
			p2.point = p1\links[idx]
			
			u1# = p1\tu
			v1# = p1\tv
			u2# = p2\tu
			v2# = p2\tv
			
			denom# = u1*v2-u2*v1
			
			If KeyHit(57) And Abs(p1\x) < 0.1 And Abs(p1\y) = 0.1
				Stop
			EndIf
			
			If Abs(denom) < 0.001
				Color 255,255,255
				p1tx = R * p1\tu + gw/2
				p1ty = R * p1\tv + gh/2
				p2tx = R * p2\tu + gw/2
				p2ty = R * p2\tv + gh/2
				
				Line p1tx,p1ty, p2tx,p2ty
			Else
		
				f# = (u1*u1+v1*v1+1)/2.
				g# = (u2*u2+v2*v2+1)/2.
			
				If Abs(v2) > 0.0001
					h# = (v2*f-v1*g)/denom
					k# = (g-u2*h)/v2
				ElseIf Abs(v1) > 0.0001
					h# = (v2*f-v1*g)/denom
					k# = (f-u1*h)/v1
				Else
					Stop
				EndIf
				
				rad# = Sqr(h*h+k*k-1)
				
				angMin# = ATan2(v1-k,u1-h)
				angMax# = ATan2(v2-k,u2-h)
				angDiff# = (angMax-angMin)
				
				If angDiff * 2 <> angDiff
					
					While angDiff < -180 Or angDiff > 180
						If angDiff < -180
							angDiff = angDiff + 360
						Else
							angDiff = angDiff - 360
						EndIf
					Wend
					
					numSteps = 20
					angStep# = angDiff/numSteps
					
					Color 255,255,255
					For n = 0 To numSteps-1
						ang# = angMin + n*angStep
						
						x1 = R*(h+rad*Cos(ang)) + GraphicsWidth()/2
						y1 = R*(k+rad*Sin(ang)) + GraphicsHeight()/2
						x2 = R*(h+rad*Cos(ang+angStep)) + GraphicsWidth()/2
						y2 = R*(k+rad*Sin(ang+angStep)) + GraphicsHeight()/2
						
						Line x1,y1, x2,y2
						
					Next
				
				EndIf
				
			EndIf
			;end of line drawing
			
		Next

	Next

	Color 255,255,0
	Line gw/2-2,gh/2, gw/2+2,gh/2
	Line gw/2,gh/2-2, gw/2,gh/2+2
	
	SaveBuffer(BackBuffer(),"Hyperbolic infinite tessellation frames/hyperbolic_infiniteTessellation_frame"+numFrames+".png")

End Function

Function offsetStep#(px#,py#, d#,theta#, which)
	k# = (Exp(d)+1/Exp(d))/2
	
	If px = 0 And py = 0
		r# = Sqr(k*k-1)
		If which = 1
			Return r*Cos(theta)
		Else
			Return r*Sin(theta)
		EndIf
	EndIf
	
	a# = 1 + px^2*Sin(theta)^2 + py^2*Cos(theta)^2 - 2*px*py*Cos(theta)*Sin(theta)
	b# = 2*(1-k)*(px*Cos(theta) + py*Sin(theta))
	c# = (1 - k*k + 2*(1-k)*(px*px + py*py) )
		
	det# = b*b-4*a*c
	rpos# = (-b + Sqr(det))/(2*a)
	rneg# = (-b - Sqr(det))/(2*a)
	
	If rpos < 0
		temp# = rpos
		rpos = rneg
		rneg = temp
	EndIf
	
	If which = 1
		Return px + rpos*Cos(theta)
	ElseIf which = 2
		Return py + rpos*Sin(theta)
	EndIf

End Function

Function translate#(px#,py#, tx#,ty#, out#[2]) ;tx,ty to origin

	pt# = Sqr(1.+px*px+py*py)
	tt# = Sqr(1.+tx*tx+ty*ty)
	
	If tx = 0 And ty = 0
		out[1] = px
		out[2] = py
		Return
	EndIf
	
	x1# = tt
	x2# = tx
	x3# = ty
	
	a# =  x1
	d# = -x2
	g# = -x3

	b# = -x2
	c# = -x3

	f# = (x1*x2*x3 - (x2*x3))/(x2*x2+x3*x3)
	h# = f
	
	e# = (x1*x2*x2 + x3*x3)/(x2*x2+x3*x3)
	i# = (x1*x3*x3 + x2*x2)/(x2*x2+x3*x3)

	;matrix is done!
	
	nt# = a*pt + b*px + c*py
	nx# = d*pt + e*px + f*py
	ny# = g*pt + h*px + i*py
	
	out[0] = nt
	out[1] = nx
	out[2] = ny
	
	Return

End Function

Function translatePoint(p.point, tx#,ty#)
	Local out#[2]
	translate(p\x,p\y, tx,ty, out)
	p\x = out[1]
	p\y = out[2]
End Function

Function transform#(x#,y#, which)

	t# = Sqr(1+x*x+y*y)

	If which = 1
		Return x/(1.+t)
	ElseIf which = 2
		Return y/(1.+t)
	EndIf

End Function

Function invTransform#(u#,v#, which)

	denom# = 1.-u*u-v*v

	If which = 1
		Return 2.*u/denom
	ElseIf which = 2
		Return 2.*v/denom
	EndIf

End Function

Function polygon(centerX#,centerY#, dis#, sides, angoffset#)

	angstep# = 360./sides
	rad# = cosh(dis)
	
	prevPoint.point = Last point
	idStart = prevPoint\id
	prevPoint = Null
	
	begPoint.point = Null

	For i = 0 To sides-1
	
		p.point = New point

		tempx# = rad*Cos(i*angstep + angoffset)
		tempy# = rad*Sin(i*angstep + angoffset)
		
		Local newXY#[2]
		translate(tempx,tempy, -centerX,-centerY, newXY)
		
		p\x = newXY[1]
		p\y = newXY[2]
		
		p\r = 255
		p\g = 255
		p\b = 255
		
		p\id = idStart + 1
		idStart = idStart + 1
		
		If prevPoint <> Null
			prevPoint\links[1] = p
			prevPoint\numLinks = 1
		EndIf
		
		prevPoint = p
	
	Next
	
	prevPoint\numLinks = 1
	For j = 1 To sides-1
		p = Before p
	Next
	prevPoint\links[1] = p
	
	p\r = 255
	p\g = 0
	p\b = 0

End Function

Function polygon_oriented(centerX#,centerY#, orientX#,orientY#, rad#, sides)

	angstep# = 360./sides
	
	prevPoint.point = Last point
	idStart = prevPoint\id
	prevPoint = Null
	
	begPoint.point = Null
	
	Local orientT_XY#[2]
	translate(orientX,orientY, -centerX,-centerY, orientT_XY)
	
	angoffset# = (ATan2(orientT_XY[2], orientT_XY[1]) + 180*((sides+0) Mod 2)) Mod 360

	For i = 0 To sides-1
	
		p.point = New point
		tempx# = rad*Cos(i*angstep + angoffset)
		tempy# = rad*Sin(i*angstep + angoffset)
		
		Local newXY#[2]
		translate(tempx,tempy, -centerX,-centerY, newXY)
		
		p\x = newXY[1]
		p\y = newXY[2]
		
		p\r = 255
		p\g = 255
		p\b = 255
		
		p\id = idStart + 1
		idStart = idStart + 1
		
		If prevPoint <> Null
			prevPoint\links[1] = p
			prevPoint\numLinks = 1
		EndIf
		
		prevPoint = p
	
	Next
	
	prevPoint\numLinks = 1
	For j = 1 To sides-1
		p = Before p
	Next
	prevPoint\links[1] = p
	
	p\r = 255
	p\g = 0
	p\b = 0

End Function

Function intersection(x1#,y1#, x2#,y2#,  u1#,v1#, u2#,v2#, XY#[2], strict=1) ;strict checks the line >segments<

	;XY[0] will be used to store failure (0) or success (1)
	XY[0] = 0

	If (x1 = u1 And y1 = v1 And x2 = u2 And y2 = v2) Or (x1 = u2 And y1 = v2 And x2 = u1 And y2 = v1)
		Return
	EndIf

	Local tXY_1#[2]
	Local tXY_2#[2]
	Local tXY_3#[2]
	
	translate(x2,y2, x1,y1, tXY_1) ;translate all points so that (x1,y1) = (0,0)
	translate(u1,v1, x1,y1, tXY_2)
	translate(u2,v2, x1,y1, tXY_3)
	
	ang# = ATan2(tXY_1[2],tXY_1[1]) ;calculate angle to (x2,y2)
	
	tU1# = tXY_2[1]*Cos(-ang) - tXY_2[2]*Sin(-ang) ;rotate (u1,v1) by -ang
	tV1# = tXY_2[1]*Sin(-ang) + tXY_2[2]*Cos(-ang)
	
	tU2# = tXY_3[1]*Cos(-ang) - tXY_3[2]*Sin(-ang) ;rotate (u2,v2) by -ang
	tV2# = tXY_3[1]*Sin(-ang) + tXY_3[2]*Cos(-ang)
	
	s# = tU1
	t# = tV1
	u# = tU2
	v# = tV2
	
	tst# = Sqr(1 + s*s + t*t)
	tuv# = Sqr(1 + u*u + v*v)
	
	p# = s*v - u*t
	q# = t*tuv - v*tst
	
	If q*q - p*p <= 0 Return
	
	k# = -p*Sgn(q)/Sqr(q*q-p*p)
	
	;intersection point is at (k,0)
	
	kX# = k*Cos(ang) ;rotate solution by ang
	kY# = k*Sin(ang)
	
	translate(kX,kY, -x1,-y1, XY) ;translate (0,0) to (x1,y1)
	
	;line segment test if desired
	If strict = 1
		If k >= 0 And k <= Sqr(tXY_1[1]^2 + tXY_1[2]^2) And tV1*tV2 <= 0
			XY[0] = 1
		EndIf
	Else
		XY[0] = 1
	EndIf
	
	Return

End Function

Function flipOverLine(px#,py#, x1#,y1#, x2#,y2#, out#[2])

	Local tXY_1#[2]
	Local tXY_2#[2]
	
	translate(x2,y2, x1,y1, tXY_1)
	translate(px,py, x1,y1, tXY_2)
	
	angx# = ATan2(tXY_1[2], tXY_1[1])
	angp# = ATan2(tXY_2[2], tXY_2[1])
	angdiff# = 2*(angx - angp)
	
	nx# = tXY_2[1]*Cos(angdiff) - tXY_2[2]*Sin(angdiff)
	ny# = tXY_2[1]*Sin(angdiff) + tXY_2[2]*Cos(angdiff)
	
	translate(nx,ny, -x1,-y1, out)

End Function

Function D#(x1#,y1#, x2#,y2#)
	Return Sqr((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
End Function

Function hyperD#(x1#,y1#, x2#,y2#)
	Return acosh( Sqr( (1+x1*x1+y1*y1)*(1+x2*x2+y2*y2) ) - x1*x2 - y1*y2 )
End Function

Function hyperAng(x1#,y1#, x2#,y2#)
	Local out#[2]
	translate(x2,y2, x1,y1, out)
	Return ATan2(out[2],out[1])
End Function

Function circ(x,y,r, fill=0)
	Oval x-r,y-r,2*r+1,2*r+1, fill
End Function

Function cosh#(x#)
	Return (Exp(x)+1./Exp(x))/2
End Function
Function sinh#(x#)
	Return (Exp(x)-1./Exp(x))/2
End Function

Function acosh#(z#)
	Return Log(z + Sqr(z*z-1))
End Function