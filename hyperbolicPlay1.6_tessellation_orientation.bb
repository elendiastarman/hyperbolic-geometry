Graphics 700,700

R = 300

Type point
	Field id
	Field x#,y#
	Field u#,v#
	Field tu#, tv#
	Field r,g,b
	Field links.point[100]
	Field numlinks
	Field fresh
End Type

Type camera
	Field id
	Field x#, y#
	Field u#, v#
	Field du#, dv#
	Field orientation#
	Field snap.point
End Type

init(3,7)
;init(7,3, 1)
init(7,3,1)
draw(R)

SetBuffer(BackBuffer())
While Not KeyHit(1)

	If getInput(R)
	
		Cls
	
		draw(R)
	
		Flip
	
	EndIf

Wend
End

Function init(n,k, offset#=0)

;	n = 3 ;polygons around a point
	theta# = 360./n
	
;	k = 7 ;sides per polygon
	phi# = 360./k
	
	numer# = (1+Cos(phi))*Cos(theta/2)
	denom# = Sin(phi)*Sin(theta/2)

	dis# = Sqr( (numer/denom)^2 - 1 )
;	rad# = cosh(dis)
;	rad2# = acosh(dis)
;	DebugLog "dis = "+dis+", rad = "+rad+", rad2 = "+rad2+", cosh(rad2) = "+cosh(rad2)
	
;	DebugLog invTransform(0.95,0,1)+","+invTransform(0.95,0,2)

;	rad# = (1+(dis^2)*(1-Cos(phi)))
;	rad# = 2*dis^2*(1-Cos(phi)) + dis^4*(1-Cos(phi))^2
;;	rad# = dis*(Sqr(1+dis^2)*(1-Cos(phi))*Cos(theta/2) + Sin(phi)*Sin(theta/2))
	rad# = Sqr( (1+(dis^2)*(1-Cos(phi)))^2 - 1 )

	distanceLimit = 8^2
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
	
	queue[1]\fresh = 1
	queue[2]\fresh = 1
	
	cam.camera = New camera
	cam\id = 1
	cam\snap = First point
	
	;draw(R)
	;WaitKey
	
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
			
			translateXY(x2,y2, -center\x,-center\y, nXY)
			
			;Cls
			;Color 255,0,0
			;circ(-transform(nXY[1],nXY[2],1)*R+GraphicsWidth()/2,transform(nXY[1],nXY[2],2)*R+GraphicsWidth()/2, 5,0)
			
			;x3# = rad*Cos(angoff)
			;y3# = rad*Sin(angoff)
			;Local nXY2#[2]
			;translateXY(x3,y3, -center\x,-center\y, nXY2)
			;Color 0,0,255
			;circ(-transform(nXY2[1],nXY2[2],1)*R+GraphicsWidth()/2,transform(nXY2[1],nXY2[2],2)*R+GraphicsWidth()/2, 7,0)
			
			coincide = 0
			For k2 = 1 To q
				If KeyHit(1) End
				
				p2 = queue[k2]
				If Abs(p2\x - nXY[1]) < tolerance And Abs(p2\y - nXY[2]) < tolerance And p2\fresh = 1
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
				np\fresh = 1
				
				Local offXY#[2]
				translateXY(center\x,center\y, np\x,np\y, offXY)
				
				q = q+1
				queue[q] = np
				ao[q] = ATan2(offXY[2],offXY[1])
				
			;	DebugLog "Point created! "+q
			EndIf
			
;			If q > limit
;				Return
;			EndIf
			
			;draw(R)
			;WaitKey
			
			If KeyDown(28) Return
			
		Next
	Wend
	
	DebugLog offset
	;Stop
	
	offdis# = -offset*Sqr( ( ((1+Cos(theta))*Cos(phi/2))/(Sin(theta)*Sin(phi/2)) )^2 - 1 )
	DebugLog offdis

	offang# = offset * theta/2
	DebugLog offang

	For p.point = Each point
		If p\fresh = 1
			rx# = p\x*Cos(offang) - p\y*Sin(offang)
			ry# = p\x*Sin(offang) + p\y*Cos(offang)
			
			;DebugLog (p\x*p\x+p\y*p\y) - (rx*rx+ry*ry)
			
			;Line p\x+GraphicsWidth()/2,p\y+GraphicsHeight()/2, rx+GraphicsWidth()/2,ry+GraphicsHeight()/2
			;Delay 10
			
			Local tXY#[2]
			translateXY(rx,ry, 0,offdis, tXY)
			
			p\x = tXY[1]
			p\y = tXY[2]
			
			p\fresh = 0
		EndIf
	Next
	
End Function

Function createPoint.point(x#,y#, r=255,g=255,b=255)

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
	
	Return p

End Function

Function getInput(R)

	active = 0

	mx = MouseX()
	my = MouseY()
	mkey = 0
	
	tx# = 0
	ty# = 0
	scrollSpeed# = 0.10
	k# = (Exp(scrollSpeed)+1/Exp(scrollSpeed))/2 ;cosh(scrollSpeed)
	
	If KeyDown(32) Or KeyDown(205)
		tx = -scrollSpeed
	ElseIf KeyDown(30) Or KeyDown(203)
		tx =  scrollSpeed
	EndIf
	If KeyDown(17) Or KeyDown(200)
		ty = -scrollSpeed
	ElseIf KeyDown(31) Or KeyDown(208)
		ty =  scrollSpeed
	EndIf
		
	cam.camera = First camera
	
	If KeyHit(49)
		;Stop
		If cam\snap = First point
			cam\snap = Last point
		Else
			cam\snap = Before cam\snap
		EndIf
		cam\x = cam\snap\x
		cam\y = cam\snap\y
		active = 1
	ElseIf KeyHit(50)
		If cam\snap = Last point
			cam\snap = First point
		Else
			cam\snap = After cam\snap
		EndIf
		cam\x = cam\snap\x
		cam\y = cam\snap\y
		active = 1
	EndIf
	
	If KeyDown(16)
		cam\orientation = cam\orientation + 10
		active = 1
	ElseIf KeyDown(18)
		cam\orientation = cam\orientation - 10
		active = 1
	EndIf
	
	If tx <> 0 Or ty <> 0
		
;		offX# = offsetStep(cam\x,cam\y, tx,90*Sgn(tx)-90, 1)
;		offY# = offsetStep(cam\x,cam\y, tx,90*Sgn(tx)-90, 2)
;		
;		offX2# = offsetStep(offX,offY, ty,90*Sgn(ty), 1)
;		offY2# = offsetStep(offX,offY, ty,90*Sgn(ty), 2)
;		
;		cam\x = offX2
;		cam\y = offY2

		d# = Sqr(tx*tx + ty*ty)
		ang# = ATan2(ty,tx) - cam\orientation
		cam\x = offsetStep(cam\x,cam\y, d,ang, 1)
		cam\y = offsetStep(cam\x,cam\y, d,ang, 2)
		
;		Local newXY#[2]
;		translateXY(cam\x,cam\y, -offX,-offY, newXY)
		
;		cam\x = newXY[1]
;		cam\y = newXY[2]
		
		active = 1
	EndIf
	
	;If MouseXSpeed() Or MouseYSpeed() Return 1
	
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
		translateXY(p\x,p\y, cam\x,cam\y, tXY)
		
		rx# = tXY[1]*Cos(cam\orientation) - tXY[2]*Sin(cam\orientation)
		ry# = tXY[1]*Sin(cam\orientation) + tXY[2]*Cos(cam\orientation)
		
		p\tu# = -transform(rx,ry, 1)
		p\tv# = transform(rx,ry, 2)
		
;		Color p\r,p\g,p\b
;		circ(R*p\tu+gw/2,R*p\tv+gh/2, 3, 1)
	Next
	
	Local vals#[3]
	
;	Return
	
	For p1.point = Each point
		For idx = 1 To p1\numlinks
			p2.point = p1\links[idx]
			
			u1# = p1\tu
			v1# = p1\tv
			u2# = p2\tu
			v2# = p2\tv
			
			denom# = u1*v2-u2*v1
			
			If KeyHit(57) And Abs(p1\x) < 0.1 And Abs(p1\y) = 0.1; And idx = 4
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
			
		Next
		
	Next
	
	Color 255,255,0
	Line gw/2-2,gh/2, gw/2+2,gh/2
	Line gw/2,gh/2-2, gw/2,gh/2+2
	
	Color 128,128,0
	cx = R*cam\u+gw/2
	cy = R*cam\v+gh/2
	Line cx-2,cy, cx+2,cy
	Line cx,cy-2, cx,cy+2

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

Function translate#(pu#,pv#, tu#,tv#, which) ;(pu,pv) = point to translate, (tu,tv) = point to move to origin

		vdotx# = pu*tu + pv*tv ;tu*p\u + tv*p\v
		mag2v# = tu*tu + tv*tv ;tu*tu + tv*tv
		mag2x# = pu*pu + pv*pv ;p\u*p\u + p\v*p\v
		
		If which = 1
			Return ( (1.+2*vdotx+mag2x)*tu + (1.-mag2v)*pu ) / (1. + 2*vdotx + mag2v*mag2x)
		ElseIf which = 2
			Return ( (1.+2*vdotx+mag2x)*tv + (1.-mag2v)*pv ) / (1. + 2*vdotx + mag2v*mag2x)
		EndIf

End Function

Function translateXY#(px#,py#, tx#,ty#, out#[2]) ;tx,ty to origin

	pt# = Sqr(1+px*px+py*py)
	tt# = Sqr(1+tx*tx+ty*ty)
	
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

Function transform#(x#,y#, which)

	;If x*x+y*y > 1
	;	Return -1
	;EndIf

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
	
;	centerU# = transform(centerX,centerY, 1)
;	centerV# = transform(centerX,centerY, 2)

	For i = 0 To sides-1
	
		p.point = New point
		;tempx# = offsetStep(0,0, rad, i*angstep + angoffset, 1)
		;tempy# = offsetStep(0,0, rad, i*angstep + angoffset, 2)
		tempx# = rad*Cos(i*angstep + angoffset)
		tempy# = rad*Sin(i*angstep + angoffset)
		
		Local newXY#[2]
		translateXY(tempx,tempy, -centerX,-centerY, newXY)
		
		p\x = newXY[1]
		p\y = newXY[2]
		
		p\r = 255
		p\g = 255
		p\b = 255
		
		p\id = idStart + 1
		idStart = idStart + 1
		
		;If begPoint <> Null
		;	begPoint = p
		;EndIf
		
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

Function polygon_oriented(centerX#,centerY#, orientX#,orientY#, dis#, sides)

	angstep# = 360./sides
	rad# = dis;cosh(dis)
	
	prevPoint.point = Last point
	idStart = prevPoint\id
	prevPoint = Null
	
	begPoint.point = Null
	
;	centerU# = transform(centerX,centerY, 1)
;	centerV# = transform(centerX,centerY, 2)

	Local orientT_XY#[2]
	translateXY(orientX,orientY, -centerX,-centerY, orientT_XY)
	
	angoffset# = (ATan2(orientT_XY[2], orientT_XY[1]) + 180*((sides+0) Mod 2)) Mod 360

	For i = 0 To sides-1
	
		p.point = New point
		;tempx# = offsetStep(0,0, rad, i*angstep + angoffset, 1)
		;tempy# = offsetStep(0,0, rad, i*angstep + angoffset, 2)
		tempx# = rad*Cos(i*angstep + angoffset)
		tempy# = rad*Sin(i*angstep + angoffset)
		
		Local newXY#[2]
		translateXY(tempx,tempy, -centerX,-centerY, newXY)
		
		p\x = newXY[1]
		p\y = newXY[2]
		
		p\r = 255
		p\g = 255
		p\b = 255
		
		p\id = idStart + 1
		idStart = idStart + 1
		
		;If begPoint <> Null
		;	begPoint = p
		;EndIf
		
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

Function polygon_old(centerX#,centerY#, rad#, sides, angoffset#)

	angstep# = 360./sides
	
	prevPoint.point = Last point
	idStart = prevPoint\id
	prevPoint = Null
	
	begPoint.point = Null
	
	centerU# = transform(centerX,centerY, 1)
	centerV# = transform(centerX,centerY, 2)

	For i = 0 To sides-1
	
		p.point = New point
		tempx# = offsetStep(0,0, rad, i*angstep + angoffset, 1)
		tempy# = offsetStep(0,0, rad, i*angstep + angoffset, 2)
		
;		p\x = tempx ;translate(tempx,tempy, -centerX,-centerY, 1)
;		p\y = tempy ;translate(tempx,tempy, -centerX,-centerY, 2)
		
		tempu# = transform(tempx,tempy, 1)
		tempv# = transform(tempx,tempy, 2)
		
		p\u = translate(tempu,tempv, centerU,centerV, 1)
		p\v = translate(tempu,tempv, centerU,centerV, 2)
		
		p\x = invTransform(p\u,p\v, 1)
		p\y = invTransform(p\u,p\v, 2)
		
		p\r = 255
		p\g = 255
		p\b = 255
		
		p\id = idStart + 1
		idStart = idStart + 1
		
		;If begPoint <> Null
		;	begPoint = p
		;EndIf
		
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

End Function

Function D#(x1,y1, x2,y2)
	Return Sqr((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
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