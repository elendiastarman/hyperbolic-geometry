Graphics 700,700

R = 300

Type point
	Field id
	Field x#,y#
	Field u#,v#
	Field tu#, tv#
	Field r,g,b
	Field links.point[4]
	Field numlinks
End Type

Type camera
	Field id
	Field x#, y#
	Field u#, v#
	Field du#, dv#
	Field snap.point
End Type

init(R);, 7)
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

Function init(R)

	n = 0

;	num = 10
	dis# = 1.15
;	ang = 0
	
;	startX = -10
;	startY = 10
	
	p.point = New point
	p\x = startX
	p\y = startY
	p\id = 1
	
	prevPoint.point = p
	lastPoint.point = Null

	For k = 0 To 3
		newX# = offsetStep(0,0, dis,k*90, 1)
		newY# = offsetStep(0,0, dis,k*90, 2)
		
;		p.point = New point
;		p\x = newX
;		p\y = newY
		
;		p\u = transform(p\x,p\y, 1)
;		p\v = transform(p\x,p\y, 2)
		
;		lastPoint = Last point
;		p\id = lastPoint\id + 1
		
;		p\r = 255
;		p\g = 255
;		p\b = 255
		
;		p\numLinks = 1
;		p\links[1] = prevPoint
;		prevPoint = p
		
		polygon(newX,newY, dis, 6, k*30)
	Next
	
	cam.camera = New camera
	cam\id = 1
	cam\snap = First point
	
End Function

Function getInput(R)

	mx = MouseX()
	my = MouseY()
	mkey = 0
	
	tx# = 0
	ty# = 0
	scrollSpeed# = 0.10
	k# = (Exp(scrollSpeed)+1/Exp(scrollSpeed))/2 ;cosh(scrollSpeed)
	
	If KeyDown(32) Or KeyDown(205)
		tx =  scrollSpeed
	ElseIf KeyDown(30) Or KeyDown(203)
		tx = -scrollSpeed
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
		Return 1
	ElseIf KeyHit(50)
		If cam\snap = Last point
			cam\snap = First point
		Else
			cam\snap = After cam\snap
		EndIf
		cam\x = cam\snap\x
		cam\y = cam\snap\y
		Return 1
	EndIf
	
	If tx <> 0 Or ty <> 0
;		DebugLog "k = "+k
		
;		DebugLog "cam\x = "+cam\x+", cam\y = "+cam\y
		
		dx# = 0
		dy# = 0
		
;		If cam\y = 0
;			ax# = 1
;		Else
;			ax# = cam\y*cam\y
;		EndIf
		ax# = 1 + cam\y*cam\y
		bx# = 2*(1-k)*cam\x
		
;		If cam\x = 0
;			ay# = 1
;		Else
;			ay# = cam\x*cam\x
;		EndIf
		ay# = 1 + cam\x*cam\x
		by# = 2*(1-k)*cam\y
		
		c# = (1-k*k + 2*(1-k) * (cam\x*cam\x + cam\y*cam\y) )
		
;		DebugLog "a = "+a+", b = "+b+", c = "+c
		
;		dxpos# = 0
;		dxneg# = 0
;		dypos# = 0
;		dyneg# = 0
		
		If tx <> 0
			dxpos# = (-bx + Sqr(bx*bx-4*ax*c))/(2*ax)
			dxneg# = (-bx - Sqr(bx*bx-4*ax*c))/(2*ax)
;			DebugLog "dxpos = "+dxpos+", dxneg = "+dxneg
		Else
			dxpos = 0
			dxneg = 0
		EndIf
		
		If ty <> 0
			dypos# = (-by + Sqr(by*by-4*ay*c))/(2*ay)
			dyneg# = (-by - Sqr(by*by-4*ay*c))/(2*ay)
;			DebugLog "dypos = "+dypos+", dyneg = "+dyneg
		Else
			dypos = 0
			dyneg = 0
		EndIf
		
;		If a*dxpos*dxpos + b*dxpos + c = 0
;			DebugLog "dxpos passes the test!"
;		ElseIf a*dxneg*dxneg + b*dxneg + c = 0
;			DebugLog "dxneg passes the test!"
;		ElseIf a*dypos*dypos + b*dypos + c = 0
;			DebugLog "dypos passes the test!"
;		ElseIf a*dyneg*dyneg + b*dyneg + c = 0
;			DebugLog "dyneg passes the test!"
;		EndIf

;		DebugLog "dxpos check: "+(ax*dxpos*dxpos + bx*dxpos + c)
;		DebugLog "dxneg check: "+(ax*dxneg*dxneg + bx*dxneg + c)
;		DebugLog "dypos check: "+(ay*dypos*dypos + by*dypos + c)
;		DebugLog "dyneg check: "+(ay*dyneg*dyneg + by*dyneg + c)
		
		If dxpos < 0
			temp# = dxpos
			dxpos = dxneg
			dxneg = temp
		EndIf
		If dypos < 0
			temp# = dypos
			dypos = dyneg
			dyneg = temp
		EndIf
		
		If tx > 0
			cam\x = cam\x + dxpos
		Else
			cam\x = cam\x + dxneg
		EndIf
		If ty > 0
			cam\y = cam\y + dypos
		Else
			cam\y = cam\y + dyneg
		EndIf

;		DebugLog "cam\x = "+cam\x+", cam\y = "+cam\y
;		DebugLog ""
		
		Return 1
	EndIf
	
	If MouseXSpeed() Or MouseYSpeed() Return 1
	
	Return 0
	
End Function

Function draw(R)

	gw = GraphicsWidth()
	gh = GraphicsHeight()

	Color 255,255,255
	circ(gw/2,gh/2, R, 0)
	
	cam.camera = First camera

	For p.point = Each point
		If KeyHit(57) Stop
		
		Local testTXY#[2]
		translateXY(p\x,p\y, cam\x,cam\y, testTXY)
		
		p\tu# = transform(testTXY[1],testTXY[2], 1)
		p\tv# = transform(testTXY[1],testTXY[2], 2)
		
		Color p\r,p\g,p\b
		circ(R*p\tu+gw/2,R*p\tv+gh/2, 3, 1)
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

Function translateXY#(px#,py#, tx#,ty#, out#[2])

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

Function polygon(centerX#,centerY#, rad#, sides, angoffset#)

	angstep# = 360./sides
	
	prevPoint.point = Last point
	idStart = prevPoint\id
	prevPoint = Null
	
	begPoint.point = Null
	
;	centerU# = transform(centerX,centerY, 1)
;	centerV# = transform(centerX,centerY, 2)

	For i = 0 To sides-1
	
		p.point = New point
		tempx# = offsetStep(0,0, rad, i*angstep + angoffset, 1)
		tempy# = offsetStep(0,0, rad, i*angstep + angoffset, 2)
		
		Local newXY#[2]
		translateXY(tempx,tempy, centerX,centerY, newXY)
		
		p\x = newXY[1]
		p\y = newXY[2]
		
;		p\x = tempx ;translate(tempx,tempy, -centerX,-centerY, 1)
;		p\y = tempy ;translate(tempx,tempy, -centerX,-centerY, 2)
		
;		tempu# = transform(tempx,tempy, 1)
;		tempv# = transform(tempx,tempy, 2)
		
;		p\u = translate(tempu,tempv, centerU,centerV, 1)
;		p\v = translate(tempu,tempv, centerU,centerV, 2)
		
;		p\x = invTransform(p\u,p\v, 1)
;		p\y = invTransform(p\u,p\v, 2)
		
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

Function linearSystem3(M#[9], R#[3], out#[3])
	;solves [M|R] and stores results in out. I.E...
	;
	;__                     __    __              __
	;| M[1] M[2] M[3] | R[1] |    | 1 0 0 | out[1] |
	;| M[4] M[5] M[6] | R[2] | => | 0 1 0 | out[2] |
	;| M[7] M[8] M[9] | R[3] |    | 0 0 1 | out[3] |
	;__                     __    __              __
	;
	;;; WARNING: this does NOT check whether the system is solvable in the first place!

	a# = M[1]
	b# = M[2]
	c# = M[3]
	d# = M[4]
	e# = M[5]
	f# = M[6]
	g# = M[7]
	h# = M[8]
	i# = M[9]
	
	u# = R[1]
	v# = R[2]
	w# = R[3]

	;

	b2# = b/a
	c2# = c/a
	u2# = u/a
	
	e2# = e - d*b2
	f2# = f - d*c2
	v2# = v - d*u2
	
	h2# = h - g*b2
	i2# = i - g*c2
	w2# = w - g*u2
	
	;
	
	f3# = f2/e2
	v3# = v2/e2
	
	i3# = i2 - h2*f3
	w3# = w2 - h2*v3
	
	w4# = w3/i3
	
	z# = w4
	y# = v3 - z*f3
	x# = u2 - y*b2 - z*c2
	
	out[1] = x
	out[2] = y
	out[3] = z
	
	Return

End Function