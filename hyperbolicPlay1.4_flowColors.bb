Graphics 700,700

gw = GraphicsWidth()
gh = GraphicsHeight()

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
End Type

init(R);, 7)

draw(R)

Color 255,255,255
For y = 0 To gh
	For x = 0 To gw
		If (x-gw/2)^2 + (y-gh/2)^2 < R*R
			;DebugLog x+y
		
			px# = (1.*x-gw/2)/R
			py# = (1.*y-gh/2)/R
			
			itx# = invTransform(px,py, 1)
			ity# = invTransform(px,py, 2)
			
			nx# = offsetStep(itx,ity, 0.1,0, 1)
			ny# = offsetStep(itx,ity, 0.1,0, 2)
			
			drawx# = transform(nx,ny,1)*R + gw/2
			drawy# = transform(nx,ny,2)*R + gh/2
			
			ang# = ATan2(drawy-y, drawx-x)
			
			red = col(ang,0,0)
			green = col(ang,0,1)
			blue = col(ang,0,2)
			
			Color red,green,blue
			Plot x,y
			
;			Color 255,0,0
;			circ(x,y, 1,1)
;			
;			Color 255,255,255
;			Line x,y, drawx,drawy
		
		EndIf
	Next
Next

Color 255,255,255
circ(gw/2,gh/2, R,0)

WaitKey
End

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

	num = 0
	mult = 4
	
	makeLinks = 1

	;For u# = -1 To 1.01 Step 0.2
	;	For v# = -1 To 1.01 Step 0.2
	For y = -num To num
		For x = -num To num
			n = n + 1
			p.point = New point
			p\id = n
			
			p\x = x*mult
			p\y = y*mult

			p\u = transform(p\x,p\y, 1)
			p\v = transform(p\x,p\y, 2)
			
			;p\x = p\u*R + GraphicsWidth()/2
			;p\y = p\v*R + GraphicsHeight()/2
			p\r = 255
			p\g = 255
			p\b = 255
			
			If makeLinks
				For p2.point = Each point
					If (p2\id = p\id - 1 And x > -num) Or (p2\id = p\id - num*2-1 And y > -num)
						p\numlinks = p\numlinks + 1
						p\links[ p\numlinks ] = p2
					EndIf
				Next
			EndIf
		Next
	Next
	
	cam.camera = New camera
	cam\id = 1
	
End Function

Function getInput(R)

	mx = MouseX()
	my = MouseY()
	mkey = 0
	
	tx# = 0
	ty# = 0
	scrollSpeed# = 0.20
	k# = (Exp(scrollSpeed)+1/Exp(scrollSpeed))/2 ;cosh(scrollSpeed)
	
	If KeyDown(32) Or KeyDown(205)
		tx = -scrollSpeed
	ElseIf KeyDown(30) Or KeyDown(203)
		tx =  scrollSpeed
	EndIf
	If KeyDown(17) Or KeyDown(200)
		ty =  scrollSpeed
	ElseIf KeyDown(31) Or KeyDown(208)
		ty = -scrollSpeed
	EndIf
	
	If tx <> 0 Or ty <> 0
;		DebugLog "k = "+k
		
		cam.camera = First camera
		
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
	;circ(GraphicsWidth()/2,GraphicsHeight()/2, 2, 1)
	
	cam.camera = First camera
	
	cam\u = transform(cam\x,cam\y, 1); + cam\dv
	cam\v = transform(cam\x,cam\y, 2); + cam\du
	
	;tempu# = translate(cam\u,cam\v, cam\du,cam\dv, 1)
	;tempv# = translate(cam\u,cam\v, cam\du,cam\dv, 2)
	
	;cam\u = tempu
	;cam\v = tempv
	
	;cam\x = invTransform(cam\u,cam\v, 1)
	;cam\y = invTransform(cam\u,cam\v, 2)

	For p.point = Each point
		Color p\r,p\g,p\b
		
;		pt# = Sqr(1 + p\x*p\x + p\y*p\y)
;		camt# = Sqr(1 + cam\x*cam\x + cam\y*cam\y)
;		
;		B# = pt*camt - p\x*cam\x - p\y*cam\y
;
;		ang# = ATan2(cam\y-p\y, cam\x-p\x)
;		slope# = (1.*cam\y-p\y)/(cam\x-p\x)
;		D# = Sqr(1 + (cam\x-p\x)*(cam\x-p\x) + (cam\y-p\y)*(cam\y-p\y))
;		mult# = B/D
;
;		pxt# = (p\x-cam\x)*mult; * Sgn(cam\x-p\x)
;		pyt# = (p\y-cam\y)*mult;*slope ;* Sgn(cam\y-p\y)
;		
;		pxt# = cam\x-p\x
;		pyt# = cam\y-p\y
;
;		p\tu = transform(pxt,pyt, 1)
;		p\tv = transform(pxt,pyt, 2)
		
		p\tu = translate(p\u,p\v, cam\u,cam\v, 1)
		p\tv = translate(p\u,p\v, cam\u,cam\v, 2)

		drawx = R * p\tu + gw/2
		drawy = R * p\tv + gh/2
		
		circ(drawx,drawy, 3, 1)
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
	
	;circle around mouse
	omx = MouseX()-gw/2
	omy = MouseY()-gh/2
	Color 0,255,255
	Line MouseX()-2,MouseY(), MouseX()+2,MouseY()
	Line MouseX(),MouseY()-2, MouseX(),MouseY()+2
	;DebugLog "omx = "+omx+", omy = "+omy
	
	If omx*omx + omy*omy < R*R
		;DebugLog "Ding!"
		mcRad# = 1
		mx# = invTransform(1.*omx/R,1.*omy/R,1)
		my# = invTransform(1.*omx/R,1.*omy/R,2)
		;DebugLog "mx = "+mx+", my = "+my
		
		;k# = (Exp(mcRad)+1/Exp(mcRad))/2
		For theta = 0 To 360-5 Step 5
;			hx2# = mx + rpos*Cos(theta)
;			hy2# = my + rpos*Sin(theta)
;			
;			hx3# = mx + rneg*Cos(theta)
;			hy3# = my + rneg*Sin(theta)

			hx2# = offsetStep(mx,my, mcRad,theta, 1)
			hy2# = offsetStep(mx,my, mcRad,theta, 2)
			
			hx3# = offsetStep(mx,my, mcRad,theta+5, 1)
			hy3# = offsetStep(mx,my, mcRad,theta+5, 2)
			
			;Stop
			
			drawx1# = transform(hx2,hy2,1)*R + gw/2
			drawy1# = transform(hx2,hy2,2)*R + gh/2
			
			drawx2# = transform(hx3,hy3,1)*R + gw/2
			drawy2# = transform(hx3,hy3,2)*R + gh/2
			
;			Color 255,0,0
;			Plot drawx1*R + gw/2, drawy1*R + gh/2
;			
;			Color 0,0,255
;			Plot drawx2*R + gw/2, drawy2*R + gh/2

			Color 255,255,255
			circ(drawx1,drawy1, 1,1)
;			Plot drawx1,drawy1

			If theta < 180
				Color 0,0,255
			Else
				Color 255,0,0
			EndIf
			Line drawx1,drawy1, drawx2,drawy2
			
			If theta = 0 Or theta = 180
				Color 0,255,0
				Line omx+gw/2,omy+gh/2, drawx1,drawy1
			EndIf
			
;			Color 255,255,255
;			Plot drawx1,drawy1
		Next
	EndIf
	
	Color 255,255,0
	Line gw/2-2,gh/2, gw/2+2,gh/2
	Line gw/2,gh/2-2, gw/2,gh/2+2

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

Function translate#(pu#,pv#, tu#,tv#, which)

		vdotx# = pu*tu + pv*tv ;tu*p\u + tv*p\v
		mag2v# = tu*tu + tv*tv ;tu*tu + tv*tv
		mag2x# = pu*pu + pv*pv ;p\u*p\u + p\v*p\v
		
		If which = 1
			Return ( (1.+2*vdotx+mag2x)*tu + (1.-mag2v)*pu ) / (1. + 2*vdotx + mag2v*mag2x)
		ElseIf which = 2
			Return ( (1.+2*vdotx+mag2x)*tv + (1.-mag2v)*pv ) / (1. + 2*vdotx + mag2v*mag2x)
		EndIf

End Function

Function getYpoints(p1.point,p2.point, x, vals[3])

	u1 = p1\x - GraphicsWidth()/2
	v1 = p1\y - GraphicsHeight()/2
	u2 = p2\x - GraphicsWidth()/2
	v2 = p2\y - GraphicsHeight()/2

	A = 1
	B = (v1*(u1*u1+u2*u2)-u1*(v1*v1+v2*v2)+v1-u1)/(u1*v2-u2*v1)
	C = 1 + x*x + x * (u2*(v1*v1+v2*v2)-v2*(u1*u1+u2*u2)+u2-v2)/(u1*v2-u2*v1)
	
	det = B*B-4*A*C ;determinant
	If det < 0
		vals[0] = 0
		Return
	ElseIf det = 0
		vals[0] = 1
		vals[1] = -B/(2*A)
		Return
	ElseIf det > 0
		vals[0] = 2
		vals[1] = (-B - Sqr(det))/(2*A)
		vals[2] = (-B + Sqr(det))/(2*A)
		Return
	EndIf

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

Function D#(x1,y1, x2,y2)
	Return Sqr((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
End Function

Function circ(x,y,r, fill=0)
	Oval x-r,y-r,2*r+1,2*r+1, fill
End Function

Function col(ang1#, ang2#, kind)

	While ang1 > 360
		ang1 = ang1 - 360
	Wend
	While ang1 < 0 
		ang1 = ang1 + 360
	Wend

	While ang2 > 180
		ang2 = ang2 - 360
	Wend
	While ang2 < -180
		ang2 = ang2 + 360
	Wend
	
	a3# = ang2/180.
	
;	DebugLog kind
	
	If ang1 > 300
		diff# = (ang1-300)/60.
		r = 255
		g = 0
		b = 255*(1-diff)
	ElseIf ang1 > 240
		diff# = (ang1-240)/60.
		r = 255*diff
		g = 0
		b = 255
	ElseIf ang1 > 180
		diff# = (ang1-180)/60.
		r = 0
		g = 255*(1-diff)
		b = 255
	ElseIf ang1 > 120
		diff# = (ang1-120)/60.
		r = 0
		g = 255
		b = 255*diff
	ElseIf ang1 > 60
		diff# = (ang1-60)/60.
		r = 255*(1-diff)
		g = 255
		b = 0
	Else
		diff# = (ang1-00)/60.
		r = 255
		g = 255*diff
		b = 0
	EndIf
	
	If a3 > 0
		r2 = r + a3*(255-r)
		g2 = g + a3*(255-g)
		b2 = b + a3*(255-b)
	Else
		r2 = r + a3*r
		g2 = g + a3*g
		b2 = b + a3*b
	EndIf
	
;	DebugLog "First: " + r2 + ":" + g2 + ":" + b2
	
	If r2 > 255
		r2 = 255
	ElseIf r2 < 0
		r2 = 0
	EndIf

	If g2 > 255
		g2 = 255
	ElseIf g2 < 0
		g2 = 0
	EndIf

	If b2 > 255
		b2 = 255
	ElseIf b2 < 0
		b2 = 0
	EndIf
	
;	DebugLog "Second: " + r2 + ":" + g2 + ":" + b2
	
	If kind = 0
		Return r2
	ElseIf kind = 1
		Return g2
	ElseIf kind = 2
		Return b2
	Else
		Return 0
	EndIf
	
	If KeyHit(1) End

End Function