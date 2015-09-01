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

Type person
	Field id
	Field x#, y#
	Field u#, v#
End Type

Global frameNum = 0

init(R);, 7)
draw(R)
;SaveBuffer(FrontBuffer(), "hyperbolicMovementDemo_frame"+n+".png")

;SetBuffer(BackBuffer())
;While Not KeyHit(1)
;
;	If getInput(R)
;	
;		Cls
;	
;		draw(R)
;	
;		Flip
;	
;	EndIf
;
;Wend

;per.person = First person

move(40, 0, 0.05, R)
move(40, -0.05, 0, R)
move(80, 0, -0.05, R)
move(80, 0.05, 0, R)
move(80, 0, 0.05, R)
move(40, -0.05, 0, R)
move(40, 0, -0.05, R)

End

Function save()
	SaveBuffer(FrontBuffer(), "hyperbolicMovementDemo_frame"+frameNum+".png")
	frameNum = frameNum+1
End Function

Function move(num, dx#, dy#, R)
	per.person = First person
	For i = 1 To num
		If KeyHit(1) End
		per\x = per\x + dx
		per\y = per\y + dy
		Cls
		draw(R)
		save()
	Next
End Function

Function init(R)

	n = 0

	;For u# = -1 To 1.01 Step 0.2
	;	For v# = -1 To 1.01 Step 0.2
	For y = -10 To 10 Step 2
		For x = -10 To 10 Step 2
			n = n + 1
			p.point = New point
			p\id = n
			
			p\x = x
			p\y = y

			p\u = transform(x,y, 1)
			p\v = transform(x,y, 2)
			
			;p\x = p\u*R + GraphicsWidth()/2
			;p\y = p\v*R + GraphicsHeight()/2
			p\r = 255
			p\g = 255
			p\b = 255
			
			For p2.point = Each point
				If (p2\id = p\id - 1 And x > -10) Or (p2\id = p\id - 11 And y > -10)
					p\numlinks = p\numlinks + 1
					p\links[ p\numlinks ] = p2
				EndIf
			Next
		Next
	Next
	
	per.person = New person
	per\id = 1
	
End Function

Function getInput(R)

	mx = MouseX()
	my = MouseY()
	mkey = 0
	
	tx# = 0
	ty# = 0
	scrollSpeed# = 0.05
	
	;If per\state = "start"
	
	Return 0
	
End Function

Function draw(R)

	gw = GraphicsWidth()
	gh = GraphicsHeight()

	Color 255,255,255
	circ(gw/2,gh/2, R, 0)
	;circ(GraphicsWidth()/2,GraphicsHeight()/2, 2, 1)
	
	per.person = First person

	For p.point = Each point
		Color p\r,p\g,p\b
		
		;is NOT correct
;		drawx = R * transform(p\x-per\x,p\y-per\y, 1) + gw/2
;		drawy = R * transform(p\x-per\x,p\y-per\y, 2) + gh/2

		p\u = transform(p\x,p\y, 1)
		p\v = transform(p\x,p\y, 2)
		
		per\u = transform(per\x,per\y, 1)
		per\v = transform(per\x,per\y, 2)
		
		p\tu = translate(p\u,p\v, per\u,per\v, 1)
		p\tv = translate(p\u,p\v, per\u,per\v, 2)

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
			
				If Abs(v2) > 0.001
					h# = (v2*f-v1*g)/denom
					k# = (g-u2*h)/v2
				ElseIf Abs(v1) > 0.001
					h# = (v2*f-v1*g)/denom
					k# = (f-u1*h)/v1
				Else
					Stop
				EndIf
				
				rad# = Sqr(h*h+k*k-1)
				
				angMin# = ATan2(v1-k,u1-h)
				angMax# = ATan2(v2-k,u2-h)
				angDiff# = (angMax-angMin)
				
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
			
		Next
		
	Next
	
	Color 255,0,0
	Line gw/2-2,gh/2, gw/2+2,gh/2
	Line gw/2,gh/2-2, gw/2,gh/2+2

End Function

;Function translatePoints(tu#,tv#, R)
;
;	For p.point = Each point
;	
;		vdotx# = tu*p\u + tv*p\v
;		mag2v# = tu*tu + tv*tv
;		mag2x# = p\u*p\u + p\v*p\v
;		
;		;Stop
;		
;		tempu# = ( (1.+2*vdotx+mag2x)*tu + (1.-mag2v)*p\u ) / (1. + 2*vdotx + mag2v*mag2x)
;		tempv# = ( (1.+2*vdotx+mag2x)*tv + (1.-mag2v)*p\v ) / (1. + 2*vdotx + mag2v*mag2x)
;		
;		;Stop
;		
;		p\u = tempu
;		p\v = tempv
;		
;		p\x = p\u*R + GraphicsWidth()/2
;		p\y = p\v*R + GraphicsHeight()/2
;	
;	Next
;
;End Function

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

