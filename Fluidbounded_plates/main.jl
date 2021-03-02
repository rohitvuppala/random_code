### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 3b4135be-6682-11eb-1766-e772017ddc9b
begin
	using Plots
	using SpecialFunctions
	using LaTeXStrings
end

# ╔═╡ 88a17f66-667b-11eb-1d4d-57d100ff452e
md"
# Question:
Fluid bounded by two parallel plates extended to infinity. The lower wall is suddenly accelerated in the x-direction while the upper wall is stationary.

Governing Equation:

$\frac{\partial u}{\partial t} = \nu \frac{\partial^2 u}{\partial y^2}$

Initial condition:\
t = 0; u = 0     for 0 < y ≤ h \
t = 0; u = $U_0$ for y = 0 	

Boundary condition: \
t ≥ 0; u = $U_0$ for y = 0 \
t ≥ 0; u = 0 for y = h      

" 

# ╔═╡ 291e4152-667e-11eb-3d13-031ad9572916
# Code for solving the problem
begin
	ν = 0.000217
	h = 0.04
	u_b = 40.0
	n = 41
    δx = h/(n-1)
	
	function makegrid(h,n,δx)		
		for i in 1:n
			x[i] = 0.0 + (i-1)*δx
		end
		return x
	end
	
	function init(u,n)
		u      = zeros(n)
		u[1]   = u_b
		return u
	end
    
	function exact_sol(x,n,t_eval)
		@info(t_eval)
		uex = zeros(n)
		for i in 1:n
			η = x[i]/(2*sqrt(ν*t_eval))
			η1= h/(2*sqrt(ν*t_eval))
			uex[i] = u_b* (erfc(η) - erfc(2*η1 - η) + erfc(2*η1 + η)  - erfc(4*η1 - η) + erfc(4*η1 + η) )
		end
		return uex
	end
	
	#FTCS Time stepping
	function FTCS(δt,tmax,δx,u,n,ν,u_b)
		u_new = zeros(n)
		
	    #BC's
		u_new[1] = u_b
		u_new[n] = 0.0
		#Update
		t = 0.0
		nt= floor(Int,tmax/δt) 
		for i in 1:nt
			for j in 2:n-1
				u_new[j] = u[j] + ν*(δt/(δx)^2) * (u[j+1] - 2*u[j] + u[j-1]) 
			end
			t = t + δt 
			u = copy(u_new)
		end
		return u
	end
	
	#Dufort-Frankel Timestepping
	function DFF(δt,tmax,δx,u,n,ν,u_b)
		u_new = zeros(n)
		u_old = zeros(n)
		
	    #BC's
		u_new[1] = u_b
		u_new[n] = 0.0
		#Update
		t = 0.0
		nt= floor(Int,tmax/δt)
		d = ν*(δt/(δx)^2)
		
		u_old = copy(u)
		for i in 1:nt
			for j in 2:n-1
				u_new[j] = ((1.0-2.0*d)*u_old[j] + 2*d*(u[j+1] + u[j-1]))/(1.0+2.0*d) 
			end
			t = t + δt
			u_old = copy(u)
			u = copy(u_new)
		end
		return u
	end
	
	function TriD(A,B,C,D,n)
		#/// n is only "n-2" in global
		a = copy(A)
		b = copy(B)
		c = copy(C)
		d = copy(D)	
		u = zeros(n)
		for i in 2:n
			w = a[i]/b[i-1]
			b[i] = b[i] - w*c[i-1]
			d[i] = d[i] - w*d[i-1]
		end
		
		#Back solve
		u[n] = d[n]/b[n]
		for i in n-1:-1:1
			u[i] = (d[i] - c[i]*u[i+1])/b[i]
		end
		#@info(u)
		return u
	end
	
	#Laason implicit solve
	function LAimpl(δt,tmax,δx,u,n,ν,u_b)
		u_new = zeros(n)
		u_old = zeros(n)
		
	    #BC's
		u_new[1] = u_b
		u_new[n] = 0.0
		#Update
		t = 0.0
		nt= floor(Int,tmax/δt)
		d = ν*(δt/(δx)^2)
		
		A = zeros(n-2)
		B = zeros(n-2)
		C = zeros(n-2)
		D = zeros(n-2)
		for i in 1:nt
			A  .= -d
			B  .= 1+2*d
			C  .= -d
			D  .= copy(u[2:n-1])
			D[1]  = D[1]   + d*u_b
			D[n-2]= D[n-2] + d*0.0
            
			u_new[2:n-1] = TriD(A,B,C,D,n-2)
			
			u_new[1] = u_b
			u_new[n] = 0.0
			u = copy(u_new)
		end
		return u
	end
	
	#Crank-Nicolson implicit solve
	function CNimpl(δt,tmax,δx,u,n,ν,u_b)
		u_new = zeros(n)
		u_old = zeros(n)
		
	    #BC's
		u_new[1] = u_b
		u_new[n] = 0.0
		
		u[1] = u_b
		u[n] = 0.0
		#Update
		t = 0.0
		nt= floor(Int,tmax/δt)
		d = ν*(δt/(δx)^2)
		
		A = zeros(n-2)
		B = zeros(n-2)
		C = zeros(n-2)
		D = zeros(n-2)
		for i in 1:nt
			A  .= -0.5*d
			B  .= 1+1*d
			C  .= -0.5*d
			D  .= copy(u[2:n-1]) + 0.5*d*(u[3:n]-2.0*u[2:n-1]+u[1:n-2])
			
			D[1]  = D[1]   + 0.5*d*u_b
			D[n-2]= D[n-2] + 0.5*d*0.0
            
			u_new[2:n-1] = TriD(A,B,C,D,n-2)
			
			u_new[1] = u_b
			u_new[n] = 0.0
			u = copy(u_new)
		end
		return u
	end
	
	x = zeros(n)
	u = zeros(n)
	uex= zeros(n)
	
	#Run necessary stuff
	x = makegrid(h,n,δx)
	u = init(u,n)
	
	
	
end

# ╔═╡ 2e306e02-6719-11eb-2755-e364e8d201cc
md"# Velocity Profiles across schemes"

# ╔═╡ 270a7404-6713-11eb-2cf1-8b17120d4405
md"# Velocity Profiles for FTCS"

# ╔═╡ 1a86f16c-6713-11eb-2f1e-53e78d8e6a15
begin
	#Velocity Profiles for FTCS explicit method
	δt = 0.002
	tmax= [0.0,0.18,0.36,0.54,0.72,0.90,1.08]


	u_1 = FTCS(δt,tmax[1],δx,u,n,ν,u_b)
	u_2 = FTCS(δt,tmax[2],δx,u,n,ν,u_b)
	u_3 = FTCS(δt,tmax[3],δx,u,n,ν,u_b)
	u_4 = FTCS(δt,tmax[4],δx,u,n,ν,u_b)
	u_5 = FTCS(δt,tmax[5],δx,u,n,ν,u_b)
	u_6 = FTCS(δt,tmax[6],δx,u,n,ν,u_b)
	u_7 = FTCS(δt,tmax[7],δx,u,n,ν,u_b)

	gr()
	plot(title="Velocity Profile for dy = 0.001, dt = 0.002, d = 0.434")
	plot!(xlabel="u",ylabel="y")
	plot!(u_1,x,markershape = :ltriangle,label="t=0.0sec")
	plot!(u_2,x,markershape = :ltriangle,label="t=0.18sec")
	plot!(u_3,x,markershape = :ltriangle,label="t=0.36sec")
	plot!(u_4,x,markershape = :ltriangle,label="t=0.54sec")
	plot!(u_5,x,markershape = :ltriangle,label="t=0.72sec")
	plot!(u_6,x,markershape = :ltriangle,label="t=0.90sec")
	plot!(u_7,x,markershape = :ltriangle,label="t=1.08sec")
	
end

# ╔═╡ fb155b98-6715-11eb-0231-2bc9fde178fa
begin
	#Velocity Profiles for FTCS explicit method
	δt1 = 0.00232
	tmax1= [0.0,0.21,0.42,0.63,0.84,1.04,1.25]


	u1_1 = FTCS(δt1,tmax1[1],δx,u,n,ν,u_b)
	u1_2 = FTCS(δt1,tmax1[2],δx,u,n,ν,u_b)
	u1_3 = FTCS(δt1,tmax1[3],δx,u,n,ν,u_b)
	u1_4 = FTCS(δt1,tmax1[4],δx,u,n,ν,u_b)
	u1_5 = FTCS(δt1,tmax1[5],δx,u,n,ν,u_b)
	u1_6 = FTCS(δt1,tmax1[6],δx,u,n,ν,u_b)
	u1_7 = FTCS(δt1,tmax1[7],δx,u,n,ν,u_b)

	gr()
	plot(title="Velocity Profile for dy = 0.001, dt = 0.00232, d = 0.50344")
	plot!(xlabel="u",ylabel="y")
	plot!(u1_1,x,markershape = :ltriangle,label="t=0.0sec")
	plot!(u1_2,x,markershape = :ltriangle,label="t=0.21sec")
	plot!(u1_3,x,markershape = :ltriangle,label="t=0.42sec")
	plot!(u1_4,x,markershape = :ltriangle,label="t=0.63sec")
	plot!(u1_5,x,markershape = :ltriangle,label="t=0.84sec")
	plot!(u1_6,x,markershape = :ltriangle,label="t=1.04sec")
	plot!(u1_7,x,markershape = :ltriangle,label="t=1.25sec")
	
end

# ╔═╡ c3e4e836-6716-11eb-1fe8-4b01a8f19efe
md"# Velocity Profiles for Dufort-Frankel"

# ╔═╡ 800e1682-6716-11eb-2968-1f922c7f834c
begin
	#Velocity Profiles for DFF explicit method
	δt2 = 0.002
	tmax2= [0.0,0.18,0.36,0.54,0.72,0.90,1.08]


	u2_1 = DFF(δt2,tmax2[1],δx,u,n,ν,u_b)
	u2_2 = DFF(δt2,tmax2[2],δx,u,n,ν,u_b)
	u2_3 = DFF(δt2,tmax2[3],δx,u,n,ν,u_b)
	u2_4 = DFF(δt2,tmax2[4],δx,u,n,ν,u_b)
	u2_5 = DFF(δt2,tmax2[5],δx,u,n,ν,u_b)
	u2_6 = DFF(δt2,tmax2[6],δx,u,n,ν,u_b)
	u2_7 = DFF(δt2,tmax2[7],δx,u,n,ν,u_b)

	gr()
	plot(title="Velocity Profile for dy = 0.001, dt = 0.002, d = 0.434")
	plot!(xlabel="u",ylabel="y")
	plot!(u2_1,x,markershape = :ltriangle,label="t=0.0sec")
	plot!(u2_2,x,markershape = :ltriangle,label="t=0.18sec")
	plot!(u2_3,x,markershape = :ltriangle,label="t=0.36sec")
	plot!(u2_4,x,markershape = :ltriangle,label="t=0.54sec")
	plot!(u2_5,x,markershape = :ltriangle,label="t=0.72sec")
	plot!(u2_6,x,markershape = :ltriangle,label="t=0.90sec")
	plot!(u2_7,x,markershape = :ltriangle,label="t=1.08sec")
	
end

# ╔═╡ 39b664cc-6717-11eb-11b1-b74350c4341d
begin
	#Velocity Profiles for FTCS explicit method
	δt3 = 0.003
	tmax3= [0.0,0.21,0.42,0.63,0.84,1.04,1.25]


	u3_1 = DFF(δt3,tmax3[1],δx,u,n,ν,u_b)
	u3_2 = DFF(δt3,tmax3[2],δx,u,n,ν,u_b)
	u3_3 = DFF(δt3,tmax3[3],δx,u,n,ν,u_b)
	u3_4 = DFF(δt3,tmax3[4],δx,u,n,ν,u_b)
	u3_5 = DFF(δt3,tmax3[5],δx,u,n,ν,u_b)
	u3_6 = DFF(δt3,tmax3[6],δx,u,n,ν,u_b)
	u3_7 = DFF(δt3,tmax3[7],δx,u,n,ν,u_b)

	gr()
	plot(title="Velocity Profile for dy = 0.001, dt = 0.003, d = 0.651")
	plot!(xlabel="u",ylabel="y")
	plot!(u3_1,x,markershape = :ltriangle,label="t=0.0sec")
	plot!(u3_2,x,markershape = :ltriangle,label="t=0.21sec")
	plot!(u3_3,x,markershape = :ltriangle,label="t=0.42sec")
	plot!(u3_4,x,markershape = :ltriangle,label="t=0.63sec")
	plot!(u3_5,x,markershape = :ltriangle,label="t=0.84sec")
	plot!(u3_6,x,markershape = :ltriangle,label="t=1.04sec")
	plot!(u3_7,x,markershape = :ltriangle,label="t=1.25sec")
	
end

# ╔═╡ 99b7da0e-6717-11eb-061d-ef8e525fa7a6
md"# Velocity Profiles for Laason Implicit"

# ╔═╡ c2e5d336-6717-11eb-04f4-d9506627a208
begin
	#Velocity Profiles for LAimplicit method
	δt4 = 0.002
	tmax4= [0.0,0.18,0.36,0.54,0.72,0.90,1.08]


	u4_1 = LAimpl(δt4,tmax4[1],δx,u,n,ν,u_b)
	u4_2 = LAimpl(δt4,tmax4[2],δx,u,n,ν,u_b)
	u4_3 = LAimpl(δt4,tmax4[3],δx,u,n,ν,u_b)
	u4_4 = LAimpl(δt4,tmax4[4],δx,u,n,ν,u_b)
	u4_5 = LAimpl(δt4,tmax4[5],δx,u,n,ν,u_b)
	u4_6 = LAimpl(δt4,tmax4[6],δx,u,n,ν,u_b)
	u4_7 = LAimpl(δt4,tmax4[7],δx,u,n,ν,u_b)

	gr()
	plot(title="Velocity Profile for dy = 0.001, dt = 0.002, d = 0.434")
	plot!(xlabel="u",ylabel="y")
	plot!(u4_1,x,markershape = :ltriangle,label="t=0.0sec")
	plot!(u4_2,x,markershape = :ltriangle,label="t=0.18sec")
	plot!(u4_3,x,markershape = :ltriangle,label="t=0.36sec")
	plot!(u4_4,x,markershape = :ltriangle,label="t=0.54sec")
	plot!(u4_5,x,markershape = :ltriangle,label="t=0.72sec")
	plot!(u4_6,x,markershape = :ltriangle,label="t=0.90sec")
	plot!(u4_7,x,markershape = :ltriangle,label="t=1.08sec")
		plot!(exact_sol(x,n,0.18),x)
end

# ╔═╡ 0ce28b94-6718-11eb-1793-c5f06bfe5686
begin
	#Velocity Profiles for LaasonImplicit method
	δt5 = 0.01
	tmax5= [0.0,0.21,0.42,0.63,0.84,1.04,1.25]


	u5_1 = LAimpl(δt5,tmax5[1],δx,u,n,ν,u_b)
	u5_2 = LAimpl(δt5,tmax5[2],δx,u,n,ν,u_b)
	u5_3 = LAimpl(δt5,tmax5[3],δx,u,n,ν,u_b)
	u5_4 = LAimpl(δt5,tmax5[4],δx,u,n,ν,u_b)
	u5_5 = LAimpl(δt5,tmax5[5],δx,u,n,ν,u_b)
	u5_6 = LAimpl(δt5,tmax5[6],δx,u,n,ν,u_b)
	u5_7 = LAimpl(δt5,tmax5[7],δx,u,n,ν,u_b)

	gr()
	plot(title="Velocity Profile for dy = 0.001, dt = 0.01, d = 2.17")
	plot!(xlabel="u",ylabel="y")
	plot!(u3_1,x,markershape = :ltriangle,label="t=0.0sec")
	plot!(u3_2,x,markershape = :ltriangle,label="t=0.21sec")
	plot!(u3_3,x,markershape = :ltriangle,label="t=0.42sec")
	plot!(u3_4,x,markershape = :ltriangle,label="t=0.63sec")
	plot!(u3_5,x,markershape = :ltriangle,label="t=0.84sec")
	plot!(u3_6,x,markershape = :ltriangle,label="t=1.04sec")
	plot!(u3_7,x,markershape = :ltriangle,label="t=1.25sec")

end

# ╔═╡ 1eb83680-6719-11eb-16c6-89de2999c61c
md"# Velocity Profiles for Crank-Nicolson Implicit"

# ╔═╡ 1d27d028-6719-11eb-3ac6-7147ea8109ef
begin
	#Velocity Profiles for CNimplicit method
	δt6 = 0.002
	tmax6= [0.0,0.18,0.36,0.54,0.72,0.90,1.08]


	u6_1 = CNimpl(δt6,tmax6[1],δx,u,n,ν,u_b)
	u6_2 = CNimpl(δt6,tmax6[2],δx,u,n,ν,u_b)
	u6_3 = CNimpl(δt6,tmax6[3],δx,u,n,ν,u_b)
	u6_4 = CNimpl(δt6,tmax6[4],δx,u,n,ν,u_b)
	u6_5 = CNimpl(δt6,tmax6[5],δx,u,n,ν,u_b)
	u6_6 = CNimpl(δt6,tmax6[6],δx,u,n,ν,u_b)
	u6_7 = CNimpl(δt6,tmax6[7],δx,u,n,ν,u_b)

	gr()
	plot(title="Velocity Profile for dy = 0.001, dt = 0.002, d = 0.434")
	plot!(xlabel="u",ylabel="y")
	plot!(u6_1,x,markershape = :ltriangle,label="t=0.0sec")
	plot!(u6_2,x,markershape = :ltriangle,label="t=0.18sec")
	plot!(u6_3,x,markershape = :ltriangle,label="t=0.36sec")
	plot!(u6_4,x,markershape = :ltriangle,label="t=0.54sec")
	plot!(u6_5,x,markershape = :ltriangle,label="t=0.72sec")
	plot!(u6_6,x,markershape = :ltriangle,label="t=0.90sec")
	plot!(u6_7,x,markershape = :ltriangle,label="t=1.08sec")
    plot!(exact_sol(x,n,0.18),x)
end
	

# ╔═╡ 04589ef0-671a-11eb-3b04-4dc8b7a0fd43
begin
	#Velocity Profiles for CNImplicit method
	δt7 = 0.01
	tmax7= [0.0,0.21,0.42,0.63,0.84,1.04,1.25]


	u7_1 = CNimpl(δt7,tmax5[1],δx,u,n,ν,u_b)
	u7_2 = CNimpl(δt7,tmax5[2],δx,u,n,ν,u_b)
	u7_3 = CNimpl(δt7,tmax5[3],δx,u,n,ν,u_b)
	u7_4 = CNimpl(δt7,tmax5[4],δx,u,n,ν,u_b)
	u7_5 = CNimpl(δt7,tmax5[5],δx,u,n,ν,u_b)
	u7_6 = CNimpl(δt7,tmax5[6],δx,u,n,ν,u_b)
	u7_7 = CNimpl(δt7,tmax5[7],δx,u,n,ν,u_b)

	gr()
	plot(title="Velocity Profile for dy = 0.001, dt = 0.01, d = 2.17")
	plot!(xlabel="u",ylabel="y")
	plot!(u7_1,x,markershape = :ltriangle,label="t=0.0sec")
	plot!(u7_2,x,markershape = :ltriangle,label="t=0.21sec")
	plot!(u7_3,x,markershape = :ltriangle,label="t=0.42sec")
	plot!(u7_4,x,markershape = :ltriangle,label="t=0.63sec")
	plot!(u7_5,x,markershape = :ltriangle,label="t=0.84sec")
	plot!(u7_6,x,markershape = :ltriangle,label="t=1.04sec")
	plot!(u7_7,x,markershape = :ltriangle,label="t=1.25sec")
	plot!(exact_sol(x,n,0.21),x)
	
end

# ╔═╡ 08bc0ee8-671d-11eb-1e49-e392f826eba5
md"# Error Distribution for t = 0.18 sec, 1.08sec"

# ╔═╡ bb3fd91e-671d-11eb-1c84-d518b84da067
#Calculating Error
function err_calc(uexact,ucalc,n)
	err = zeros(n)
	for i in 1:n
		err[i] = (uexact[i] - ucalc[i]) 
	end
	return err
end

# ╔═╡ 24c371dc-671d-11eb-2611-c39d01c2ecf2
begin
	u_ex = exact_sol(x,n,0.18)
	u_ftcs = FTCS(0.002,0.18,δx,u,n,ν,u_b)
	u_dff  = DFF(0.002,0.18,δx,u,n,ν,u_b)
	u_la   = LAimpl(0.002,0.18,δx,u,n,ν,u_b)
	u_cn   = CNimpl(0.002,0.18,δx,u,n,ν,u_b)
	
	
	ftcs_err = err_calc(u_ex,u_ftcs,n)
	dff_err  = err_calc(u_ex,u_dff,n)
	la_err  = err_calc(u_ex,u_la,n)
	cn_err  = err_calc(u_ex,u_cn,n)
	
	plot(title="Error comparison for various schmes dt = 0.002 at t=0.18")
	plot!(xlabel="Error",ylabel="y",xlim=(-0.11,0.11))
	plot!(ftcs_err,x,markershape = :ltriangle,label="FTCS")
	plot!(dff_err,x,markershape = :ltriangle,label="Dufort-Frankel")
	plot!(la_err,x,markershape = :ltriangle,label="Laasonen-implicit")
	plot!(cn_err,x,markershape = :ltriangle,label="Crank-Nicolson-implicit")
	
end

# ╔═╡ 615146cc-672c-11eb-097d-41f289d3f896
begin
	u1_ex = exact_sol(x,n,1.08)
	u1_ftcs = FTCS(0.002,1.08,δx,u,n,ν,u_b)
	u1_dff  = DFF(0.002,1.08,δx,u,n,ν,u_b)
	u1_la   = LAimpl(0.002,1.08,δx,u,n,ν,u_b)
	u1_cn   = CNimpl(0.002,1.08,δx,u,n,ν,u_b)
	
	
	ftcs1_err = err_calc(u1_ex,u1_ftcs,n)
	dff1_err  = err_calc(u1_ex,u1_dff,n)
	la1_err  = err_calc(u1_ex,u1_la,n)
	cn1_err  = err_calc(u1_ex,u1_cn,n)
	
	plot(title="Error comparison for various schmes dt = 0.002 and t=1.08")
	plot!(xlabel="Error",ylabel="y",xlim=(-0.05,0.05))
	plot!(ftcs1_err,x,markershape = :ltriangle,label="FTCS")
	plot!(dff1_err,x,markershape = :ltriangle,label="Dufort-Frankel")
	plot!(la1_err,x,markershape = :ltriangle,label="Laasonen-implicit")
	plot!(cn1_err,x,markershape = :ltriangle,label="Crank-Nicolson-implicit")
	
end

# ╔═╡ dc70983e-672d-11eb-101b-47d6c05f02e2
md"# Error Distribution for different time-steps at t = 1.0sec for Lassonen"

# ╔═╡ 251d0ac8-672d-11eb-1d64-8daa66f59ffa
begin
	u2_ex = exact_sol(x,n,1.0)
	ula_1 = LAimpl(0.005,1.0,δx,u,n,ν,u_b)
	ula_2 = LAimpl(0.01,1.0,δx,u,n,ν,u_b)
	ula_3 = LAimpl(0.1,1.0,δx,u,n,ν,u_b)
	ula_4 = LAimpl(0.2,1.0,δx,u,n,ν,u_b)
	
	
	ula1_err  = err_calc(u2_ex,ula_1,n)
	ula2_err  = err_calc(u2_ex,ula_2,n)
	ula3_err  = err_calc(u2_ex,ula_3,n)
	ula4_err  = err_calc(u2_ex,ula_4,n)
	
	plot(title="Error comparison for different timesteps at t=1.0 
		for Lassonen")
	plot!(xlabel="Error",ylabel="y")
	plot!(ula1_err,x,markershape = :ltriangle,label="dt=0.005")
	plot!(ula2_err,x,markershape = :ltriangle,label="dt=0.01")
	plot!(ula3_err,x,markershape = :ltriangle,label="dt=0.1")
	plot!(ula4_err,x,markershape = :ltriangle,label="dt=0.2")
	
end

# ╔═╡ Cell order:
# ╟─88a17f66-667b-11eb-1d4d-57d100ff452e
# ╠═3b4135be-6682-11eb-1766-e772017ddc9b
# ╠═291e4152-667e-11eb-3d13-031ad9572916
# ╟─2e306e02-6719-11eb-2755-e364e8d201cc
# ╟─270a7404-6713-11eb-2cf1-8b17120d4405
# ╟─1a86f16c-6713-11eb-2f1e-53e78d8e6a15
# ╟─fb155b98-6715-11eb-0231-2bc9fde178fa
# ╟─c3e4e836-6716-11eb-1fe8-4b01a8f19efe
# ╟─800e1682-6716-11eb-2968-1f922c7f834c
# ╟─39b664cc-6717-11eb-11b1-b74350c4341d
# ╟─99b7da0e-6717-11eb-061d-ef8e525fa7a6
# ╟─c2e5d336-6717-11eb-04f4-d9506627a208
# ╟─0ce28b94-6718-11eb-1793-c5f06bfe5686
# ╟─1eb83680-6719-11eb-16c6-89de2999c61c
# ╟─1d27d028-6719-11eb-3ac6-7147ea8109ef
# ╟─04589ef0-671a-11eb-3b04-4dc8b7a0fd43
# ╟─08bc0ee8-671d-11eb-1e49-e392f826eba5
# ╠═bb3fd91e-671d-11eb-1c84-d518b84da067
# ╟─24c371dc-671d-11eb-2611-c39d01c2ecf2
# ╟─615146cc-672c-11eb-097d-41f289d3f896
# ╟─dc70983e-672d-11eb-101b-47d6c05f02e2
# ╟─251d0ac8-672d-11eb-1d64-8daa66f59ffa
