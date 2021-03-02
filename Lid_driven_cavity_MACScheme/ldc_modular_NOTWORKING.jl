using PyPlot
using OffsetArrays
"""
Program to solve 2D lid driven cavity steady flow problem
"""
# Initial settings
xmin = 0.0
xmax = 1.0
ymin = 0.0
ymax = 1.0

global u_lid = 1.0
global Re = 100.0

#Total points in each direction
nx = 101
ny = 101

global dx = (xmax - xmin) / (nx - 1)
global dy = (ymax - ymin) / (ny - 1)

cfl = 0.2

#
# We are using Staggered grid
#

#u,v needs ghosts
uloc = OffsetArray(zeros(nx, ny + 1, 2), 0:nx-1, 0:ny, 1:2)
vloc = OffsetArray(zeros(nx + 1, ny, 2), 0:nx, 0:ny-1, 1:2)

# p needs ghost cells
global ploc = OffsetArray(zeros(nx + 1, ny + 1, 2), 0:nx, 0:ny, 1:2)


# Get the locations for debugging purposes
for i = 0:nx-1
    for j = 1:ny-1
        uloc[i, j, 1] = xmin + (i) * dx           #x loc
        uloc[i, j, 2] = ymin + dy / 2 + (j - 1) * dy  #y loc
    end
end

for i = 1:nx-1
    for j = 0:ny-1
        vloc[i, j, 1] = xmin + dx / 2 + (i - 1) * dx
        vloc[i, j, 2] = ymin + (j) * dy
    end
end


# Allocate arrays for storage
u = OffsetArray(zeros(nx, ny + 1), 0:nx-1, 0:ny)
v = OffsetArray(zeros(nx + 1, ny), 0:nx, 0:ny-1)
p = OffsetArray(zeros(nx + 1, ny + 1), 0:nx, 0:ny)

uα = OffsetArray(zeros(nx, ny + 1), 0:nx-1, 0:ny)
vα = OffsetArray(zeros(nx + 1, ny), 0:nx, 0:ny-1)

#Initialize arrays
u .= 0.0
v .= 0.0
p .= 0.0

uα .= 0.0
vα .= 0.0

#Function to calculate dt
function calc_dt(u, v, dx, dy, cfl)
    umax = maximum(u)
    vmax = maximum(v)
    dt = minimum([cfl * dx / umax, cfl * dy / vmax])
    return dt
end

#Function for expbc
function expbc!(u, v, p)
    #Lower BC
    u[:, 0] = -u[:, 1]
    v[:, 0] .= 0.0
    p[:, 0] = p[:, 1]

    #Upper BC
    u[:, end] = 2 * u_lid .- u[:, end-1]
    v[:, end] .= 0.0
    p[:, end] = p[:, end-1]

    #Left BC
    u[0, :] .= 0.0
    v[0, :] = -v[1, :]
    p[0, :] = p[1, :]

    #Right BC
    u[end, :] .= 0.0
    v[end, :] = -v[end-1, :]
    p[end, :] = p[end-1, :]

    return nothing
end

#Function find uα and vα
function predictor!(u, v, dt, uα, vα)
    #calculate u_approx
    for i = 1:nx-2
        for j = 1:ny-1
            u1 = -0.5 * (u[i+1, j]^2 - u[i-1, j]^2) / (dx)

            u21 = (u[i, j+1] + u[i, j]) * (v[i, j] + v[i+1, j])
            u22 = (u[i, j] + u[i, j-1]) * (v[i, j-1] + v[i+1, j-1])
            u2 = -0.25 * (u21 - u22) / dy

            u31 = (u[i+1, j] - 2 * u[i, j] + u[i-1, j]) / (dx^2)
            u32 = (u[i, j+1] - 2 * u[i, j] + u[i, j+1]) / (dy^2)
            u3 = 1 / Re * (u31 + u32)

            uα[i, j] = u[i, j] + dt * (u1 + u2 + u3)
        end
    end

    #calculate v_approx
    for i = 1:nx-1
        for j = 1:ny-2
            v11 = (u[i, j+1] + u[i, j]) * (v[i, j] + v[i+1, j])
            v12 = (u[i-1, j+1] + u[i-1, j]) * (v[i, j] + v[i-1, j])
            v1 = -0.25 * (v11 - v12) / dx

            v2 = -0.5 * (v[i, j+1^2] - v[i, j-1]^2) / dy

            v31 = (v[i, j+1] - 2 * v[i, j] + v[i, j-1]) / (dy^2)
            v32 = (v[i+1, j] - 2 * v[i, j] + v[i-1, j]) / (dx^2)
            v3 = 1 / Re * (v31 + v32)

            vα[i, j] = v[i, j] + dt * (v1 + v2 + v3)
        end
    end

    #Function for expbc
    #Lower BC
    uα[:, 0] = -uα[:, 1]
    vα[:, 0] .= 0.0

    #Upper BC
    uα[:, end] = 2 * u_lid .- uα[:, end-1]
    vα[:, end] .= 0.0

    #Left BC
    uα[0, :] .= 0.0
    vα[0, :] = -vα[1, :]

    #Right BC
    uα[end, :] .= 0.0
    vα[end, :] = -vα[end-1, :]

    return nothing
end


#Function for pressure solver- Gauss Seidel
function psolve!(uα, vα, dt, p)
    tol = 1.0E-4
    nmax = 10
    for n = 1:nmax
        for i = 1:nx-1
            for j = 1:ny-1
                rhs =
                    0.5 / dt * ((uα[i, j] - uα[i-1, j]) / dx + (vα[i, j] - vα[i, j-1]) / dy)
                p[i, j] =
                    -0.5 / (1 / dx^2 + 1 / dy^2) * (
                        rhs - 1 / dx^2 * (p[i+1, j] + p[i-1, j]) -
                        1 / dy^2 * (p[i, j+1] + p[i, j-1])
                    )
            end
        end

    end

    return nothing
end

#Function for corrector
function corrector!(u, v, uα, vα, p, dt)

    #calculate u_approx
    for i = 1:nx-2
        for j = 1:ny-1
            u[i, j] = uα[i, j] - dt / dx * (p[i+1, j] - p[i, j])
        end
    end

    #calculate v_approx
    for i = 1:nx-1
        for j = 1:ny-2
            v[i, j] = vα[i, j] - dt / dy * (p[i, j+1] - p[i, j])
        end
    end

    return nothing
end

function main(nmax)
    for n = 1:nmax
        println(n)
        expbc!(u, v, p)
        dt = calc_dt(u, v, dx, dy, cfl)
        predictor!(u, v, dt, uα, vα)
        expbc!(uα, vα, p)
        psolve!(uα, vα, dt, p)
        expbc!(uα, vα, p)
        corrector!(u, v, uα, vα, p, dt)
    end
    return nothing
end
#%%
@btime main(10)
contourf(transpose(u), levels = 10)
