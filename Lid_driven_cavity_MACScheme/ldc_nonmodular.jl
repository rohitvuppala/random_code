using PyPlot
using OffsetArrays
using Printf
using DelimitedFiles
"""
	Program to solve 2D lid driven cavity steady flow problem
"""
function ldcav(Re = 100.0, nx = 100, ny = 100)

    #nx = 100
    #ny = 100

    dx = 1.0 / nx
    dy = 1.0 / ny

    dtc = min(dx, dy)
    dtv = 0.25 * Re * min(dx^2, dy^2)
    sigma = 0.5
    dt = sigma * min(dtc, dtv)

    eps = 1.0E-6
    ptol = 1.0E-2
    #%%
    #Set initial values
    u = OffsetArray(zeros(nx + 1, ny + 2), 0:nx, 0:ny+1)
    v = OffsetArray(zeros(nx + 2, ny + 1), 0:nx+1, 0:ny)
    p = OffsetArray(zeros(nx + 2, ny + 2), 0:nx+1, 0:ny+1)

    u0 = zeros(nx, ny)
    v0 = zeros(nx, ny)
    p0 = zeros(nx, ny)

    kmax = 100000
    iimax = 100

    resu = Float64[]
    resv = Float64[]
    resp = Float64[]

    for k = 1:kmax
        # copy for future residual plotting
        u0 = copy(u)
        v0 = copy(v)
        p0 = copy(p)

        #BCs
        u[0:nx, 0] = -copy(u[0:nx, 1])
        u[0:nx, ny+1] = -copy(u[0:nx, ny]) .+ 2.0

        v[0, 0:ny] = copy(v[1, 0:ny])
        v[nx+1, 0:ny] = -copy(v[nx, 0:ny])

        #Compute the predicted velocities
        for i = 1:nx-1
            for j = 1:ny
                c3 =
                    1 / Re * (
                        (u[i+1, j] - 2.0 * u[i, j] + u[i-1, j]) / (dx^2) +
                        (u[i, j+1] - 2.0 * u[i, j] + u[i, j-1]) / (dy^2)
                    )
                c1 =
                    (0.25 / dx) *
                    ((u[i+1, j] + u[i, j])^2 - (u[i, j] + u[i-1, j])^2)
                c2 =
                    (0.25 / dy) * (
                        (u[i, j+1] + u[i, j]) * (v[i+1, j] + v[i, j]) -
                        (u[i, j] + u[i, j-1]) * (v[i+1, j-1] + v[i, j-1])
                    )

                u[i, j] = copy(u[i, j]) + dt * (-c1 - c2 + c3)
            end
        end

        for i = 1:nx
            for j = 1:ny-1
                c3 =
                    (1.0 / Re) * (
                        (v[i+1, j] - 2.0 * v[i, j] + v[i-1, j]) / (dx^2) +
                        (v[i, j+1] - 2.0 * v[i, j] + v[i, j-1]) / (dy * dy)
                    )
                c2 =
                    (0.25 / dy) *
                    ((v[i, j+1] + v[i, j])^2 - (v[i, j] + v[i, j-1])^2)
                c1 =
                    (0.25 / dx) * (
                        (u[i, j+1] + u[i, j]) * (v[i+1, j] + v[i, j]) -
                        (u[i-1, j+1] + u[i-1, j]) * (v[i, j] + v[i-1, j])
                    )

                v[i, j] = copy(v[i, j]) + dt * (-c1 - c2 + c3)
            end
        end

        #BCs
        u[0:nx, 0] = -copy(u[0:nx, 1])
        u[0:nx, ny+1] = -copy(u[0:nx, ny]) .+ 2.0

        v[0, 0:ny] = copy(v[1, 0:ny])
        v[nx+1, 0:ny] = -copy(v[nx, 0:ny])

        #B.Cs for pressure (Neumann condition)
        p[1:nx, 0] = copy(p[1:nx, 1])
        p[1:nx, ny+1] = copy(p[1:nx, ny])

        p[0, 0:ny+1] = copy(p[1, 0:ny+1])
        p[nx+1, 0:ny+1] = copy(p[nx, 0:ny+1])

        #Pressure poisson solver
        a = -2.0 / (dx^2) - 2.0 / (dy^2)
        omega = 1.0
        for i1 = 1:iimax
            pold = copy(p)
            pnorm = 0.0
            for i = 1:nx
                for j = 1:ny
                    f =
                        (
                            (u[i, j] - u[i-1, j]) / dx +
                            (v[i, j] - v[i, j-1]) / dy
                        ) / dt
                    r = (
                        f -
                        (p[i+1, j] - 2.0 * p[i, j] + p[i-1, j]) / (dx * dx) -
                        (p[i, j+1] - 2.0 * p[i, j] + p[i, j-1]) / (dy * dy)
                    )

                    p[i, j] = copy(p[i, j]) + omega * r / a
                    pnorm = copy(pnorm) + (p[i, j] - pold[i, j])^2
                end
            end
            pnorm = sqrt(pnorm / (nx * ny))
            if pnorm <= ptol
                break
            end
        end

        #Corrector
        for i = 1:nx-1
            for j = 1:ny
                u[i, j] = copy(u[i, j]) - dt * (p[i+1, j] - p[i, j]) / dx
            end
        end

        for i = 1:nx
            for j = 1:ny-1
                v[i, j] = copy(v[i, j]) - dt * (p[i, j+1] - p[i, j]) / dy
            end
        end

        #Compute residuals
        ru = 0.0
        rv = 0.0
        rp = 0.0
        for i = 1:nx
            for j = 1:ny
                ru = copy(ru) + (u[i, j] - u0[i, j])^2
                rv = copy(rv) + (v[i, j] - v0[i, j])^2
                rp = copy(rp) + (p[i, j] - p0[i, j])^2
            end
        end
        ru = sqrt(ru / (nx * ny))
        rv = sqrt(rv / (nx * ny))
        rp = sqrt(rp / (nx * ny))

        #println(k,' ',ru,' ',rv,' ',rp)
        @printf("%6d %.3E %.3E %.3E \n", k, ru, rv, rp)

        #Make array of residuals
        push!(resu, ru)
        push!(resv, rv)
        push!(resp, rp)

        if ru <= eps && rv <= eps && rp <= eps
            break
        end
    end
    return u, v, p, resu, resv, resp
end

#%%
u, v, p, resu, resv, resp = ldcav(100.0, 50, 50)
#%%
#Writing the fields
#open("compareu.txt","w") do io
#    writedlm(io,u[])
#end
#%%
#Plotting
if isdir("figs")
else
    mkdir("figs")
end
figure()
title("U-vel")
xlabel("x")
ylabel("y")
contourf(transpose(u), levels = 30, cmap = "jet")
colorbar()
savefig("figs/uvel.jpg")

figure()
title("V-vel")
xlabel("x")
ylabel("y")
contourf(transpose(v), levels = 30, cmap = "jet")
colorbar()
savefig("figs/vvel.jpg")

figure()
title("Pressure")
xlabel("x")
ylabel("y")
contourf(transpose(p), levels = 30, cmap = "jet")
colorbar()
savefig("figs/pres.jpg")

figure()
semilogx()
semilogy()
title("Residuals")
xlabel("Iterations")
ylabel("Residual")
plot(resu, label = "res u ")
plot(resv, label = "res v")
plot(resp, label = "res p ")
legend()
savefig("figs/res.jpg")
