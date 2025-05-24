using Plots

"""
Решение одномерного уравнения переноса
с использованием различных численных методов.
"""

# Параметры задачи
const c = 1.0  # скорость переноса
const L = 10.0  # длина области
const T = 1.0   # время моделирования
const dx = 0.01  # шаг по пространству
const dt = 0.005 # шаг по времени
const x = 0:dx:L
const N = length(x)
const Nt = Int(round(T/dt))

# Начальное условие (гаусс)
function initial_condition(x)
    exp(-(x - L/4)^2 / 0.5)
end

# Метод Upwind
function upwind(u, c, dt, dx)
    u_new = similar(u)
    if c > 0
        for i in 2:N
            u_new[i] = u[i] - c*dt/dx * (u[i] - u[i-1])
        end
        u_new[1] = u_new[2]  
    else
        for i in 1:N-1
            u_new[i] = u[i] - c*dt/dx * (u[i+1] - u[i])
        end
        u_new[N] = u_new[N-1]  
    end
    return u_new
end

# Метод Лакса-Вендроффа
function lax_wendroff(u, c, dt, dx)
    u_new = similar(u)
    for i in 2:N-1
        u_new[i] = u[i] - c*dt/(2*dx)*(u[i+1] - u[i-1]) + 
                   c^2*dt^2/(2*dx^2)*(u[i+1] - 2*u[i] + u[i-1])
    end
    u_new[1] = u_new[2]
    u_new[N] = u_new[N-1]
    return u_new
end

# Аналитическое решение
function exact_solution(x, t)
    x_shift = mod.(x .- c*t, L)
    return initial_condition.(x_shift)
end

function compare_methods()
    u0 = initial_condition.(x)
    u_upwind = copy(u0)
    u_lw = copy(u0)
    u_exact = copy(u0)
    
    for i in 1:Nt
        u_upwind = upwind(u_upwind, c, dt, dx)
        u_lw = lax_wendroff(u_lw, c, dt, dx)
    end
    u_exact = exact_solution(x, T)
    
    # Решение
    plot(x, u_upwind, label="Upwind", linewidth=2, xlims=(0, L), ylims=(0, 1.2))
    plot!(x, u_lw, label="Лакса-Вендроффа", linewidth=2)
    plot!(x, u_exact, label="Точное решение", linewidth=2, linestyle=:dash)
    plot!(x, u0, label="Начальное условие", linewidth=2, linestyle=:dash)
    title!("Решение t=$(round(T, digits=2))")
    xlabel!("x")
    ylabel!("u(x,t)")
    png(raw"C:\Users\Пользователь\Desktop\distr\решение")

    # Разность
    delta_u_upwind = u_upwind.-u_exact
    delta_u_lw = u_lw.-u_exact
    plot(x, delta_u_upwind, label="Upwind", linewidth=2, xlims=(0, L))
    plot!(x, delta_u_lw, label="Лакса-Вендроффа", linewidth=2)
    title!("Решение t=$(round(T, digits=2))")
    xlabel!("x")
    ylabel!("Δu(x,t)")
    png(raw"C:\Users\Пользователь\Desktop\distr\разность")
end

compare_methods()