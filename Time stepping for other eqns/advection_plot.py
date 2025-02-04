import numpy as np
import matplotlib.pyplot as plt

# Domain and discretization parameters
N = 1000  # Number of grid points
L = 2 * np.pi  # Domain length
dx = L / N  # Grid spacing
mu = 2 * np.pi / 9990  # Viscosity coefficient
T = 70.0  # Extended simulation time for steady state
dt = 0.001  # Smaller time step for better resolution
omega = 1.0  # Frequency of forcing term
c = 1.0  # Advection speed

# Discretized spatial domain
x = np.linspace(0, L, N, endpoint=False)  # Periodic domain

# First derivative matrix (central difference)
A = np.zeros((N, N))
for i in range(N):
    A[i, (i + 1) % N] = 0.5 / dx
    A[i, (i - 1) % N] = -0.5 / dx

# Second derivative matrix (central difference)
B = np.zeros((N, N))
for i in range(N):
    B[i, (i + 1) % N] = 1.0 / dx**2
    B[i, i] = -2.0 / dx**2
    B[i, (i - 1) % N] = 1.0 / dx**2

# Initial condition
u = 1.1*np.cos(x-1)

# Time stepping loop
t = 0.0
times = [t]
sol = [u.copy()]

while t < T:
    f = np.sin(omega * t)  # Forcing term
    du_dt = c * (A @ u) + mu * (B @ u) + f  # Corrected advection term
    u = u + dt * du_dt  # Forward Euler step

    t += dt
    times.append(t)
    sol.append(u.copy())

# Convert lists to numpy arrays for plotting and subtract 1 from each solution element
times = np.array(times)
Z = np.array(sol) - 1
...
# # Create contour plot
# plt.figure(figsize=(10, 6))
# # Set vmin and vmax to force the colorbar range from -2 to 2
# cf = plt.contourf(times, x, Z.T, levels=20, cmap='inferno', vmin=-2, vmax=2)
# # Add contour lines
# contours = plt.contour(times, x, Z.T, levels=20, colors='black', linewidths=0.5)
# # Attach colorbar with customized tick labels from 2 down to -2 
# cbar = plt.colorbar(cf, label='u(x,t)', ticks=[2.0, 1.5, 1.0, 0.5, 0.0, -0.5, -1.0, -1.5, -2.0])
# plt.xlabel('Time')
# plt.ylabel('Space')
# plt.title('Contour plot of u(x,t), nx = 1000, Central Scheme using Time stepping ')
# plt.show()


t_cut = 1.0  # for example, ignore data before t = 10.0 seconds
cut_idx = np.where(times >= t_cut)[0][0]

# Extract the steady state part at x=pi (time series), ignoring early transients.
u_pi_steady = Z[cut_idx:, N//2]

# Also extract the corresponding times (if needed)
times_steady = times[cut_idx:]

# Compute the FFT on the steady state portion using one-sided FFT
U_fft = np.fft.rfft(u_pi_steady) * 2 / len(u_pi_steady)
# Get frequencies (radians/s): np.fft.rfftfreq returns frequencies in Hz so multiply by 2π
freq = np.fft.rfftfreq(len(u_pi_steady), d=dt) * 2 * np.pi

plt.figure(figsize=(10, 6))
plt.plot(freq, np.abs(U_fft))
plt.xlim(0, 6)
plt.xlabel('Frequency (rad/s)')
plt.ylabel('Amplitude')
plt.title('FFT of u(x=pi, t) ')
plt.show()
