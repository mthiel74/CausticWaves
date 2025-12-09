#!/usr/bin/env python3
"""
Pool Caustics Simulator (Geometric Optics)
==========================================

This script simulates the optical phenomenon of water caustics—the bright patterns
formed on the bottom of a pool by sunlight refracting through surface ripples.

Methodology:
1.  **Surface Generation**: We use a spectral method (Tessendorf waves) to generate
    realistic, band-limited water ripples. This is preferred over PDE-based methods
    (like Shallow Water Equations) for this visual application because it allows
    for precise control over the wave spectrum and dispersion relations, resulting
    in more aesthetically pleasing "pool-like" ripples.
2.  **Geometric Optics**: We trace light rays from the sun, refract them at the
    water surface using Snell's Law, and intersect them with the pool floor.
3.  **Intensity accumulation**: We use a 2D histogram to accumulate the light
    energy falling on the bottom, effectively rendering the caustic intensity.

Dependencies:
    - numpy
    - matplotlib
    - scipy (optional, for Gaussian blur)
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Try to import scipy for high-quality blurring; fall back if missing
try:
    from scipy.ndimage import gaussian_filter
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


# ----------------------------
# Physics / Math Helpers
# ----------------------------

def sun_direction(elev_deg: float, az_deg: float) -> np.ndarray:
    """
    Calculate the unit direction vector of sunlight traveling *downward*.
    
    Args:
        elev_deg: Sun elevation above horizon (90 = zenith, 0 = horizon).
        az_deg:   Sun azimuth in x-y plane (0 along +x, 90 along +y).
        
    Returns:
        3-element numpy array representing the normalized direction vector.
    """
    elev = np.deg2rad(elev_deg)
    az = np.deg2rad(az_deg)
    
    # Spherical to Cartesian conversion
    sx = np.cos(elev) * np.cos(az)
    sy = np.cos(elev) * np.sin(az)
    sz = -np.sin(elev)  # Negative z implies downward direction
    
    v = np.array([sx, sy, sz], dtype=np.float64)
    return v / np.linalg.norm(v)


def refract_air_to_water(I: np.ndarray,
                         nx: np.ndarray, ny: np.ndarray, nz: np.ndarray,
                         n_air: float = 1.0, n_water: float = 1.333) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Vectorized Snell's Law refraction.
    
    Calculates the transmitted ray direction for many surface normals at once.
    
    Args:
        I:  Incident vector (3,), pointing downward.
        nx, ny, nz: Arrays of surface normal components (pointing upward).
        n_air:   Refractive index of air (approx 1.0).
        n_water: Refractive index of water (approx 1.333).
        
    Returns:
        (tx, ty, tz, cosi):
            tx, ty, tz: Components of the transmitted (refracted) ray vector.
            cosi: Cosine of the incident angle.
    """
    sx, sy, sz = I
    
    # Dot product N . I
    # Since N points up and I points down, the dot product is negative.
    # We want cos(theta_i) to be positive.
    cosi = -(nx * sx + ny * sy + nz * sz)
    cosi = np.clip(cosi, 0.0, 1.0) # Numerical stability

    eta = n_air / n_water
    
    # Snell's law in vector form:
    # T = eta * I + (eta * cos(theta_i) - sqrt(1 - eta^2 * sin^2(theta_i))) * N
    # k represents the term under the square root
    k = 1.0 - eta * eta * (1.0 - cosi * cosi)
    
    # Total Internal Reflection check (k < 0). 
    # For air->water this shouldn't happen, but good for robustness.
    sqrtk = np.sqrt(np.maximum(k, 0.0))

    c2 = eta * cosi - sqrtk
    
    tx = eta * sx + c2 * nx
    ty = eta * sy + c2 * ny
    tz = eta * sz + c2 * nz
    
    return tx, ty, tz, cosi


def fresnel_transmission_schlick(cosi: np.ndarray, n1: float, n2: float) -> np.ndarray:
    """
    Approximates the fraction of light transmitted into the water using Schlick's approximation.
    
    Returns T = 1 - R (Transmittance = 1 - Reflectance).
    """
    # Reflection coefficient at normal incidence
    r0 = ((n1 - n2) / (n1 + n2)) ** 2
    
    # Schlick's approximation for R
    R = r0 + (1.0 - r0) * (1.0 - cosi) ** 5
    
    return 1.0 - R


# ----------------------------
# Spectral Ripple Surface
# ----------------------------

class SpectralRipples:
    """
    Generates a band-limited height field using FFT-based synthesis.
    
    This follows the Tessendorf approach often used in oceanography graphics:
    H(k,t) = h0(k) * exp(i*w*t) + h0*(-k) * exp(-i*w*t)
    
    This method guarantees a real-valued spatial field that tiles seamlessly.
    """

    def __init__(self,
                 Nx: int = 256, Ny: int = 256,
                 Lx: float = 2.0, Ly: float = 2.0,
                 depth: float = 1.2,
                 target_rms: float = 0.004,
                 lambda0: float = 0.15, bandwidth: float = 0.07,
                 seed: int = 0):
        """
        Initialize the wave spectrum.
        
        Args:
            Nx, Ny: Grid resolution.
            Lx, Ly: Physical dimensions of the pool surface [meters].
            depth:  Water depth [meters]. Used for dispersion relation.
            target_rms: Desired root-mean-square wave height [meters].
            lambda0: Dominant wavelength of the ripples.
            bandwidth: Width of the frequency band (controls how chaotic the waves look).
        """
        self.Nx, self.Ny = Nx, Ny
        self.Lx, self.Ly = Lx, Ly
        self.depth = depth
        self.dx = Lx / Nx
        self.dy = Ly / Ny
        self.g = 9.81 # Gravity

        rng = np.random.default_rng(seed)

        # 1. Generate Wavevectors (k) on the FFT grid
        kx = 2.0 * np.pi * np.fft.fftfreq(Nx, d=self.dx)
        ky = 2.0 * np.pi * np.fft.fftfreq(Ny, d=self.dy)
        KX, KY = np.meshgrid(kx, ky)
        k = np.sqrt(KX**2 + KY**2)

        # 2. Dispersion Relation
        # Shallow/Intermediate water dispersion: w^2 = g * k * tanh(k * depth)
        # This controls the speed of waves based on their size.
        omega = np.sqrt(self.g * k * np.tanh(k * depth))
        omega[k == 0.0] = 0.0 # Handle DC component
        self.omega = omega

        # 3. Define the Spectrum S(k)
        # We use a Gaussian band-pass centered at k0 corresponding to lambda0.
        k0 = 2.0 * np.pi / lambda0
        sigma_k = 2.0 * np.pi / max(bandwidth, 1e-6)
        
        S = np.exp(-0.5 * ((k - k0) / sigma_k)**2)
        S[k == 0.0] = 0.0

        # 4. Generate Initial Amplitudes h0(k)
        # Random phases with magnitude derived from the spectrum S(k)
        h0 = (rng.standard_normal((Ny, Nx)) + 1j * rng.standard_normal((Ny, Nx))) * np.sqrt(S / 2.0)

        # 5. Enforce Hermitian Conjugate Symmetry
        # This ensures that when we take the inverse FFT, the result is purely real.
        # We need h0(-k) to construct the conjugate part.
        iy_neg = (-np.arange(Ny)) % Ny
        ix_neg = (-np.arange(Nx)) % Nx
        h0_neg = h0[iy_neg[:, None], ix_neg[None, :]]
        h0_star_neg = np.conj(h0_neg)

        # 6. Normalize Amplitudes to match target RMS height
        H_hat_0 = h0 + h0_star_neg
        h_spatial_0 = np.fft.ifft2(H_hat_0).real
        rms = h_spatial_0.std()
        scale = (target_rms / rms) if rms > 0 else 1.0
        
        self.h0 = h0 * scale
        self.h0_star_neg = h0_star_neg * scale

        # Precompute spatial grid for ray tracing later
        x = np.linspace(0.0, Lx, Nx, endpoint=False)
        y = np.linspace(0.0, Ly, Ny, endpoint=False)
        self.X, self.Y = np.meshgrid(x, y)

    def height(self, t: float) -> np.ndarray:
        """
        Calculate the water surface height field z = h(x,y) at time t.
        Uses the analytic time evolution in Fourier domain.
        """
        # Evolve phases: exp(i * w * t)
        exp_pos = np.exp(1j * self.omega * t)
        exp_neg = np.conj(exp_pos)
        
        # Combine forward and backward traveling waves
        H_hat = self.h0 * exp_pos + self.h0_star_neg * exp_neg
        
        # Convert back to spatial domain
        return np.fft.ifft2(H_hat).real


# ----------------------------
# Caustic Rendering Pipeline
# ----------------------------

def gaussian_blur_fft(img: np.ndarray, sigma_px: float) -> np.ndarray:
    """Fallback periodic Gaussian blur using FFT (if scipy is missing)."""
    H, W = img.shape
    ky = np.fft.fftfreq(H) * 2.0 * np.pi
    kx = np.fft.fftfreq(W) * 2.0 * np.pi
    KX, KY = np.meshgrid(kx, ky)
    # Gaussian filter in frequency domain
    G = np.exp(-0.5 * (sigma_px**2) * (KX*KX + KY*KY))
    return np.fft.ifft2(np.fft.fft2(img) * G).real


def caustics_frame(surface: SpectralRipples,
                   t: float,
                   bottom_res: int = 700,
                   sun_elev_deg: float = 60.0,
                   sun_az_deg: float = 20.0,
                   n_air: float = 1.0,
                   n_water: float = 1.333,
                   exposure: float = 60.0,
                   gamma: float = 0.7,
                   wrap: bool = True,
                   blur_sigma_px: float | None = None) -> np.ndarray:
    """
    Compute a single frame of the caustic pattern.
    
    1. Generate surface H(x,y,t).
    2. Compute gradients (slope) to get normals N.
    3. Refract sun rays I through N to get transmitted rays T.
    4. Intersect T with the bottom plane (z = -depth).
    5. Accumulate ray hits into a histogram (the caustic image).
    """
    # 1. Get Surface
    h = surface.height(t)
    X, Y = surface.X, surface.Y
    dx, dy = surface.dx, surface.dy
    Lx, Ly = surface.Lx, surface.Ly
    H = surface.depth

    # 2. Calculate Normals
    # Use finite differences for gradient
    hy, hx = np.gradient(h, dy, dx, edge_order=2)
    
    # Normalize to get unit vectors (-dh/dx, -dh/dy, 1) / magnitude
    denom = np.sqrt(1.0 + hx*hx + hy*hy)
    nx = -hx / denom
    ny = -hy / denom
    nz =  1.0 / denom

    I = sun_direction(sun_elev_deg, sun_az_deg)

    # 3. Refract Rays
    tx, ty, tz, cosi = refract_air_to_water(I, nx, ny, nz, n_air=n_air, n_water=n_water)

    # Filter out rays that reflect upwards (tz >= 0) or total internal reflection
    mask = tz < -1e-8
    if not np.any(mask):
        return np.zeros((bottom_res, bottom_res), dtype=np.float64)

    # 4. Intersect with Bottom
    # Plane equation: z = -H
    # Ray equation: P(s) = P_surf + s * T
    # Solve for s where P(s).z = -H => -H = h + s * tz  =>  s = (-H - h) / tz
    s = (-H - h) / tz
    bx = X + s * tx
    by = Y + s * ty

    # Wrap coordinates for seamless tiling (periodic boundary conditions)
    if wrap:
        bx = np.mod(bx, Lx)
        by = np.mod(by, Ly)

    # 5. Ray Intensity Weighting
    # T_fresnel: Amount of light transmitting through surface
    # cosi: Lambertian term (projected area)
    # 1/nz: Correction for converting flux per surface area to flux per horizontal area
    T = fresnel_transmission_schlick(cosi, n_air, n_water)
    w = T * cosi / np.maximum(nz, 1e-8)

    # Flatten arrays for histogram
    bx1 = bx[mask].ravel()
    by1 = by[mask].ravel()
    w1  = w[mask].ravel()

    # 6. Accumulate into 2D Histogram (The "Render")
    img, _, _ = np.histogram2d(
        by1, bx1,
        bins=(bottom_res, bottom_res),
        range=[[0.0, Ly], [0.0, Lx]],
        weights=w1
    )

    # 7. Post-Processing
    # Blur simulates the finite angular size of the sun (soft shadows)
    if blur_sigma_px is not None and blur_sigma_px > 0:
        if SCIPY_AVAILABLE:
            img = gaussian_filter(img, blur_sigma_px, mode="wrap" if wrap else "nearest")
        else:
            img = gaussian_blur_fft(img, blur_sigma_px)

    # Tonemapping (Compress dynamic range)
    # Normalize -> Log -> Gamma
    img = img / (img.mean() + 1e-12)
    img = np.log1p(exposure * img)
    img = img / (img.max() + 1e-12)
    img = img ** gamma
    
    return img


# ----------------------------
# Main Execution
# ----------------------------

def main():
    # --- Configuration ---
    Nx = Ny = 256          # Simulation grid resolution
    Lx = Ly = 2.0          # Simulation physical size [m]
    depth = 1.2            # Pool depth [m]
    bottom_res = 700       # Resolution of the output image
    
    # Wave parameters
    wave_rms = 0.004       # Wave height RMS [m] (~4mm)
    lambda0 = 0.14         # Main wavelength [m]
    bandwidth = 0.06       # Frequency spread
    seed = 3

    # Lighting
    sun_elev_deg = 62.0
    sun_az_deg = 25.0
    
    # Rendering
    exposure = 70.0
    gamma = 0.75
    blur_sigma_px = 1.0
    
    # Animation
    fps = 30
    frames = 300
    dt = 1.0 / fps

    print("Initializing Spectral Ripples...")
    surface = SpectralRipples(
        Nx=Nx, Ny=Ny, Lx=Lx, Ly=Ly, depth=depth,
        target_rms=wave_rms,
        lambda0=lambda0, bandwidth=bandwidth,
        seed=seed
    )

    # Setup Plot
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_title("Pool Caustics Simulation")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    
    # Initial Frame
    print("Rendering initial frame...")
    img0 = np.zeros((bottom_res, bottom_res), dtype=np.float64)
    im = ax.imshow(
        img0,
        origin="lower",
        extent=(0.0, Lx, 0.0, Ly),
        interpolation="nearest",
        cmap="gray",
        vmin=0.0, vmax=1.0
    )

    def update(k: int):
        t = k * dt
        img = caustics_frame(
            surface, t,
            bottom_res=bottom_res,
            sun_elev_deg=sun_elev_deg, sun_az_deg=sun_az_deg,
            exposure=exposure, gamma=gamma,
            blur_sigma_px=blur_sigma_px
        )
        im.set_data(img)
        ax.set_title(f"Pool Caustics | t = {t:6.3f} s")
        return (im,)

    print("Starting Animation...")
    ani = FuncAnimation(fig, update, frames=frames, interval=1000/fps, blit=True)
    plt.show()


if __name__ == "__main__":
    main()
