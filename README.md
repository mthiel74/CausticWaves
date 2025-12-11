# Caustic Waves Simulator

A Python-based simulation of optical caustics on a pool bottom, powered by geometric optics and spectral surface synthesis.

![Caustics Simulation](Output/caustics_simulation.gif)

*Video generated using `CausticWaves original.ipynb`*

## Project Overview

This project simulates the beautiful "dancing light" patterns (caustics) seen at the bottom of swimming pools. It models the interaction of sunlight with a dynamic water surface to compute the intensity distribution on the pool floor.

### Core Technologies

*   **Spectral Synthesis (Tessendorf Waves)**: The water surface is generated using Fast Fourier Transforms (FFT) on a band-limited spectrum (Phillips/Tessendorf spectrum).
    *   *Note*: This approach ("Good Waves") was selected over PDE-based fluid simulations (e.g., Shallow Water Equations) because it produces far more realistic and aesthetically pleasing ripple patterns for optical rendering.
*   **Geometric Optics**: Light rays are traced from a directional light source (Sun), refracted at the air-water interface using **Snell's Law**, and intersected with the pool bottom.
*   **Accumulation Buffer**: The intensity on the bottom is computed by accumulating ray hits into a high-resolution grid (histogram), creating the characteristic sharp bright lines of caustics.

## Features

*   **Realistic Dispersion**: Simulates wave dispersion where phase velocity depends on wavelength ($\omega^2 = gk \tanh(kH)$).
*   **Fresnel Effects**: Approximates surface transmission using Schlick's approximation (rays are weighted by how much light actually enters the water).
*   **Tonemapping**: Applies logarithmic exposure and gamma correction to handle the high dynamic range of the focused light rays.
*   **Animation**: Real-time visualization using `matplotlib.animation`.

## Code Structure

The core logic is contained in the Jupyter Notebook `CausticWaves original.ipynb` and the standalone script `simulation.py`.
*   `SpectralRipples`: Class handling the FFT generation of the height field.
*   `refract_air_to_water`: Vectorized implementation of 3D Snell's refraction.
*   `caustics_frame`: Main rendering pipeline for a single time step.

## Usage

**Recommended: Start with the Jupyter Notebook**

1.  Open `CausticWaves original.ipynb` in Jupyter or VS Code.
2.  Run the cell to generate the simulation.
3.  The animation will render showing the evolving caustic patterns with blue color grading.
4.  Video output will be saved to `Output/pool_caustics_clean.mp4`.

**Alternative: Standalone Script**

Run `python simulation.py` for a command-line version with detailed comments and grayscale output.

## Dependencies

*   `numpy`
*   `matplotlib`
*   `scipy` (optional, for Gaussian blur post-processing)