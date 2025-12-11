(* ::Package:: *)

(* ::Title:: *)
(*Pool Caustics Simulation*)


(* ::Subtitle:: *)
(*Geometric Optics Simulation of Light Patterns on Swimming Pool Bottoms*)


(* ::Text:: *)
(*This notebook implements a physically-based simulation of caustic patterns - *)
(*the beautiful "dancing light" formations seen at the bottom of swimming pools.*)
(**)
(*Author: Translated from Python implementation*)
(*Date: December 2025*)
(**)
(*The simulation combines:*)
(*  - Spectral synthesis (FFT-based Tessendorf waves) for realistic water surface*)
(*  - Geometric optics (Snell's Law) for light refraction*)
(*  - Ray tracing and accumulation for caustic intensity computation*)


(* ::Section:: *)
(*1. Physical Constants and Helper Functions*)


(* ::Text:: *)
(*We begin by defining the core physics functions that govern light behavior.*)
(*These functions are vectorized to work efficiently on 2D arrays representing*)
(*the water surface.*)


(* ::Subsection:: *)
(*1.1 Sun Direction Vector*)


(* ::Text:: *)
(*The sun direction is specified in spherical coordinates:*)
(*  - Elevation: angle above the horizon (0\[Degree] = horizon, 90\[Degree] = directly overhead)*)
(*  - Azimuth: compass direction in the horizontal plane*)
(**)
(*We convert this to a 3D Cartesian unit vector pointing DOWNWARD (from sun to surface).*)


sunDirection[elevDeg_, azDeg_] := Module[{elev, az, sx, sy, sz},
  elev = elevDeg * Degree;
  az = azDeg * Degree;
  sx = Cos[elev] * Cos[az];
  sy = Cos[elev] * Sin[az];
  sz = -Sin[elev];  (* Negative because rays travel downward *)
  Normalize[{sx, sy, sz}]
];


(* ::Subsection:: *)
(*1.2 Surface Gradient Computation*)


(* ::Text:: *)
(*To calculate how light refracts at the water surface, we need the surface normals,*)
(*which are derived from the height field gradients. We use central finite differences*)
(*with periodic boundary conditions to ensure seamless tiling.*)
(**)
(*Mathematical formulation:*)
(*  \[PartialD]h/\[PartialD]x \[TildeTilde] (h[i+1,j] - h[i-1,j]) / (2\[CapitalDelta]x)*)
(*  \[PartialD]h/\[PartialD]y \[TildeTilde] (h[i,j+1] - h[i,j-1]) / (2\[CapitalDelta]y)*)


gradient2D[h_, dx_, dy_] := Module[{ny, nx, hPadX, hPadY, gradX, gradY},
  {ny, nx} = Dimensions[h];

  (* ArrayPad with "Periodic" ensures seamless wrapping at domain boundaries *)
  hPadX = ArrayPad[h, {{0, 0}, {1, 1}}, "Periodic"];  (* Pad in x-direction *)
  hPadY = ArrayPad[h, {{1, 1}, {0, 0}}, "Periodic"];  (* Pad in y-direction *)

  (* Compute centered finite differences *)
  gradX = Table[
    (hPadX[[j, i+1]] - hPadX[[j, i-1]])/(2*dx),
    {j, ny}, {i, 2, nx+1}
  ];
  gradY = Table[
    (hPadY[[j+1, i]] - hPadY[[j-1, i]])/(2*dy),
    {j, 2, ny+1}, {i, nx}
  ];

  {gradX, gradY}
];


(* ::Section:: *)
(*2. Spectral Wave Surface Generation*)


(* ::Text:: *)
(*The water surface is generated using Fourier synthesis, following the method*)
(*popularized by Tessendorf (2001) for ocean rendering. This approach is preferred*)
(*over PDE-based simulations because it:*)
(*  1. Produces more aesthetically realistic pool-scale ripples*)
(*  2. Allows precise control over wavelength spectrum*)
(*  3. Is computationally efficient via FFT*)
(*  4. Guarantees periodic boundaries (seamless tiling)*)


(* ::Subsection:: *)
(*2.1 Physics Background: Dispersion Relation*)


(* ::Text:: *)
(*Water waves exhibit DISPERSION: their phase velocity depends on wavelength.*)
(*For water of finite depth H, the dispersion relation is:*)
(**)
(*  \[Omega]\.b2 = g\[CenterDot]k\[CenterDot]tanh(k\[CenterDot]H)*)
(**)
(*where:*)
(*  \[Omega] = angular frequency (rad/s)*)
(*  g = gravitational acceleration (9.81 m/s\.b2)*)
(*  k = wavenumber = 2\[Pi]/\[Lambda] (rad/m)*)
(*  H = water depth (m)*)
(**)
(*This reduces to:*)
(*  - Deep water (kH >> 1): \[Omega]\.b2 \[TildeTilde] g\[CenterDot]k  (short wavelengths)*)
(*  - Shallow water (kH << 1): \[Omega]\.b2 \[TildeTilde] g\[CenterDot]k\.b2\[CenterDot]H  (long wavelengths)*)


(* ::Subsection:: *)
(*2.2 Create Spectral Surface Function*)


createSurface[Nx_, Ny_, Lx_, Ly_, depth_, rms_, lambda0_, seed_] := Module[
  {dx, dy, g, kx, ky, KX, KY, k, omega, k0, sigma, S, h0, h0conj,
   scale, h0test, x, y, X, Y},

  (* Grid spacing and physical constants *)
  dx = Lx/Nx;
  dy = Ly/Ny;
  g = 9.81;

  SeedRandom[seed];

  (* ===== FREQUENCY DOMAIN SETUP ===== *)

  (* Create wavenumber grids corresponding to FFT frequency bins *)
  (* Note: FFT uses specific ordering: [0, 1, ..., N/2-1, -N/2, ..., -1] *)
  kx = 2*Pi*Join[Range[0, Floor[Nx/2]-1], Range[-Ceiling[Nx/2], -1]]/Lx;
  ky = 2*Pi*Join[Range[0, Floor[Ny/2]-1], Range[-Ceiling[Ny/2], -1]]/Ly;

  (* Create 2D wavenumber grids *)
  KX = ConstantArray[kx, Ny];
  KY = Transpose[ConstantArray[ky, Nx]];
  k = Sqrt[KX^2 + KY^2];

  (* ===== DISPERSION RELATION ===== *)

  (* Calculate angular frequency from wavenumber using shallow/intermediate water formula *)
  omega = Sqrt[g * k * Tanh[k * depth]];

  (* Handle singularities at k=0 *)
  omega = omega /. {Indeterminate -> 0, ComplexInfinity -> 0};

  (* ===== WAVE SPECTRUM S(k) ===== *)

  (* We use a Gaussian spectrum centered at k\:2080 corresponding to dominant wavelength \[Lambda]\:2080 *)
  k0 = 2*Pi/lambda0;
  sigma = 0.16 * k0;  (* Spectral width: 16% of central wavenumber *)

  (* Narrow-band Gaussian around k\:2080 *)
  S = Exp[-0.5*((k - k0)/sigma)^2];

  (* Add high-k suppression to prevent aliasing *)
  S = S * Exp[-(k/(3*k0))^4];

  (* Zero out any problematic values *)
  S = S /. {Indeterminate -> 0};

  (* ===== RANDOM FOURIER AMPLITUDES ===== *)

  (* Generate complex Gaussian random amplitudes h\:2080(k) *)
  (* The variance at each k is proportional to the spectrum S(k) *)
  h0 = (RandomVariate[NormalDistribution[], {Ny, Nx}] +
        I*RandomVariate[NormalDistribution[], {Ny, Nx}]) * Sqrt[S/2];

  (* ===== HERMITIAN SYMMETRY ===== *)

  (* For the height field to be real-valued, we need h\:2080(-k) = h\:2080*(k) *)
  (* This ensures that h(x,y,t) = IFFT[H(k,t)] is purely real *)
  h0conj = Conjugate[h0[[
    Table[Mod[-i, Ny]+1, {i, 0, Ny-1}],
    Table[Mod[-i, Nx]+1, {i, 0, Nx-1}]
  ]]];

  (* ===== AMPLITUDE SCALING ===== *)

  (* Scale amplitudes to achieve target RMS wave height *)
  h0test = Re[InverseFourier[h0 + h0conj, FourierParameters -> {1, -1}]];
  scale = rms / StandardDeviation[Flatten[h0test]];
  h0 = h0 * scale;
  h0conj = h0conj * scale;

  (* ===== SPATIAL GRID ===== *)

  (* Create coordinate meshgrid for ray tracing *)
  x = Range[0, Lx-dx, dx];
  y = Range[0, Ly-dy, dy];
  X = ConstantArray[x, Ny];
  Y = Transpose[ConstantArray[y, Nx]];

  (* Return surface data as an Association (like Python dictionary) *)
  <|"omega" -> omega, "h0" -> h0, "h0conj" -> h0conj,
    "X" -> X, "Y" -> Y, "Lx" -> Lx, "Ly" -> Ly,
    "dx" -> dx, "dy" -> dy, "depth" -> depth|>
];


(* ::Subsection:: *)
(*2.3 Time Evolution of Surface*)


(* ::Text:: *)
(*The height field at time t is computed via inverse FFT of the time-evolved spectrum:*)
(**)
(*  H(k,t) = h\:2080(k)\[CenterDot]exp(i\[Omega]t) + h\:2080*(-k)\[CenterDot]exp(-i\[Omega]t)*)
(*  h(x,y,t) = IFFT[H(k,t)]*)
(**)
(*This analytical time evolution in Fourier space ensures the field remains real-valued*)
(*and obeys the correct dispersion relation.*)


getHeight[surf_, t_] := Re[InverseFourier[
  surf["h0"] * Exp[I*surf["omega"]*t] +
  surf["h0conj"] * Exp[-I*surf["omega"]*t],
  FourierParameters -> {1, -1}
]];


(* ::Section:: *)
(*3. Caustics Rendering via Ray Tracing*)


(* ::Text:: *)
(*This is the core rendering function that simulates light refraction through*)
(*the water surface and accumulates ray hits on the pool bottom.*)
(**)
(*Physical process:*)
(*  1. Sun rays approach the surface from direction I (incident)*)
(*  2. Rays refract at the air-water interface following Snell's Law*)
(*  3. Refracted rays propagate through water until hitting the bottom*)
(*  4. Each ray contributes intensity weighted by Fresnel transmission*)
(*  5. Ray hits are accumulated into a 2D histogram (the caustic image)*)


(* ::Subsection:: *)
(*3.1 Main Rendering Function*)


renderFrame[surf_, t_, res_] := Module[
  {h, X, Y, dx, dy, Lx, Ly, depth, gradX, gradY, denom, nx, ny, nz,
   sunDir, sx, sy, sz, cosi, eta, k, sqrtk, c2, tx, ty, tz,
   s, bx, by, r0, R, T, w, good, xv, yv, wv, ix, iy, img, mean, maxval},

  (* ===== GET SURFACE GEOMETRY ===== *)

  h = getHeight[surf, t];
  X = surf["X"]; Y = surf["Y"];
  dx = surf["dx"]; dy = surf["dy"];
  Lx = surf["Lx"]; Ly = surf["Ly"];
  depth = surf["depth"];

  (* ===== COMPUTE SURFACE NORMALS ===== *)

  (* Get surface gradients *)
  {gradX, gradY} = gradient2D[h, dx, dy];

  (* Surface normal N = (-\[PartialD]h/\[PartialD]x, -\[PartialD]h/\[PartialD]y, 1) / ||\[CenterDot]|| (pointing upward) *)
  denom = Sqrt[1 + gradX^2 + gradY^2];
  nx = -gradX/denom;
  ny = -gradY/denom;
  nz = 1/denom;

  (* ===== SNELL'S LAW REFRACTION ===== *)

  (* Sun direction vector (incident ray) *)
  sunDir = sunDirection[62, 25];  (* 62\[Degree] elevation, 25\[Degree] azimuth *)
  {sx, sy, sz} = sunDir;

  (* Cosine of incident angle: cos(\[Theta]\:1d62) = -N\[CenterDot]I *)
  (* (Negative because N points up, I points down) *)
  cosi = Clip[-(nx*sx + ny*sy + nz*sz), {0, 1}];

  (* Relative refractive index: \[Eta] = n_air / n_water \[TildeTilde] 1/1.333 *)
  eta = 1/1.333;

  (* Vector form of Snell's Law: *)
  (* T = \[Eta]\[CenterDot]I + (\[Eta]\[CenterDot]cos(\[Theta]\:1d62) - \[Sqrt](1 - \[Eta]\.b2\[CenterDot]sin\.b2(\[Theta]\:1d62)))\[CenterDot]N *)
  k = 1 - eta^2*(1 - cosi^2);
  sqrtk = Sqrt[Clip[k, {0, Infinity}]];  (* Handle total internal reflection *)
  c2 = eta*cosi - sqrtk;

  (* Transmitted (refracted) ray direction *)
  tx = eta*sx + c2*nx;
  ty = eta*sy + c2*ny;
  tz = eta*sz + c2*nz;

  (* ===== RAY-PLANE INTERSECTION ===== *)

  (* Intersect refracted rays with pool bottom at z = -depth *)
  (* Ray equation: P(s) = (X,Y,h) + s\[CenterDot](tx,ty,tz) *)
  (* Bottom plane: z = -depth *)
  (* Solve: h + s\[CenterDot]tz = -depth  =>  s = (-depth - h)/tz *)
  s = (-depth - h)/(tz - 1*^-20);  (* Small offset prevents division by zero *)

  (* Bottom hit coordinates (with periodic wrapping) *)
  bx = Mod[X + s*tx, Lx];
  by = Mod[Y + s*ty, Ly];

  (* ===== FRESNEL TRANSMISSION WEIGHTING ===== *)

  (* Schlick's approximation for Fresnel reflectance *)
  r0 = ((1 - 1.333)/(1 + 1.333))^2;
  R = r0 + (1 - r0)*(1 - cosi)^5;
  T = 1 - R;  (* Transmission coefficient *)

  (* Ray weight combines: *)
  (*   - Fresnel transmission T *)
  (*   - Projected solid angle cos(\[Theta]\:1d62) *)
  (*   - Surface area correction 1/nz *)
  w = T * cosi / Clip[nz, {1*^-8, Infinity}];

  (* ===== FILTER VALID RAYS ===== *)

  (* Only keep rays traveling downward in water (tz < 0) *)
  good = Flatten[Position[Flatten[tz], _?((# < 0)&)]];
  xv = Flatten[bx][[good]];
  yv = Flatten[by][[good]];
  wv = Flatten[w][[good]];

  (* ===== ACCUMULATE INTO IMAGE ===== *)

  (* Create caustic intensity map via weighted histogram *)
  img = ConstantArray[0.0, {res, res}];

  (* Convert physical coordinates to pixel indices *)
  ix = Clip[Floor[(xv/Lx)*res] + 1, {1, res}];
  iy = Clip[Floor[(yv/Ly)*res] + 1, {1, res}];

  (* Accumulate ray contributions *)
  Do[img[[iy[[i]], ix[[i]]]] += wv[[i]], {i, Length[wv]}];

  (* ===== POST-PROCESSING ===== *)

  (* Simple 5-point blur to simulate finite sun angular size *)
  img = (img +
         RotateRight[img, {0, 1}] + RotateLeft[img, {0, 1}] +
         RotateRight[img, {1, 0}] + RotateLeft[img, {1, 0}])/5.0;

  (* ===== TONEMAPPING ===== *)

  (* Map high dynamic range to displayable [0,1] range *)
  (* 1. Normalize by mean (exposure compensation) *)
  mean = Mean[Flatten[img]];
  img = img / (mean + 1*^-12);

  (* 2. Logarithmic compression (like human eye response) *)
  img = Log[1 + 70*img];  (* 70 is exposure parameter *)

  (* 3. Normalize to [0,1] *)
  maxval = Max[Flatten[img]];
  img = img / (maxval + 1*^-12);

  (* 4. Gamma correction for display *)
  img^0.8
];


(* ::Section:: *)
(*4. Visualization and Color Mapping*)


(* ::Text:: *)
(*We use a custom color map that mimics the appearance of water:*)
(*  - Dark blue in shadows (low intensity)*)
(*  - Turquoise in mid-tones*)
(*  - Bright cyan to white in caustic highlights*)


waterColor = Blend[{
  {0, RGBColor[0.01, 0.05, 0.12]},    (* Very dark blue - deep shadows *)
  {0.25, RGBColor[0.02, 0.2, 0.4]},   (* Deep blue *)
  {0.5, RGBColor[0, 0.45, 0.7]},      (* Turquoise - medium intensity *)
  {0.75, RGBColor[0.6, 0.9, 0.95]},   (* Pale cyan - bright areas *)
  {1, RGBColor[1, 1, 1]}              (* White - caustic peaks *)
}, #]&;


(* ::Section:: *)
(*5. Simulation Parameters*)


(* ::Text:: *)
(*Physical parameters chosen to match realistic swimming pool conditions:*)


(* ::Subsection:: *)
(*5.1 Geometry*)


(* ::Item:: *)
(*Domain: 2m * 2m pool section*)


(* ::Item:: *)
(*Grid: 128 * 128 points (sufficient for smooth caustics)*)


(* ::Item:: *)
(*Pool depth: 1.2 m (typical shallow end)*)


(* ::Subsection:: *)
(*5.2 Wave Characteristics*)


(* ::Item:: *)
(*Dominant wavelength \[Lambda]\:2080 = 0.45 m (45 cm) - realistic pool ripples*)


(* ::Item:: *)
(*Wave height RMS = 8 mm - gentle wind-driven ripples*)


(* ::Item:: *)
(*Spectral width = 16% - natural wave variability*)


(* ::Subsection:: *)
(*5.3 Lighting*)


(* ::Item:: *)
(*Sun elevation: 62\[Degree] (mid-afternoon)*)


(* ::Item:: *)
(*Sun azimuth: 25\[Degree] (sun from slightly east of north)*)


(* ::Item:: *)
(*Refractive indices: n_air = 1.0, n_water = 1.333*)


(* ::Section:: *)
(*6. Run Simulation*)


(* ::Text:: *)
(*Now we execute the simulation in three steps:*)
(*  1. Create the spectral wave surface*)
(*  2. Render a test frame to verify correctness*)
(*  3. Generate full animation sequence*)


Print["="*70];
Print["POOL CAUSTICS SIMULATION"];
Print["="*70];
Print[""];
Print["Creating spectral wave surface..."];
Print["  Domain: 2m \[Times] 2m"];
Print["  Grid: 128 \[Times] 128 points"];
Print["  Wavelength \[Lambda]\:2080 = 0.45 m (45 cm)"];
Print["  Wave RMS height = 8 mm"];
Print["  Pool depth = 1.2 m"];
Print[""];

surface = createSurface[128, 128, 2.0, 2.0, 1.2, 0.008, 0.45, 3];

Print["Surface created successfully!"];
Print[""];
Print["Rendering test frame at t=0.5s..."];

test = renderFrame[surface, 0.5, 400];

Print[""];
Print["Displaying test frame...");
Print["You should see caustic patterns with NO vertical artifacts.");
Print[""];

ArrayPlot[test,
  ColorFunction -> waterColor,
  ColorFunctionScaling -> False,
  PlotLabel -> "Pool Caustics - Test Frame (t = 0.5 s)",
  FrameLabel -> {"x position [m]", "y position [m]"},
  ImageSize -> 600,
  Frame -> True
]


(* ::Section:: *)
(*7. Generate Animation*)


(* ::Text:: *)
(*Generate 60 frames spanning 2 seconds (30 fps) to visualize the time evolution*)
(*of caustic patterns as waves propagate across the surface.*)


Print[""];
Print["="*70];
Print["GENERATING ANIMATION"];
Print["="*70];
Print[""];
Print["Rendering 60 frames (this will take 2-3 minutes)..."];

frames = Table[
  If[Mod[k, 10] == 0,
    Print["  Progress: Frame ", k, "/60 (", Round[100*k/60], "%)"]
  ];
  renderFrame[surface, k/30.0, 400],
  {k, 0, 59}
];

Print[""];
Print["Rendering complete! Creating animation..."];
Print[""];

ListAnimate[
  Table[
    ArrayPlot[frames[[k]],
      ColorFunction -> waterColor,
      ColorFunctionScaling -> False,
      PlotLabel -> Row[{
        "Pool Caustics  |  ",
        "t = ", NumberForm[k/30.0, {4, 2}], " s  |  ",
        "\[Lambda]\:2080 = 0.45 m"
      }],
      FrameLabel -> {"x [m]", "y [m]"},
      ImageSize -> 600,
      Frame -> True
    ],
    {k, 1, 60}
  ],
  AnimationRunning -> True,
  AnimationRate -> 30
]


(* ::Section:: *)
(*8. Export Options*)


(* ::Text:: *)
(*To export the animation as a GIF or video, uncomment and run the code below:*)


(* ::Input:: *)
(*Export["~/Desktop/pool_caustics.gif",*)
(*  Table[Image[frames[[k]], "Real"], {k, 1, Length[frames]}],*)
(*  "AnimationRepetitions" -> Infinity,*)
(*  "DisplayDurations" -> 1/30*)
(*];*)


(* ::Section:: *)
(*9. Scientific Notes*)


(* ::Text:: *)
(*This simulation demonstrates several key physical phenomena:*)


(* ::Subsection:: *)
(*9.1 Why Spectral Synthesis?*)


(* ::Text:: *)
(*We chose FFT-based spectral synthesis over solving PDEs (like the shallow water*)
(*equations) because:*)
(*  - It produces more aesthetically realistic pool ripples*)
(*  - Allows precise control over the wavelength spectrum*)
(*  - Is computationally efficient via Fast Fourier Transform*)
(*  - Guarantees seamless periodic boundaries for animation loops*)


(* ::Subsection:: *)
(*9.2 Caustic Formation*)


(* ::Text:: *)
(*Caustics are NOT images of the surface ripples - they are concentrations of light*)
(*caused by wave focusing. The bright lines occur where:*)
(*  1. Surface curvature is locally concave (acts like a converging lens)*)
(*  2. Many rays are redirected to converge at the same bottom location*)
(*  3. The Jacobian determinant of the ray mapping approaches zero*)


(* ::Subsection:: *)
(*9.3 Key Approximations*)


(* ::Text:: *)
(*This simulation makes several standard approximations:*)
(*  - Geometric optics (ignores wave nature of light, valid for \[Lambda]_light << \[Lambda]_water)*)
(*  - Single-scattering (ignores multiple reflections between surface and bottom)*)
(*  - Unpolarized light (uses scalar Fresnel coefficients)*)
(*  - Flat bottom (no bottom topology)*)
(*  - Clear water (no absorption or scattering in volume)*)


(* ::Subsection:: *)
(*9.4 Comparison with Python Version*)


(* ::Text:: *)
(*This Wolfram Language implementation is a faithful translation of the Python code,*)
(*maintaining the same:*)
(*  - Physical models (dispersion, refraction, Fresnel)*)
(*  - Numerical methods (FFT synthesis, ray tracing)*)
(*  - Tonemapping and color grading*)
(**)
(*Key differences:*)
(*  - Wolfram uses 1-based indexing vs Python's 0-based*)
(*  - Fourier transforms use different conventions (handled via FourierParameters)*)
(*  - Array operations syntax differs but mathematics is identical*)


(* ::Section:: *)
(*End of Simulation*)


Print[""];
Print["="*70];
Print["SIMULATION COMPLETE"];
Print["="*70];
Print[""];
Print["The animation above shows caustic patterns evolving over 2 seconds."];
Print["Notice how the bright lines move and dance as waves propagate!");
Print[""];
Print["To run again with different parameters, edit Section 5 and re-evaluate."];
